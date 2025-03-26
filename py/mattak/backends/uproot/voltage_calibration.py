import mattak.Dataset

from typing import Union, Optional, List
import numpy
import uproot


class VoltageCalibration(object):

    def __init__(
            self, path : str,
            caching : bool = True, caching_mode : str = "lookup",
            fit_min : float = -1.3, fit_max : float = 0.7,
            table_step_size : float = 0):
        """ Helper class to apply the voltage calibration for the uproot backend.

        This class reads-in either the "calibration root files" which contain the parameter
        of the 9th order polynominals + residuals, or directly the entire bias scan.

        It caches calibration "tables" to increase performance. Two different caching methods
        are available: "lookup" and "interpolate".

        With "lookup" (default) we cache the voltage value (!) associate to each ADC bit [-2048, 2048).
        This requires the cached table to store 24 (channels) x 4096 (LAB4D samples) x 4096
        (ADC-bits) values. Those tables are created using linear interpolation in between measured bits.
        To keep memory consumption managable the voltage values are encoded as 16-bit unsigned integers.
        This will require ~ 0.8 Gb but will provide the best performance per event.

        With "interpolate" we only cache the ADC counts (!) which correspond to the sampled voltages^*.
        This results in the following table: 24 (channels) x 4096 (LAB4D samples) x N
        with N being the number of sampled voltages. From those tables we linearly interpolate the voltage
        matching the requested ADC bit.

        ^* in practice we can upsample the voltages to exploit the 9th-order polynominal instead of linear
        interpolation.

        Parameters
        ----------

        path : str
            Path to the calibration root file. You can either pass a bias scan file or a
            voltage calibration constant file.

        caching : bool (Default: True)
            If True, the adc tables for the calibration are cached.

        caching_mode : "lookup" (default) or "interpolate"
            See explanation above. TL;DR:

                * "lookup": fastest, memory intensive (~0.8Gb)
                * "interpolate": slower, less memory (~0.3Gb) (depending of `table_step_size`)

        fit_min : float (Default: -1.3 V)
            Lower bound for the upsample voltage range (in Volt).

        fit_max : float (Default: 0.7 V)
            Upper bound for the upsample voltage range (in Volt).

        table_step_size : float (Default: 0.001)
            Step size of the voltage used for the cached tables. If 0, do not upsample the measured voltages.
            Typical step size of measurements is 16 mV (155 entries)
        """

        self.NUM_CHANNELS = mattak.Dataset.AbstractDataset.NUM_CHANNELS
        self.NUM_WF_SAMPLES = mattak.Dataset.AbstractDataset.NUM_WF_SAMPLES
        self.NUM_DIGI_SAMPLES = mattak.Dataset.AbstractDataset.NUM_DIGI_SAMPLES
        self.NUM_BITS = 4096

        self.table_step_size = table_step_size
        self.__caching = caching
        self.__caching_mode = caching_mode

        self.fit_min = fit_min
        self.fit_max = fit_max

        self.__voltage = None
        if self.__caching:
            if self.__caching_mode not in ["lookup", "interpolate"]:
                raise ValueError("Argument `caching_mode` can either be \"lookup\" or \"interpolate\" but is: "
                                 f"{caching_mode}")

            self.__cached_tables = numpy.array([None] * self.NUM_CHANNELS, dtype=object)

        self.__path = path
        if path.endswith(".root"):
            with uproot.open(path) as cal_file:
                if "coeffs_tree" in cal_file:
                    self.full_bias_scan = False
                    self.__cal_param = unpack_cal_parameters(cal_file)
                    self.__cal_residuals_v, self.__cal_residuals_adc = unpack_cal_residuals(cal_file)


                    if numpy.any(numpy.diff(self.__cal_residuals_v, axis=0) != 0):
                        raise ValueError("The pedestal voltage of the bias scan (residual) is different for the two DAC, "
                                        "the code expects them to be the same!")

                    if numpy.any(self.__cal_residuals_v < self.fit_min) or numpy.any(self.__cal_residuals_v > self.fit_max):
                        raise ValueError("The pedestal voltage of the bias scan (residual) exceeds fit limits!")

                    self.fit_min = max([self.fit_min, self.__cal_residuals_v[0][0]])
                    self.fit_max = min([self.fit_max, self.__cal_residuals_v[0][-1]])

                    self.time = int(cal_file["general_tree"]["startTime"].array(library="np")[0])
                    self.station = int(cal_file["general_tree"]["stationNumber"].array(library="np")[0])

                elif "pedestals" in cal_file:
                    self.full_bias_scan = True
                    self.__vbias, self.__adc = unpack_raw_bias_scan(cal_file)

                    if numpy.any(self.__vbias[:, 0] != self.__vbias[:, 1]):
                        raise ValueError("The pedestal voltage of the bias scan is different for the two DAC, "
                                            "the code expects them to be the same!")

                    self.time = int(cal_file["pedestals"]["when"].array(library="np")[0])
                    self.station = int(cal_file["pedestals"]["station_number"].array(library="np")[0])

                    # No need to check if `self.__vbias`` is within `fit_min``, `fit_max`` here, that happens
                    # in `get_adcs_from_biasscan`.
                else:
                    raise ValueError("No 'coeffs_tree' or 'pedestals' keys found in the root file")
        else:
            raise ValueError(f"{path} is not recognized as a root file")


    def getStationNumber(self):
        """ To copy the pyroot interface """
        return self.station

    def getStartTime(self):
        """ To copy the pyroot interface """
        return self.time

    def __set_voltage(self, value):
        """ Set the voltage vector which correspond to the `adc_table` for calibration """
        if self.__voltage is None:
            self.__voltage = value
        else:
            if not numpy.all(self.__voltage == value):
                raise ValueError("The voltage vector for the adc table changed!")


    def get_lookup_table_per_adc(self, adcs, voltage):
        """ Interpolates the table ADC -> Voltage to every possible ADC count to create a lookup table

        Parameters
        ----------
        adcs : array
            ADC counts corresponding to bias voltage vector.
            Can be of shape (channel, samples, counts) or (samples, counts).

        voltage : array
            Bias voltage corresponding to ADC counts. Has shape counts.

        Returns
        -------
        voltage : array
            Lookup table for voltage values associated to each ADC bit [0 .. self.NUM_BITS).
            The voltage is encoded as 16-bit unsigned ints to save memory.
        """
        if adcs.shape[0] == self.NUM_DIGI_SAMPLES:
            adcs = numpy.expand_dims(adcs, axis=0)

        adcs = numpy.squeeze([[
            numpy.interp(
                numpy.arange(self.NUM_BITS) - self.NUM_BITS // 2,
                adcs_sample, voltage, left=self.fit_min, right=self.fit_max)
            for adcs_sample in adcs_channel] for adcs_channel in adcs])

        return self.convert_to_ints(adcs)


    def get_adcs_from_parameters(self, channel=None, add_residual=True):
        """ Create ADC count tables form the 9-pol parameter of the calibration  """

        assert not self.full_bias_scan, "No calibration available"

        if self.table_step_size == 0:
            voltage = self.__cal_residuals_v[0]
        else:
            voltage = numpy.linspace(
                self.fit_min, self.fit_max, int((self.fit_max - self.fit_min) // self.table_step_size))

        self.__set_voltage(voltage)

        if channel is None:
            parameters = self.__cal_param.reshape(self.NUM_CHANNELS, self.NUM_DIGI_SAMPLES, 10)
        else:
            parameters = numpy.expand_dims(
                self.__cal_param[self.NUM_DIGI_SAMPLES * channel:self.NUM_DIGI_SAMPLES * (channel + 1)], axis=0)

        adc_samples = numpy.array([
            [numpy.polyval(param_sample[::-1], voltage) for param_sample in param_channel]
            for param_channel in parameters])

        if add_residual:
            if self.table_step_size == 0:
                residuals = self.__cal_residuals_adc
            else:
                # we have to interpolate the residuals such that they correspond to the correct voltage
                # w.r.t adc = pol(voltage)
                residuals = [numpy.interp(voltage, self.__cal_residuals_v[i], self.__cal_residuals_adc[i]) for i in range(len(self.__cal_residuals_adc))]

            channels = range(self.NUM_CHANNELS) if channel is None else [channel]
            for idx, ch in enumerate(channels):
                # first residual for ch 0 .. 11, second for 12 .. 23
                adc_samples[idx] = adc_samples[idx] + residuals[idx]

        if self.__caching_mode == "lookup":
            full_lookup_tables = self.get_lookup_table_per_adc(adc_samples, voltage)
            return voltage, full_lookup_tables
        else:
            return voltage, numpy.squeeze(adc_samples)


    def convert_to_ints(self, adc_table):
        return numpy.asarray(
            (adc_table - self.fit_min) / (self.fit_max - self.fit_min) * (2 ** 16 - 1),
            dtype=numpy.uint16)


    def convert_to_floats(self, adc_table):
        return numpy.asarray(
            adc_table, dtype=float) / (2 ** 16 - 1) * (self.fit_max - self.fit_min) + self.fit_min


    def get_adcs_from_biasscan(self):
        """ Create ADC count tables directly from the bias scans  """

        assert self.full_bias_scan, "No bias scan available"

        vbias, adc = rescale_adc(self.__vbias, self.__adc)
        adc_cut = [[] for i in range(self.NUM_CHANNELS)]
        adc_cut[:12] = adc[:12, :, numpy.all([self.fit_min < vbias[:, 0], vbias[:, 0] < self.fit_max], axis=0)]
        adc_cut[12:] = adc[12:, :, numpy.all([self.fit_min < vbias[:, 1], vbias[:, 1] < self.fit_max], axis=0)]
        adc = numpy.array(adc_cut)
        vbias =  numpy.array([[v for v in vbias[:, DAC] if self.fit_min < v < self.fit_max] for DAC in range(2)])

        assert numpy.all(vbias[0] == vbias[1]), "Bias voltage of the two DAC is not equal"

        self.__set_voltage(vbias[0])
        self.fit_min = max([self.fit_min, vbias[0][0]])
        self.fit_max = min([self.fit_max, vbias[0][-1]])

        if self.__caching_mode == "lookup":
            # Running that per channel avoids peaks in the memory consumption and is also faster
            # Weirdly also the per event time to apply the calibration decreased ...
            full_lookup_tables = numpy.array([
                self.get_lookup_table_per_adc(adc_per_channel, vbias[0])
                for adc_per_channel in adc])
            return vbias[0], full_lookup_tables
        else:
            # There is no need to "upsample", i.e., linearly interpolate since we do this in calibrate() anyway
            return vbias[0], adc


    def __get_cached_table(self, channel, sample=None):
        """ Returns the cached table of ADC counts for a given voltage """
        if self.__cached_tables[channel] is None:

            if not self.full_bias_scan:
                _, adc = self.get_adcs_from_parameters(channel=channel)
                self.__cached_tables[channel] = adc
            else:
                _, adc = self.get_adcs_from_biasscan()
                self.__cached_tables = adc

        if sample is None:
            return self.__cached_tables[channel]
        else:
            return self.__cached_tables[channel][sample]


    def calibrate(
            self, waveform_array : numpy.ndarray, starting_window : Union[float, int],
            channels : Optional[List[int]] = None) -> numpy.ndarray:
        """
        The calibration function that transforms waveforms from ADC to voltage. Uses caching.

        Parameters
        ----------
        waveform_array : array of shape (num_channels, num_wf_samples)
            array of one waveform
        starting_window : int | float
            the sample on which the run started
        channels : list(int) (Default: None -> range(self.num_channels))
            List of channels. Length need to match with first dimension of waveform_array

        Returns
        -------
        waveform_volt : array of shape (num_channels, num_wf_samples)
            calibrated waveform in volt
        """

        channels = channels or list(range(self.NUM_CHANNELS))
        waveform_volt = numpy.zeros_like(waveform_array, dtype=float)

        assert numpy.iinfo(type(starting_window[0])).max >= 4096, f"Data type of starting window has not enought precission:\n({numpy.iinfo(type(starting_window[0]))})"

        for idx, (ch, starting_window_channel, wf_channel) in enumerate(zip(channels, starting_window, waveform_array)):
            # residuals split over DACs
            samples_idx = (128 * starting_window_channel + numpy.arange(self.NUM_WF_SAMPLES)) % self.NUM_WF_SAMPLES
            if starting_window_channel >= 16:
                samples_idx += self.NUM_WF_SAMPLES

            if self.__caching_mode == "lookup":
                # Accuracy is "exact" -> we have a voltage value for every possible ADC count (-2048 .. 2047)
                voltage_lookup_table = self.__get_cached_table(ch)
                trace = voltage_lookup_table[samples_idx, wf_channel + self.NUM_BITS // 2]
                waveform_volt[idx] = self.convert_to_floats(trace)
            else:
                # Accuracy is not "exact" -> we have to interpolate between the values in the lookup table
                for sample_wf, (sample_lab, adc) in enumerate(zip(samples_idx, wf_channel)):
                    adc_samples = self.__get_cached_table(ch, sample_lab)
                    volt = numpy.interp(adc, adc_samples, self.__voltage, left=self.fit_min, right=self.fit_max)
                    waveform_volt[idx, sample_wf] = volt

        return waveform_volt


    def __call__(self, waveform_array : numpy.ndarray, starting_window : Union[float, int]) -> numpy.ndarray:
        if self.__caching:
            return self.calibrate(waveform_array, starting_window)
        else:
            return calibrate(
                waveform_array, self.__cal_param, self.__cal_residuals_v,
                self.__cal_residuals_adc, starting_window, upsampling=False)  # already upsampled

    def plot_ch(self, ax1=None, xs=numpy.linspace(-1000, 1000, 100), ch=0):

        from matplotlib import pyplot as plt

        if ax1 is None:
            fig, (ax, ax2) = plt.subplots(2, 1, height_ratios=[3, 1], sharex=True)
        else:
            ax2 = ax1

        # tables = numpy.array([self.__get_adc_table(ch, i) for i in range(self.NUM_DIGI_SAMPLES)])
        out = []
        out2 = []
        for x in xs:
            # vs = [numpy.interp(x, t, self.__voltage) * 1000 for t in tables]
            # out.append([numpy.mean(vs), numpy.std(vs)])
            vs2 = self.calibrate([numpy.full(self.NUM_WF_SAMPLES, int(x))], [0], channels=[ch]) * 1000
            out2.append([numpy.mean(vs2), numpy.std(vs2)])

        out = numpy.array(out2)

        if ax1 is None:
            ax.errorbar(xs, out[:, 0], out[:, 1], ls="",  marker="s", markersize=0.5, label=f"Ch {ch} (mean over all samples)")
            ax.plot(xs, xs * 2500 / 4096, "k--", lw=1, label="2500 / 4096")
            ax.set_ylabel("output voltage / mV")
            ax.grid()
            ax.legend()
            print(out[:, 0])

        ax2.errorbar(xs, out[:, 0] - xs * 2500 / 4096, out[:, 1], marker="s", markersize=0.5, ls="", label=f"Ch {ch}")
        ax2.set_ylabel("residual / mV")
        ax2.set_xlabel("input ADC")
        ax2.grid()

        up = numpy.argmin(numpy.abs(out[:, 0] - 700))
        print(xs[up])
        ax2.axvspan(xs[up], xs[-1], color="grey", alpha=0.5)

        low = numpy.argmin(numpy.abs(out[:, 0] + 1300))
        print(xs[low])
        ax2.axvspan(xs[0], xs[low], color="grey", alpha=0.5)

        if ax1 is None:
            title = self.__path.replace('.root', f'_ch{ch}')
            ax.set_title(title)
            fig.tight_layout()
            # plt.show()
            name = self.__path.replace('.root', f'_ch{ch}_test.png')
            plt.savefig(name, transparent=False)
        else:
            return ax2


def unpack_cal_parameters(cal_file : uproot.ReadOnlyDirectory) -> numpy.ndarray:
    """
    Function that reads out 9th order polynomial parameters from the calibration file

    Returns
    -------
    coef : numpy.ndarray of shape (24 * 4096, 10)
    """
    # stack is needed to convert ndarray of ndarrays to 'normally shaped' array, otherwise you have nested ndarrays
    coef = numpy.stack(cal_file["coeffs_tree/coeff"].array(library = 'np'))
    return coef


def unpack_cal_residuals(cal_file : uproot.ReadOnlyDirectory) -> numpy.ndarray:
    """
    Function that reads out the residuals from the calibration file

    Returns
    -------
    (v_residuals, residuals) : tuple of np arrays
        both v_residuals and residuals have shape (points, 24)

    """
    res = numpy.array([[ele.member("fX"), ele.member("fY")] for ele in 
            cal_file["aveResidGraph_tree/aveResidGraph"].array(library="np")])
    return res[:, 0], res[:, 1]



def unpack_raw_bias_scan(
        bias_scan : uproot.ReadOnlyDirectory,
        num_samples : Optional[int] = mattak.Dataset.AbstractDataset.NUM_DIGI_SAMPLES,
        num_channels : Optional[int] = mattak.Dataset.AbstractDataset.NUM_CHANNELS) -> tuple:
    """
    Parser for the raw bias scans, used when performing a "raw" voltage calibration
    (for testing purposes)

    Parameters
    ----------

    bias_scan : uproot.ReadOnlyDirectory
        Uproot file containing the bias scan data

    num_samples : int
        Number of samples of the digitizer

    num_channels : int
        Number of channels

    Returns
    -------
    vbias : numpy.ndarray of shape (points, 2)
        The voltage pedestals used when taking the bias scan
    adc : numpy.ndarray of shape (channels, samples, points)
        The measured adc counts
    """

    vbias = bias_scan["pedestals/vbias[2]"].array(library = "np")
    adc = bias_scan[f"pedestals/pedestals[{num_channels}][{num_samples}]"].array(library = "np").astype(numpy.float32)
    # (v_ped, channel, sample) -> (channel, sample, v_ped)
    adc = numpy.moveaxis(adc, (0, 1, 2), (2, 0, 1))
    return vbias, adc


def rescale_adc(
        vbias : numpy.ndarray, adc : numpy.ndarray, Vref : float = 1.5,
        num_samples : Optional[int] = mattak.Dataset.AbstractDataset.NUM_DIGI_SAMPLES,
        num_channels : Optional[int] = mattak.Dataset.AbstractDataset.NUM_CHANNELS) -> tuple:

    """
    Rescaling function to set the base pedestal ( 1.5 V ) as the origin

    Parameters
    ----------
    vbias : numpy.ndarray of shape (points, 2)
        The voltage pedestals used when taking the bias scan
    adc : numpy.ndarray of shape (channels, samples, points)
        The measured adc counts

    Returns
    -------
    vbias_rescaled, adc_rescaled : tuple of numpy.ndarrays
        The rescaled arrays
    """
    # The mattak src rescaled according to the first value GREATER than Vref, hence this function does the same
    vidx = [min([idx for idx, _ in enumerate(vbias[:, dac]) if vbias[idx, dac] >= Vref]) for dac in range(2)]
    adc_rescaled = numpy.zeros_like(adc)
    for ch in range(num_channels):
        dac = int(ch / 12)
        for s in range(num_samples):
            two_bins_around_pedestal = [vidx[dac] - 1, vidx[dac]]
            adc_rescaled[ch, s, :] = \
                adc[ch, s, :] - numpy.interp(
                    Vref, vbias[two_bins_around_pedestal, dac], adc[ch, s, two_bins_around_pedestal])

    vbias_rescaled = vbias - Vref
    return vbias_rescaled, adc_rescaled


def raw_calibrate(
        waveform_array : numpy.ndarray, vbias : numpy.ndarray, adc : numpy.ndarray,
        starting_window : Union[float, int],
        num_wf_samples : Optional[int] = mattak.Dataset.AbstractDataset.NUM_WF_SAMPLES,
        num_phys_samples : Optional[int] = mattak.Dataset.AbstractDataset.NUM_DIGI_SAMPLES
        ) -> numpy.ndarray:
    """
    Function that interpolates raw bias scans to perform ADC to voltage conversion
    (for testing purposes)
    """

    waveform_volt = numpy.zeros_like(waveform_array, dtype=float)

    vbias, adc = rescale_adc(vbias, adc)
    adc_cut = [[] for i in range(len(waveform_volt))]
    adc_cut[:12] = adc[:12, :, numpy.all([-1.3 < vbias[:, 0], vbias[:, 0] < 0.7], axis=0)]
    adc_cut[12:] = adc[12:, :, numpy.all([-1.3 < vbias[:, 1], vbias[:, 1] < 0.7], axis=0)]
    adc = adc_cut
    vbias =  numpy.array([[v for v in vbias[:, DAC] if -1.3 < v < 0.7] for DAC in range(2)]).T

    for c, wf_channel in enumerate(waveform_array):
        starting_window_channel = starting_window[c]
        # Reordering the parameters to match the correct starting window
        # (calibration is done per sample)
        adc_channel = adc[c]

        samples_idx = (128 * starting_window_channel + numpy.arange(num_wf_samples)) % num_wf_samples
        if starting_window_channel >= 16:
            samples_idx += num_wf_samples
        adc_channel = adc_channel[samples_idx]

        for s, (adc_wf, adc_bias) in enumerate(zip(wf_channel, adc_channel)):
            volt = numpy.interp(adc_wf, adc_bias, vbias[:, int(c/12)])
            waveform_volt[c, s] = volt

    return waveform_volt


def calibrate(
        waveform_array : numpy.ndarray, param : numpy.ndarray,
        vres : numpy.ndarray, res : numpy.ndarray, starting_window : Union[float, int],
        upsampling : bool = True,
        fit_min : float = -1.3, fit_max : float = 0.7, accuracy : float = 0.005,
        num_wf_samples : Optional[int] = mattak.Dataset.AbstractDataset.NUM_WF_SAMPLES,
        num_phys_samples : Optional[int] = mattak.Dataset.AbstractDataset.NUM_DIGI_SAMPLES
        ) -> numpy.ndarray:
    """
    The calibration function that transforms waveforms from ADC to voltage. Uses no caching

    Parameters
    ----------
    waveform_array : array of shape (24, 2048)
        array of one waveform
    param : array of shape (24 * 4096, 10)
        the parameters found in the calibration file
    vres : array of shape (points, 24)
        the voltage points of the residuals
    res : array of shape (points, 24)
        the ADC values of the residuals
    starting_window : int | float
        the sample on which the run started
    fit_min : float
        lower bound of original fit used on the bias scan
    fit_max : float
        upper bound of original fit used on the bias scan
    accuracy : float
        nr of voltage steps in table of the calibration function

    Return
    ------
    waveform_volt : array of shape (24, 2048)
        calibrated waveform in volt
    """
    waveform_volt = numpy.zeros_like(waveform_array, dtype=float)

    if upsampling:
        # "discrete" inverse
        vsamples = numpy.arange(fit_min, fit_max, accuracy)
        # residuals split over DACs
        res = [numpy.interp(vsamples, vres[i], res[i]) for i in range(len(res))]
    else:
        vsamples = vres[0]  # assuming the same for all channels

    for c, wf_channel in enumerate(waveform_array):
        starting_window_channel = starting_window[c]
        # Reordering the parameters to match the correct starting window
        # (calibration is done per sample)
        param_channel = param[num_phys_samples * c : num_phys_samples * (c + 1)]

        samples_idx = (128 * starting_window_channel + numpy.arange(num_wf_samples)) % num_wf_samples
        if starting_window_channel >= 16:
            samples_idx += num_wf_samples

        param_channel = param_channel[samples_idx]
        for s, (adc, p) in enumerate(zip(wf_channel, param_channel)):
            # discrete inverse
            adcsamples = numpy.polyval(p[::-1], vsamples) + res[c]
            volt = numpy.interp(adc, adcsamples, vsamples, left = fit_min, right = fit_max)
            waveform_volt[c, s] = volt

    return waveform_volt


if __name__ == "__main__":
    from mattak.Dataset import Dataset
    import os
    run = 1144
    station_id = 23
    channel_id = 0
    run_path = os.environ["RNO_G_DATA"] + "/" + f"station{station_id}/" + f"run{run}"
    vc_name = [file for file in os.listdir(run_path) if file.startswith("volCal") and file.endswith(".root")][0]

    vc = VoltageCalibration(run_path + "/" + vc_name)
    vc.plot_ch(ch=channel_id)

    ds = Dataset(0, 0, run_path, backend="uproot")
    wf = ds.wfs(calibrated=False)
    vc = VoltageCalibration(run_path + "/" + vc_name, caching=False)
    vc(wf, numpy.zeros(24, dtype=int))
    print(wf)
