import mattak.Dataset

from typing import Union, Optional, List
import numpy
import uproot


class VoltageCalibration(object):

    def __init__(self, path : str, upsample_residuals : bool = True,
                 caching : bool = True,
                 fit_min : float = -1.3, fit_max : float = 0.7, accuracy : float = 0.005):
        """ Helper class to apply the voltage calibration for the uproot backend.

        Parameters
        ----------

        path : str
            Path to the calibration root file. You can either pass a bias scan file or a
            voltage calibration constant file.

        upsample_residuals : bool (Default: True)
            If true, upsample the residuals used in the calibration. Not used when passing a
            bias scan file.

        caching : bool (Default: True)
            If True, the adc tables for the calibration are cached.

        fit_min : float (Default: -1.3 V)
            Lower bound for the upsample voltage range (in Volt).

        fit_max : float (Default: 0.7 V)
            Upper bound for the upsample voltage range (in Volt).

        accuracy : float (Default: 0.005)
            Step width of the upsampled voltage vector.
        """

        self.NUM_CHANNELS = mattak.Dataset.AbstractDataset.NUM_CHANNELS
        self.NUM_WF_SAMPLES = mattak.Dataset.AbstractDataset.NUM_WF_SAMPLES
        self.NUM_DIGI_SAMPLES = mattak.Dataset.AbstractDataset.NUM_DIGI_SAMPLES
        self.__upsample_residuals = upsample_residuals
        self.__caching = caching

        self.__adc_table_voltage = None
        if self.__caching:
            self.__adc_table = numpy.array([None] * self.NUM_CHANNELS, dtype=object)

        self.__adc_table_voltage = None
        self.__adc_table = numpy.array([None] * self.NUM_CHANNELS, dtype=object)

        self.__path = path
        if path.endswith(".root"):
            self.cal_file = uproot.open(path)
            if "coeffs_tree" in self.cal_file:
                self.full_bias_scan = False
                self.__cal_param = unpack_cal_parameters(self.cal_file)
                self.__cal_residuals_v, self.__cal_residuals_adc = unpack_cal_residuals(self.cal_file)

                self.__cal_residuals_v = self.__cal_residuals_v.T
                self.__cal_residuals_adc = self.__cal_residuals_adc.T

                if numpy.any(self.__cal_residuals_v[0] != self.__cal_residuals_v[1]):
                    raise ValueError("The pedestal voltage of the bias scan is different for the two DAC, "
                                     "the code expects them to be the same!")

                if self.__upsample_residuals:
                    vsamples = numpy.arange(fit_min, fit_max, accuracy)
                    # residuals split over DACs
                    ressamples = (numpy.interp(vsamples, self.__cal_residuals_v[0], self.__cal_residuals_adc[0]),
                                  numpy.interp(vsamples, self.__cal_residuals_v[1], self.__cal_residuals_adc[1]))

                    self.__set_adc_table_voltage(vsamples)
                    self.__cal_residuals_adc = ressamples
                else:
                    self.__set_adc_table_voltage(self.__cal_residuals_v[0])

            elif "pedestals" in self.cal_file:
                self.full_bias_scan = True
                self.__vbias, self.__adc = unpack_raw_bias_scan(self.cal_file)

                if numpy.any(self.__vbias[:, 0] != self.__vbias[:, 1]):
                    raise ValueError("The pedestal voltage of the bias scan is different for the two DAC, "
                                        "the code expects them to be the same!")

            else:
                raise ValueError("No 'coeffs_tree' or 'pedestals' keys found in the root file")
        else:
            raise ValueError(f"{path} is not recognized as a root file")


    def __set_adc_table_voltage(self, value):
        """ Set the voltage vector which correspond to the `adc_table` for calibration """
        if self.__adc_table_voltage is None:
            self.__adc_table_voltage = value
        else:
            if self.__adc_table_voltage != value:
                raise ValueError("The voltage vector for the adc table changed!")


    def __get_adc_table(self, channel, sample):
        """ Returns the cached table of ADC counts for a given voltage """

        if self.__adc_table[channel] is None:

            if not self.full_bias_scan:
                param_channel = self.__cal_param[self.NUM_DIGI_SAMPLES * channel:self.NUM_DIGI_SAMPLES * (channel + 1)]
                # checked that self.__cal_residuals_v is equal for both DACs

                adcsamples = numpy.array([numpy.polyval(p[::-1], self.__adc_table_voltage) for p in param_channel])

                adcsamples[:2048] += self.__cal_residuals_adc[0]
                adcsamples[2048:] += self.__cal_residuals_adc[1]

                self.__adc_table[channel] = numpy.asarray(adcsamples, dtype=float)

            else:
                vbias, adc = rescale_adc(self.__vbias, self.__adc)
                adc_cut = [[] for i in range(self.NUM_CHANNELS)]
                adc_cut[:12] = adc[:12, :, numpy.all([-1.3 < vbias[:, 0], vbias[:, 0] < 0.7], axis=0)]
                adc_cut[12:] = adc[12:, :, numpy.all([-1.3 < vbias[:, 1], vbias[:, 1] < 0.7], axis=0)]
                adc = numpy.array(adc_cut)
                vbias =  numpy.array([[v for v in vbias[:, DAC] if -1.3 < v < 0.7] for DAC in range(2)])

                self.__set_adc_table_voltage(vbias[0])
                self.__adc_table = adc

        return self.__adc_table[channel][sample]


    def calibrate(self, waveform_array : numpy.ndarray, starting_window : Union[float, int],
                  channels : Optional[List[int]] = None,
                  fit_min : float = -1.3, fit_max : float = 0.7) -> numpy.ndarray:
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
        fit_min : float
            lower bound of original fit used on the bias scan
        fit_max : float
            upper bound of original fit used on the bias scann

        Returns
        -------
        waveform_volt : array of shape (num_channels, num_wf_samples)
            calibrated waveform in volt
        """

        channels = channels or list(range(self.NUM_CHANNELS))
        waveform_volt = numpy.zeros_like(waveform_array, dtype=float)

        for idx, (ch, starting_window_channel, wf_channel) in enumerate(zip(channels, starting_window, waveform_array)):
            # residuals split over DACs
            samples_idx = (128 * starting_window_channel + numpy.arange(self.NUM_WF_SAMPLES)) % self.NUM_WF_SAMPLES
            if starting_window_channel >= 16:
                samples_idx += self.NUM_WF_SAMPLES

            for sample_wf, (sample_lab, adc) in enumerate(zip(samples_idx, wf_channel)):
                adcsamples = self.__get_adc_table(ch, sample_lab)
                volt = numpy.interp(adc, adcsamples, self.__adc_table_voltage, left = fit_min, right = fit_max)
                waveform_volt[idx, sample_wf] = volt

        return waveform_volt


    def __call__(self, waveform_array : numpy.ndarray, starting_window : Union[float, int]) -> numpy.ndarray:
        if self.__caching:
            return self.calibrate(waveform_array, starting_window)
        else:
            return calibrate(waveform_array, self.__cal_param, [self.__adc_table_voltage, self.__adc_table_voltage],
                             self.__cal_residuals_adc, starting_window, upsampling=False)  # already upsampled

    def plot_ch(self, xs=numpy.linspace(-1000, 1000, 100), ch=0):

        from matplotlib import pyplot as plt

        fig, (ax, ax2) = plt.subplots(2, 1, height_ratios=[3, 1], sharex=True)

        # tables = numpy.array([self.__get_adc_table(ch, i) for i in range(self.NUM_DIGI_SAMPLES)])
        out = []
        out2 = []
        for x in xs:
            # vs = [numpy.interp(x, t, self.__adc_table_voltage) * 1000 for t in tables]
            # out.append([numpy.mean(vs), numpy.std(vs)])
            vs2 = self.calibrate([numpy.full(self.NUM_WF_SAMPLES, x)], [0], channels=[ch]) * 1000
            out2.append([numpy.mean(vs2), numpy.std(vs2)])

        out = numpy.array(out2)

        ax.errorbar(xs, out[:, 0], out[:, 1], ls="", label=f"Ch {ch} (mean over all samples)")
        ax.plot(xs, xs * 2500 / 4096, "k--", lw=1, label="2500 / 4096")
        ax2.set_xlabel("input ADC")
        ax.set_ylabel("output voltage / mV")

        ax2.errorbar(xs, out[:, 0] - xs * 2500 / 4096, out[:, 1], ls="", label=f"Ch {ch}")
        ax2.set_ylabel("residual / mV")

        ax.grid()
        ax2.grid()
        ax.legend()
        title = self.__path.replace('.root', f'_ch{ch}')
        ax.set_title(title)
        fig.tight_layout()
        plt.show()
        # name = self.__path.replace('.root', f'_ch{ch}_test.png')
        # plt.savefig(name, transparent=False)


def unpack_cal_parameters(cal_file : uproot.ReadOnlyDirectory) -> numpy.ndarray:
    """
    Function that reads out 9th order polynomial parameters from the calibration file

    Returns
    -------
    coef : np.ndarray of shape (24 * 4096, 10)
    """
    # stack is needed to convert ndarray of ndarrays to 'normally shaped' array, otherwise you have nested ndarrays
    coef = numpy.stack(cal_file["coeffs_tree/coeff"].array(library = 'np'))
    return coef


def unpack_cal_residuals(cal_file : uproot.ReadOnlyDirectory) -> numpy.ndarray:
    """
    Function that reads out the residuals from the calibration file

    Returns
    -------
    (v_residuals, residuals) : tuple of numpy arrays
        both v_residuals and residuals have shape (points, 2)

    """
    vres_dac1 = cal_file["aveResid_dac1"].member("fX")
    vres_dac2 = cal_file["aveResid_dac2"].member("fX")
    residual_dac1 = cal_file["aveResid_dac1"].member("fY")
    residual_dac2 = cal_file["aveResid_dac2"].member("fY")

    return numpy.stack(numpy.array([vres_dac1, vres_dac2]), axis = -1), \
        numpy.stack(numpy.array([residual_dac1, residual_dac2]), axis = -1)


def unpack_raw_bias_scan(bias_scan : uproot.ReadOnlyDirectory,
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


def rescale_adc(vbias : numpy.ndarray, adc : numpy.ndarray, Vref : float = 1.5,
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
                adc[ch, s, :] - numpy.interp(Vref, vbias[two_bins_around_pedestal, dac], adc[ch, s, two_bins_around_pedestal])

    vbias_rescaled = vbias - Vref
    return vbias_rescaled, adc_rescaled


def raw_calibrate(waveform_array : numpy.ndarray, vbias : numpy.ndarray, adc : numpy.ndarray,
                  starting_window : Union[float, int],
                  num_wf_samples : Optional[int] = mattak.Dataset.AbstractDataset.NUM_WF_SAMPLES,
                  num_phys_samples : Optional[int] = mattak.Dataset.AbstractDataset.NUM_DIGI_SAMPLES
                  )-> numpy.ndarray:
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


def calibrate(waveform_array : numpy.ndarray, param : numpy.ndarray,
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
    vres : array of shape (points, 2)
        the voltage points of the residuals, shape
    res : array of shape (points, 2)
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
        res = (numpy.interp(vsamples, vres[0], res[0]), numpy.interp(vsamples, vres[1], res[1]))
    else:
        vsamples = vres[0]  # assuming both are the same

    for c, wf_channel in enumerate(waveform_array):
        starting_window_channel = starting_window[c]
        # Reordering the parameters to match the correct starting window
        # (calibration is done per sample)
        param_channel = param[num_phys_samples * c : num_phys_samples * (c + 1)]

        samples_idx = (128 * starting_window_channel + numpy.arange(num_wf_samples)) % num_wf_samples
        if starting_window_channel >= 16:
            samples_idx += num_wf_samples

        dac = int(c / 12)
        param_channel = param_channel[samples_idx]
        for s, (adc, p) in enumerate(zip(wf_channel, param_channel)):
            # discrete inverse
            adcsamples = numpy.polyval(p[::-1], vsamples) + res[dac]
            volt = numpy.interp(adc, adcsamples, vsamples, left = fit_min, right = fit_max)
            waveform_volt[c, s] = volt

    return waveform_volt
