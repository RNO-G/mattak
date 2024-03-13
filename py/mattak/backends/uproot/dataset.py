import os
import uproot
import configparser

import mattak.Dataset
from typing import Union, Optional, Tuple, Generator, Callable, Sequence
import numpy
import math
import logging


# Dublicated from Dataset.cc
waveform_tree_names = ["waveforms", "wfs", "wf", "waveform"]
header_tree_names = ["hdr", "header", "hd", "hds", "headers"]
daqstatus_tree_names = ["daqstatus", "ds", "status"]


def read_tree(ur_file, tree_names):
    """
    Parameters
    ----------

    ur_file: uproot.open()
        Input file

    tree_names: list(str)
        List of possible tree names, the first found is returned

    Returns
    -------

    tree: uproot.tree
    """

    for tree_name in tree_names:
        if tree_name in ur_file:
            tree = ur_file[tree_name]
            return tree, tree_name


class Dataset(mattak.Dataset.AbstractDataset):

    def __init__(self, station : int, run : int, data_path : str, verbose : bool = False,
                 skip_incomplete : bool = True, read_daq_status : bool = True,
                 read_run_info : bool = True, preferred_file : Optional[str] = None,
                 voltage_calibration : Optional[str] = None):
        """
        Uproot backend for the python interface of the mattak Dataset. See further information in
        `mattak.Dataset.Dataset` about the arguments `station`, `run`, `data_path` (called `data_dir` there),
        and `preferred_file`.

        Parameters
        ----------

        station: int
            Station Id of the data to be read.

        run: int
            Run number of the data to be read

        data_path: str
            Path of the run directory or root file

        verbose: bool
            Enable verbose debug log

        skip_incomplete: bool
            Default: True

        read_daq_status: bool
            Enable reading of daqstatus.root. Only possible  when a complete run directory is passed. (Default: True)

        read_run_info: bool
            Enable reading of the run information (sampling rate). Only possible when a run directory is passed
            which contains the runinfo.txt file. (Default: True, unless you pass a specific file as input)

        preferred_file: str
            Specify a prefered file name to load.

        voltage_calibration : str
            Path to a voltage calibration file. If None, check for file in run directory.
        """

        self.backend = "uproot"
        self.__verbose = verbose
        self.__read_daq_status = read_daq_status
        self.__read_run_info = read_run_info

        # special case where data_dir is a file
        if os.path.isfile(data_path):
            self.data_path_is_file = True
            self.rundir = os.path.dirname(data_path)
        else:
            if data_path.endswith(".root"):  # catch the case where it was intented to pass a file
                raise FileNotFoundError(f"Could not find the file {data_path}")

            self.data_path_is_file = False
            # special case where we load a directory instead of a station/run
            if station == 0 and run == 0:
                self.rundir = data_path
            else:
                self.rundir = f"{data_path}/station{station}/run{run}"

        if skip_incomplete is False and self.data_path_is_file:
            logging.warining("`skip_incomplete == False` is incompatible with data_dir as file. "
                             "Set `skip_incomplete == True`")
            skip_incomplete = True

        self.skip_incomplete = skip_incomplete

        # this duplicates a bunch of C++ code in mattak::Dataset
        # check for full or partial run by looking for waveforms.root
        # Only open files/trees, do not access data.
        # First we try to load a combined file if the `data_path`` is a file or
        # `preferred_file` is set. Otherwise we try to load a full run (i.e.,
        # waveforms.root, headers.root, ...). If that does not work we try to read
        # combined.root

        self.combined_tree = None

        if self.data_path_is_file:
            self.combined_tree = uproot.open(f"{data_path}:combined")
        elif preferred_file not in [None, ""]:
            if preferred_file.endswith(".root"):
                preferred_file = preferred_file[:-5]  # strip ".root"

            if os.path.exists(f"{self.rundir}/{preferred_file}.root"):
                self.combined_tree = uproot.open(f"{self.rundir}/{preferred_file}.root:combined")
            else:
                logging.warning(f"Could not find prefered file {self.rundir}/{preferred_file}.root. "
                                "Revert to default behaviour ...")

        # if we didn't load the combined_tree already, try to load full tree
        if self.combined_tree is None:
            try:
                self.wf_file = uproot.open("%s/waveforms.root" % (self.rundir))
                if self.__verbose:
                    print ("Open waveforms.root (Found full run folder) ...")

                self.full = True

                self.wf_tree, self.wf_branch = read_tree(self.wf_file, waveform_tree_names)
                self._wfs = self.wf_tree[self.wf_branch]

                self.hd_file = uproot.open("%s/headers.root" % (self.rundir))
                self.hd_tree, self.hd_branch = read_tree(self.hd_file, header_tree_names)
                self._hds = self.hd_tree[self.hd_branch]

                if self.__read_daq_status:
                    self.ds_file = uproot.open("%s/daqstatus.root" % (self.rundir))
                    self.ds_tree, self.ds_branch = read_tree(self.ds_file, daqstatus_tree_names)
                    self._dss = self.ds_tree[self.ds_branch]
            except Exception:
                self.full = False
        else:
            self.full = False  # because combined_tree is not None

        # we haven't already loaded the full tree
        if not self.full:
            if self.combined_tree is None: # we didn't already load our preference
                self.combined_tree = uproot.open(f"{self.rundir}/combined.root:combined")
                if self.__verbose:
                    print("Found combined file")

            self._wfs, self.wf_branch = read_tree(self.combined_tree, waveform_tree_names)

            # build an index of the waveforms we do have
            if not self.skip_incomplete:
                wfs_included = self._wfs['event_number'].array()
                self.events_with_waveforms = {ev: idx for idx, ev in enumerate(wfs_included)}

                self.full_head_file = uproot.open(f"{self.rundir}/headers.root")
                self.full_head_tree,_ = read_tree(self.full_head_file, header_tree_names)

                if self.__read_daq_status:
                    self.full_daq_file = uproot.open(f"{self.rundir}/daqstatus.root")
                    self.full_daq_tree, _ = read_tree(self.full_daq_file, daqstatus_tree_names)


            # Get header and daq information from combined file or from full run
            hd_tree = self.combined_tree if skip_incomplete else self.full_head_tree
            self._hds, self.hd_branch = read_tree(hd_tree, header_tree_names)

            if self.__read_daq_status:
                ds_tree = self.combined_tree if skip_incomplete else self.full_daq_tree
                self._dss, self.ds_branch =  read_tree(ds_tree, daqstatus_tree_names)

        if station == 0 and run == 0 or self.data_path_is_file:
            self.station = self._hds['station_number'].array(entry_start=0, entry_stop=1)[0]
            self.run = self._hds['run_number'].array(entry_start=0, entry_stop=1)[0]
        else:
            self.station = station
            self.run = run

        self.data_path = data_path
        self.setEntries(0)

        self.run_info = None
        self.runfile = None
        if self.__read_run_info and self.runfile is None:
            # try to get the run info, if we're using combined tree, try looking in there
            # doh, uproot can't read the runinfo ROOT files... let's parse the text files instead

            if os.path.exists("%s/aux/runinfo.txt" % (self.rundir)):
                with open("%s/aux/runinfo.txt" % (self.rundir)) as fruninfo:
                    # we'll abuse configparser to read the runinfo, but we have to add a dummy
                    # section to properly abuse it
                    config = configparser.ConfigParser()
                    config.read_string('[dummy]\n' + fruninfo.read())
                    self.run_info = config['dummy']

        if voltage_calibration is None:
            # do find_VC here
            time = self._hds['trigger_time'].array()[0]
            cal = mattak.Dataset.find_voltage_calibration(self.rundir, self.station, time)
            calibration_files = [cal] if cal is not None else []
        elif isinstance(voltage_calibration, str):
            calibration_files = [voltage_calibration]
        else:
            raise TypeError(f"Unknown type for voltage calibration in uproot backend ({voltage_calibration})")

        self.__adc = None
        self.__cal_param = None

        if len(calibration_files):
            if calibration_files[0].endswith(".root"):
                self.cal_file = uproot.open(calibration_files[0])
                if "coeffs_tree" in self.cal_file:
                    self.__cal_param = unpack_cal_parameters(self.cal_file)
                    self.__cal_residuals_v, self.__cal_residuals_adc = unpack_cal_residuals(self.cal_file)

                    self.__cal_residuals_v = self.__cal_residuals_v.T
                    self.__cal_residuals_adc = self.__cal_residuals_adc.T

                    if numpy.any(self.__cal_residuals_v[0] != self.__cal_residuals_v[1]):
                        raise ValueError("The pedestal voltage of the bias scan is different for the two DAC, "
                                         "the code expects them to be the same!")

                elif "pedestals" in self.cal_file:
                    self.__vbias, self.__adc = unpack_raw_bias_scan(self.cal_file)

                    if numpy.any(self.__vbias[:, 0] != self.__vbias[:, 1]):
                        raise ValueError("The pedestal voltage of the bias scan is different for the two DAC, "
                                         "the code expects them to be the same!")

                else:
                    raise ValueError("No 'coeffs_tree' or 'pedestals' keys found in the root file")
            else:
                raise ValueError(f"{calibration_files[0]} is not recognized as a root file")

        self.has_calib = self.__cal_param is not None or self.__adc is not None

        self.__adc_table_voltage = None
        self.__adc_table = numpy.array([None] * 24, dtype=object)

        # keep it hard-coded for the moment
        self.__upsample_residuals = True

        if self.__cal_param is not None:
            if self.__upsample_residuals:
                vsamples = numpy.arange(-1.3, 0.7, 0.005)
                # residuals split over DACs
                ressamples = (numpy.interp(vsamples, self.__cal_residuals_v[0], self.__cal_residuals_adc[0]),
                            numpy.interp(vsamples, self.__cal_residuals_v[1], self.__cal_residuals_adc[1]))

                self.__set_adc_table_voltage(vsamples)
                self.__cal_residuals_adc = ressamples
            else:
                self.__set_adc_table_voltage(self.__cal_residuals_v[0])


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
            if self.__cal_param is not None:
                param_channel = self.__cal_param[4096 * channel : 4096 * (channel + 1)]

                # checked that self.__cal_residuals_v is equal for both DACs
                adcsamples = numpy.array([numpy.polyval(p[::-1], self.__adc_table_voltage) for p in param_channel])

                adcsamples[:2048] += self.__cal_residuals_adc[0]
                adcsamples[2048:] += self.__cal_residuals_adc[1]

                self.__adc_table[channel] = numpy.asarray(adcsamples, dtype="float32")

            elif self.__adc is not None:

                vbias, adc = rescale_adc(self.__vbias, self.__adc)
                adc_cut = [[] for i in range(24)]
                adc_cut[:12] = adc[:12, :, numpy.all([-1.3 < vbias[:, 0], vbias[:, 0] < 0.7], axis=0)]
                adc_cut[12:] = adc[12:, :, numpy.all([-1.3 < vbias[:, 1], vbias[:, 1] < 0.7], axis=0)]
                adc = numpy.array(adc_cut)
                vbias =  numpy.array([[v for v in vbias[:, DAC] if -1.3 < v < 0.7] for DAC in range(2)])

                self.__set_adc_table_voltage(vbias[0])
                self.__adc_table = adc

            else:
                raise ValueError("No calibration data available.")

        return self.__adc_table[channel][sample]

    def eventInfo(self) -> Union[Optional[mattak.Dataset.EventInfo], Sequence[Optional[mattak.Dataset.EventInfo]]]:
        kw = dict(entry_start = self.first, entry_stop = self.last)

        station = self._hds['station_number'].array(**kw)
        run = self._hds['run_number'].array(**kw)
        eventNumber = self._hds['event_number'].array(**kw)
        readoutTime = self._hds['readout_time'].array(**kw)
        triggerTime = self._hds['trigger_time'].array(**kw)
        triggerInfo = self._hds['trigger_info'].array(**kw)
        pps = self._hds['pps_num'].array(**kw)
        sysclk = self._hds['sysclk'].array(**kw)
        sysclk_lastpps = self._hds['sysclk_last_pps'].array(**kw)
        sysclk_lastlastpps = self._hds['sysclk_last_last_pps'].array(**kw)

        if self.__read_daq_status:
            radiantThrs = numpy.array(self._dss['radiant_thresholds[24]'])
            lowTrigThrs = numpy.array(self._dss['lt_trigger_thresholds[4]'])

        if self.run_info is not None:
            sampleRate = float(self.run_info['radiant-samplerate']) / 1000
        else:
            sampleRate = None

        # um... yeah, that's obvious
        radiantStartWindows = self.get_windows(kw)

        infos = []
        info = None  # if range(0)
        for i in range(self.last-self.first):

            triggerType  = "UNKNOWN"
            if triggerInfo[i]['trigger_info.radiant_trigger']:
                which = triggerInfo[i]['trigger_info.which_radiant_trigger']
                if which == -1:
                    which = "X"
                triggerType = "RADIANT" + str(which)
            elif triggerInfo[i]['trigger_info.lt_trigger']:
                triggerType = "LT"
            elif triggerInfo[i]['trigger_info.force_trigger']:
                triggerType = "FORCE"
            elif triggerInfo[i]['trigger_info.pps_trigger']:
                triggerType = "PPS"

            info = mattak.Dataset.EventInfo(
                eventNumber = eventNumber[i],
                station = station[i],
                run = run[i],
                readoutTime = readoutTime[i],
                triggerTime = triggerTime[i],
                triggerType = triggerType,
                sysclk = sysclk[i],
                sysclkLastPPS = (sysclk_lastpps[i], sysclk_lastlastpps[i]),
                pps = pps[i],
                radiantStartWindows = radiantStartWindows[i],
                sampleRate = sampleRate,
                radiantThrs=radiantThrs[i] if self.__read_daq_status else None,
                lowTrigThrs=lowTrigThrs[i] if self.__read_daq_status else None
            )

            infos.append(info)

        if not self.multiple:
            return info
        else:
            return infos


    def N(self) -> int:
        return self._hds.num_entries

    def get_windows(self, kw):
        return self._hds['trigger_info/trigger_info.radiant_info.start_windows[24][2]'].array(**kw)

    def wfs(self, calibrated : bool = False, raw_calibration = False) -> Optional[numpy.ndarray]:
        if calibrated and not self.has_calib:
            raise ValueError("You requested a calibrated waveform but no calibration is available")

        # assert(not calibrated) # not implemented yet
        kw = dict(entry_start=self.first, entry_stop=self.last, library='np')

        w = None
        if self.full or self.skip_incomplete:
            w = self._wfs['radiant_data[24][2048]'].array(**kw)
            starting_window = self.get_windows(kw)
        elif not self.multiple:
            # if you only selected one event and have an incomplete dataset
            if self.first in self.events_with_waveforms:
                idx = self.events_with_waveforms[self.first]
                w = self._wfs['radiant_data[24][2048]'].array(entry_start=idx, entry_stop=idx+1, library='np')
                starting_window = self.get_windows(dict(entry_start=idx, entry_stop=idx+1, library='np'))
        else:
            # so ... we need to loop through and find which things we have actually have waveforms
            # start by allocating the output
            w = numpy.zeros((self.last - self.first, 24, 2048), dtype='float64' if calibrated else 'int16')
            starting_window = numpy.zeros((self.last - self.first))
            # now figure out how much of the data array we need
            wf_start = None
            wf_end = None

            # store the indices that will be non-zero
            wf_idxs = []
            for i in range(self.first, self.last):
                if i in self.events_with_waveforms:
                    wf_idxs.append(i - self.first)
                    # these are the start and stop of our array we need to load
                    if wf_start is None:
                        wf_start = self.events_with_waveforms[i]

                    wf_end = self.events_with_waveforms[i] + 1

            if len(wf_idxs):
                w[wf_idxs] = self._wfs['radiant_data[24][2048]'].array(entry_start=wf_start, entry_stop=wf_end, library='np')
                starting_window[wf_idxs] = self.get_windows(dict(entry_start=wf_start, entry_stop=wf_end, library='np'))

        # calibration
        starting_window = starting_window[:, :, 0]
        if calibrated:
            if self.__cal_param is None and self.__adc is None:
                raise ValueError("Calibration not available")
            # this can run now both normal and raw calibration
            w = numpy.array([self.calibrate(ele, starting_window[i]) for i, ele in enumerate(w)])

        elif raw_calibration:
            # This still uses the slow version of the raw calibration
            if self.__adc is None:
                raise ValueError("Calibration not available")
            w = numpy.array([
                raw_calibrate(ele, self.__vbias, self.__adc, starting_window[i]) for i, ele in enumerate(w)])

        w = numpy.asarray(w, dtype=float)

        if self.multiple:
            return w

        return None if w is None else w[0]

    def _iterate(self, start : int , stop : int , calibrated: bool,  max_in_mem : int,
                 selector: Optional[Callable[[mattak.Dataset.EventInfo],bool]] = None) \
                    -> Generator[Tuple[Optional[mattak.Dataset.EventInfo], Optional[numpy.ndarray]],None,None]:

       # cache current values given by setEntries(..)
        original_entry : Union[int, Tuple[int,int]] = (self.first, self.last) if self.multiple else self.entry

        # determine in how many batches we want to access the data given how much events we want to load into the RAM at once
        n_batches = math.ceil((stop - start) / max_in_mem)

        for i_batch in range(n_batches):

            # looping over the batches defining the start and stop index
            batch_start = start + i_batch * max_in_mem
            batch_stop = min(stop, batch_start + max_in_mem)

            self.setEntries((batch_start, batch_stop))

            # load events from file
            w = self.wfs(calibrated)
            e = self.eventInfo()

            # we modified the internal data pointers with the prev. call of self.setEntries(...)
            # this is intransparent for the outside world and has to be reverted
            self.setEntries(original_entry)

            for idx in range(batch_stop - batch_start):
                if selector is not None:
                    if selector(e[idx]):
                        yield e[idx], w[idx]
                else:
                    yield e[idx], w[idx]

    def calibrate(self, waveform_array : numpy.ndarray, starting_window : Union[float, int],
                fit_min : float = -1.3, fit_max : float = 0.7, accuracy : float = 0.005) -> numpy.ndarray:
        """
        The calibration function that transforms waveforms from ADC to voltage

        Parameters
        ----------
        waveform_array : array of shape (24, 2048)
            array of one waveform
        starting_window : int | float
            the sample on which the run started
        fit_min : float
            lower bound of original fit used on the bias scan
        fit_max : float
            upper bound of original fit used on the bias scan
        accuracy : float
            nr of voltage steps in table of the calibration function

        Returns
        -------
        waveform_volt : array of shape (24, 2048)
            calibrated waveform in volt
        """

        waveform_volt = numpy.zeros((24, 2048))
        for c, wf_channel in enumerate(waveform_array):
            # residuals split over DACs

            starting_window_channel = starting_window[c]
            samples_idx = (128 * starting_window_channel + numpy.arange(2048)) % 2048
            if starting_window_channel >= 16:
                samples_idx += 2048

            for sample_wf, (sample_lab, adc) in enumerate(zip(samples_idx, wf_channel)):
                adcsamples = self.__get_adc_table(c, sample_lab)
                volt = numpy.interp(adc, adcsamples, self.__adc_table_voltage, left = fit_min, right = fit_max)
                waveform_volt[c, sample_wf] = volt

        return waveform_volt


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

    vres_dac1 = cal_file["aveResid_dac1"].values(axis = 0)
    vres_dac2 = cal_file["aveResid_dac2"].values(axis = 0)
    residual_dac1 = cal_file["aveResid_dac1"].values(axis = 1)
    residual_dac2 = cal_file["aveResid_dac2"].values(axis = 1)
    return numpy.stack(numpy.array([vres_dac1, vres_dac2]), axis = -1), \
        numpy.stack(numpy.array([residual_dac1, residual_dac2]), axis = -1)


def unpack_raw_bias_scan(bias_scan : uproot.ReadOnlyDirectory) -> tuple:
    """
    Parser for the raw bias scans, used when performing a "raw" voltage calibration
    (for testing purposes)

    Returns
    -------
    vbias : numpy.ndarray of shape (points, 2)
        The voltage pedestals used when taking the bias scan
    adc : numpy.ndarray of shape (channels, samples, points)
        The measured adc counts
    """

    vbias = bias_scan["pedestals/vbias[2]"].array(library = "np")
    adc = bias_scan["pedestals/pedestals[24][4096]"].array(library = "np").astype(numpy.float32)
    # (v_ped, channel, sample) -> (channel, sample, v_ped)
    adc = numpy.moveaxis(adc, (0, 1, 2), (2, 0, 1))
    return vbias, adc


def rescale_adc(vbias : numpy.ndarray, adc : numpy.ndarray, Vref = 1.5) -> tuple:
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
    for ch in range(24):
        dac = int(ch / 12)
        for s in range(4096):
            two_bins_around_pedestal = [vidx[dac] - 1, vidx[dac]]
            adc_rescaled[ch, s, :] = \
                adc[ch, s, :] - numpy.interp(Vref, vbias[two_bins_around_pedestal, dac], adc[ch, s, two_bins_around_pedestal])

    vbias_rescaled = vbias - Vref
    return vbias_rescaled, adc_rescaled


def raw_calibrate(waveform_array : numpy.ndarray, vbias : numpy.ndarray, adc : numpy.ndarray,
                        starting_window : Union[float, int]) -> numpy.ndarray:
    """
    Function that interpolates raw bias scans to perform ADC to voltage conversion
    (for testing purposes)
    """

    waveform_volt = numpy.zeros((24, 2048))

    vbias, adc = rescale_adc(vbias, adc)
    adc_cut = [[] for i in range(24)]
    adc_cut[:12] = adc[:12, :, numpy.all([-1.3 < vbias[:, 0], vbias[:, 0] < 0.7], axis=0)]
    adc_cut[12:] = adc[12:, :, numpy.all([-1.3 < vbias[:, 1], vbias[:, 1] < 0.7], axis=0)]
    adc = adc_cut
    vbias =  numpy.array([[v for v in vbias[:, DAC] if -1.3 < v < 0.7] for DAC in range(2)]).T

    for c, wf_channel in enumerate(waveform_array):
        starting_window_channel = starting_window[c]
        # Reordering the parameters to match the correct starting window
        # (calibration is done per sample)
        adc_channel = adc[c]

        samples_idx = (128 * starting_window_channel + numpy.arange(2048)) % 2048
        if starting_window_channel >= 16:
            samples_idx += 2048
        adc_channel = adc_channel[samples_idx]

        for s, (adc_wf, adc_bias) in enumerate(zip(wf_channel, adc_channel)):
            volt = numpy.interp(adc_wf, adc_bias, vbias[:, int(c/12)])
            waveform_volt[c, s] = volt

    return waveform_volt


def calibrate(waveform_array : numpy.ndarray, param : numpy.ndarray,
              vres : numpy.ndarray, res : numpy.ndarray, starting_window : Union[float, int],
              fit_min : float = -1.3, fit_max : float = 0.7, accuracy : float = 0.005) -> numpy.ndarray:
    """
    The calibration function that transforms waveforms from ADC to voltage

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

    # "discrete" inverse
    vsamples = numpy.arange(fit_min, fit_max, accuracy)
    waveform_volt = numpy.zeros((24, 2048))
    # residuals split over DACs
    ressamples = (numpy.interp(vsamples, vres[0], res[0]), numpy.interp(vsamples, vres[1], res[1]))

    for c, wf_channel in enumerate(waveform_array):
        starting_window_channel = starting_window[c]
        # Reordering the parameters to match the correct starting window
        # (calibration is done per sample)
        param_channel = param[4096 * c : 4096 * (c + 1)]

        samples_idx = (128 * starting_window_channel + numpy.arange(2048)) % 2048
        if starting_window_channel >= 16:
            samples_idx += 2048
        param_channel = param_channel[samples_idx]

        for s, (adc, p) in enumerate(zip(wf_channel, param_channel)):
            # discrete inverse
            adcsamples = numpy.polyval(p[::-1], vsamples) + ressamples[int(c/12)]
            volt = numpy.interp(adc, adcsamples, vsamples, left = fit_min, right = fit_max)
            waveform_volt[c, s] = volt

    return waveform_volt
