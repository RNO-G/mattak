import os
import uproot
import configparser

import mattak.Dataset
from .voltage_calibration import VoltageCalibration
from typing import Union, Optional, Tuple, Generator, Callable, Sequence, List
import numpy
import math
import logging


# Dublicated from Dataset.cc
waveform_tree_names = ["waveforms", "wfs", "wf", "waveform"]
header_tree_names = ["hdr", "header", "hd", "hds", "headers"]
daqstatus_tree_names = ["daqstatus", "ds", "status"]


def find_daq_status_index(event_readout_time, daq_readout_times):
    """
    Find the corresponding entry in the daq status for a particular event (readout time).
    Looking or the cloest entry before the event.
    """

    closest_idx = numpy.argmin(numpy.abs(daq_readout_times - event_readout_time))

    if daq_readout_times[closest_idx] > event_readout_time:
        closest_idx -= 1

    return closest_idx



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
                 voltage_calibration : Optional[str] = None, cache_calibration : Optional[bool] = True):
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

        self._radiantThrs = None
        self._lowTrigThrs = None

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
            voltage_calibration = mattak.Dataset.find_voltage_calibration_for_dataset(self)
        elif isinstance(voltage_calibration, str):
            pass
        elif isinstance(voltage_calibration, VoltageCalibration):
            self.vc = voltage_calibration
            self.has_calib = True

        else:
            raise TypeError(f"Unknown type for voltage calibration in uproot backend ({voltage_calibration})")

        if voltage_calibration is not None:
            if not self.has_calib:
                self.vc = VoltageCalibration(voltage_calibration, caching=cache_calibration)
                self.has_calib = True
        else:
            self.has_calib = False


    def eventInfo(self, override_skip_incomplete : Optional[bool] = None) -> Union[Optional[mattak.Dataset.EventInfo], Sequence[Optional[mattak.Dataset.EventInfo]]]:
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

        if self.__read_daq_status and self._radiantThrs is None:
            # The daq status is read asynchronously w.r.t. the events. Hence,
            # we always read all information at once and associate the to the
            # event below.
            self._radiantThrs = numpy.array(self._dss[f'radiant_thresholds[{self.NUM_CHANNELS}]'])
            self._lowTrigThrs = numpy.array(self._dss['lt_trigger_thresholds[4]'])
            self._readout_time_radiant = numpy.array(self._dss['readout_time_radiant'])
            self._readout_time_lt = numpy.array(self._dss['readout_time_lt'])

        try:
            sampleRate = self._wfs["mattak::IWaveforms/radiant_sampling_rate"].array(**kw) / 1000
            if not self.skip_incomplete:
                rate = numpy.unique(sampleRate)

                if len(rate) == 0:  # no waveforms available
                    raise uproot.exceptions.KeyInFileError("")  # HACK: let the except block handle it

                assert len(rate) == 1, "Sampling rate derived from waveforms in not unique. Can not extend to incomplete events ..."
                sampleRate = [rate[0]] * (self.last - self.first)

        except uproot.exceptions.KeyInFileError:
            if self.run_info is not None:
                sampleRate = float(self.run_info['radiant-samplerate']) / 1000
            else:
                sampleRate = 3.2  # GHz

            sampleRate = [sampleRate] * (self.last - self.first)

        # um... yeah, that's obvious
        radiantStartWindows = self._get_windows(kw)

        infos = []
        info = None  # if range(0)
        for i in range(self.last - self.first):

            if override_skip_incomplete is not None and override_skip_incomplete:
                if eventNumber[i] not in self.events_with_waveforms.keys():
                    continue

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

            radiantThrs = None
            lowTrigThrs = None
            if self.__read_daq_status:
                # associate daq infomation of event based on readout times
                readout_time = readoutTime[i]
                radiant_idx = find_daq_status_index(readout_time, self._readout_time_radiant)
                lt_idx = find_daq_status_index(readout_time, self._readout_time_lt)
                radiantThrs = self._radiantThrs[radiant_idx]
                lowTrigThrs = self._lowTrigThrs[lt_idx]

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
                sampleRate = sampleRate[i],
                radiantThrs = radiantThrs,
                lowTrigThrs = lowTrigThrs,
                hasWaveforms = eventNumber[i] in self.events_with_waveforms.keys() if not self.skip_incomplete else True
            )

            infos.append(info)

        if not self.multiple:
            return info
        else:
            return infos


    def N(self) -> int:
        return self._hds.num_entries

    def _get_windows(self, kw):
        """ Helper to access uproot file """
        return self._hds[f'trigger_info/trigger_info.radiant_info.start_windows[{self.NUM_CHANNELS}][2]'].array(**kw)

    def _get_waveforms(self, kw):
        """ Helper to access uproot file """
        return self._wfs[f'radiant_data[{self.NUM_CHANNELS}][{self.NUM_WF_SAMPLES}]'].array(**kw)

    def wfs(self, calibrated : bool = False, channels : Optional[Union[int, List[int]]] = None, override_skip_incomplete : Optional[bool] = None) -> Optional[numpy.ndarray]:
        """
        Returns the waveform data for the selected event(s).

        Parameters
        ----------

        calibrated : bool (default: False)
            If True, return the calibrated waveform data in units of volt. If False, return the raw waveform data in units of ADC counts.

        channels : int or list(int) (default: None)
            If None, return all channels. If int, return only the specified channel. If list(int), return only the specified channels.
            This is in particular useful if you are only interested in a subset of calibrated channels.

        override_skip_incomplete : bool (default: None)
            If not None override behaviour set in the contructor.

        Returns
        -------

        wfs : numpy.ndarray
            An array containing the waveform data with the shape (n_events, n_channels, n_samples).
        """

        if calibrated and not self.has_calib:
            raise ValueError("You requested a calibrated waveform but no calibration is available")

        # assert(not calibrated) # not implemented yet
        kw = dict(entry_start=self.first, entry_stop=self.last, library='np')

        w = None
        if self.full or self.skip_incomplete:
            w = self._get_waveforms(kw)
            starting_window = self._get_windows(kw)
        elif not self.multiple:
            # if you only selected one event and have an incomplete dataset
            if self.first in self.events_with_waveforms:
                idx = self.events_with_waveforms[self.first]
                w = self._get_waveforms(dict(entry_start=idx, entry_stop=idx+1, library='np'))
                starting_window = self._get_windows(dict(entry_start=idx, entry_stop=idx+1, library='np'))
        else:
            # so ... we need to loop through and find which things we have actually have waveforms
            # now figure out how much of the data array we need
            wf_start = None
            wf_end = None

            # store the indices that will be non-zero
            wf_idxs = []
            for idx in range(self.first, self.last):
                if idx in self.events_with_waveforms:
                    wf_idxs.append(idx - self.first)
                    # these are the start and stop of our array we need to load
                    if wf_start is None:
                        wf_start = self.events_with_waveforms[idx]

                    wf_end = self.events_with_waveforms[idx] + 1

            if len(wf_idxs):
                # take the overwrite if it is not None
                skip_incomplete = override_skip_incomplete or self.skip_incomplete
                if skip_incomplete:
                    w = self._get_waveforms(dict(entry_start=wf_start, entry_stop=wf_end, library='np'))
                    starting_window = self._get_windows(dict(entry_start=wf_start, entry_stop=wf_end, library='np'))
                else:
                    w = numpy.empty((self.last - self.first, self.NUM_CHANNELS, self.NUM_WF_SAMPLES), dtype='float64' if calibrated else 'int16')
                    starting_window = numpy.empty((self.last - self.first, self.NUM_CHANNELS, 2))
                    w[wf_idxs] = self._get_waveforms(dict(entry_start=wf_start, entry_stop=wf_end, library='np'))
                    starting_window[wf_idxs] = self._get_windows(dict(entry_start=wf_start, entry_stop=wf_end, library='np'))

        if channels is not None:
            if isinstance(channels, int):
                channels = [channels]

            w = w[:, channels]
            starting_window = starting_window[:, channels]

        if w is not None:
            # calibration
            starting_window = starting_window[:, :, 0]
            if calibrated:
                # this can run now both normal and raw calibration
                w = numpy.array([self.vc(ele, starting_window[i]) for i, ele in enumerate(w)])

            w = numpy.asarray(w, dtype=float)

            if self.multiple:
                return w

        return None if w is None else w[0]

    def _iterate(
            self, start : int, stop : int, calibrated: bool,  max_in_mem : int,
            selectors: Optional[Union[Callable[[mattak.Dataset.EventInfo], bool],
            Sequence[Callable[[mattak.Dataset.EventInfo], bool]]]] = None,
            override_skip_incomplete : Optional[bool] = None) \
            -> Generator[Tuple[Optional[mattak.Dataset.EventInfo], Optional[numpy.ndarray]], None, None]:

        # cache current values given by setEntries(..)
        original_entry : Union[int, Tuple[int, int]] = (self.first, self.last) if self.multiple else self.entry

        # determine in how many batches we want to access the data given how much events we want to load into the RAM at once
        n_batches = math.ceil((stop - start) / max_in_mem)

        if not isinstance(selectors, (list, numpy.ndarray)) and selectors is not None:
            selectors = [selectors]

        for i_batch in range(n_batches):

            # looping over the batches defining the start and stop index
            batch_start = start + i_batch * max_in_mem
            batch_stop = min(stop, batch_start + max_in_mem)
            self.setEntries((batch_start, batch_stop))

            # load events from file
            wfs = self.wfs(calibrated, override_skip_incomplete=override_skip_incomplete)
            es = self.eventInfo(override_skip_incomplete=override_skip_incomplete)

            # we modified the internal data pointers with the prev. call of self.setEntries(...)
            # this is intransparent for the outside world and has to be reverted
            self.setEntries(original_entry)

            # can happen if we skip incomplete events
            if wfs is None:
                continue

            for e, w in zip(es, wfs):
                if selectors is not None:
                    if numpy.all([selector(e) for selector in selectors]):
                        yield e, w
                else:
                    yield e, w
