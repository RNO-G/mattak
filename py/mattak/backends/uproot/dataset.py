import os 
import uproot
import os.path
import configparser
import glob

import mattak.Dataset
import typing
from typing import Union, Tuple, Optional, Callable, Sequence
import numpy
import math

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

    def __init__(self, station : int, run : int, data_dir_file : str, verbose : bool = False,
                 skip_incomplete : bool = True, read_daq_status : bool = True,
                 read_run_info : bool = True):
        """
        Uproot backend for the python interface of the mattak Dataset
        
        Parameters
        ----------
        
        station: int
            Station Id of the data to be read. 
            
        run: int
            Run number of the data to be read
            
        data_dir_file: str
            Path of the data directory or file
            
        verbose: bool
            Enable verbose debug log
            
        skip_incomplete: bool
            Default: True
        
        read_daq_status: bool 
            Enable reading of daqstatus.root. Only possible  when a complete run directory is passed. (Default: True)
        
        read_run_info: bool
            Enable reading of the run information (sampling rate). Only possible when a run directory is passed 
            which contains the runinfo.txt file. (Default: True, unless you pass a specific file as input)
        """
        
        self.backend = "uproot"
        self.__verbose = verbose
        self.__read_daq_status = read_daq_status
        self.__read_run_info = read_run_info
        
        # default: consider each run as complete
        self.full = True
        self.runfile = None
        # special case where we load a directory instead of a station/run
        if station == 0 and run == 0:
            if os.path.isfile(data_dir_file):
                self.runfile = data_dir_file
                self.full = False
            else:
                self.rundir = data_dir_file
        else: 
            self.rundir = "%s/station%d/run%d" % (data_dir_file, station, run)

        self.skip_incomplete = skip_incomplete 

        # this duplicates a bunch of C++ code in mattak::Dataset
        # check for full or partial run by looking for waveforms.root
        # Only open files/trees, do not access data
        if self.full and os.path.exists("%s/waveforms.root" % (self.rundir)):

            self.wf_file = uproot.open("%s/waveforms.root" % (self.rundir))
            if self.__verbose: 
                print("Open waveforms.root (Found full run folder) ...")
            
            self.combined_tree = None
            
            self.wf_tree, self.wf_branch = read_tree(self.wf_file, waveform_tree_names)
            self._wfs = self.wf_tree[self.wf_branch]
            
            self.hd_file = uproot.open("%s/headers.root" % (self.rundir))
            self.hd_tree, self.hd_branch = read_tree(self.hd_file, header_tree_names)
            self._hds = self.hd_tree[self.hd_branch]

            # try finding a calibration file in the run directory
            try:
                self.cal_file = uproot.open(glob.glob(f"{self.rundir}/volCalConst*.root")[0])
                if self.__verbose:
                    print("Opening calibration file")
            except:  
                self.cal_file = None
                if self.__verbose:
                    print("No calibration file found, continuing without calibration")
            
            if self.__read_daq_status:
                self.ds_file = uproot.open("%s/daqstatus.root" % (self.rundir))
                self.ds_tree, self.ds_branch = read_tree(self.ds_file, daqstatus_tree_names)
                self._dss = self.ds_tree[self.ds_branch]

        elif self.runfile is not None or os.path.exists("%s/combined.root" % (self.rundir)):
            self.full = False
            
            if self.runfile is not None:
                self.combined_tree = uproot.open("%s:combined" %(self.runfile))
            else:
                self.combined_tree = uproot.open("%s/combined.root:combined" %(self.rundir))
            
            if self.__verbose: 
                print("Open combined.root ...")
                
            self._wfs, self.wf_branch = read_tree(self.combined_tree, waveform_tree_names)

            if not skip_incomplete:
                
                # build an index of the waveforms we do have
                wfs_included = self._wfs['event_number'].array() 
                self.events_with_waveforms = {ev: idx for idx, ev in enumerate(wfs_included)}
                
                # Open header and daqstatus files to get information of all events!
                self.full_head_file = uproot.open("%s/headers.root" % (self.rundir))
                self.full_head_tree, _ = read_tree(self.full_head_file, header_tree_names)
                
                if self.__read_daq_status:
                    self.full_daq_file = uproot.open("%s/daqstatus.root" % (self.rundir))
                    self.full_daq_tree, _ = read_tree(self.full_daq_file, daqstatus_tree_names)
                    
            # Get header and daq information from combined file or from full run
            hd_tree = self.combined_tree if skip_incomplete else self.full_head_tree
            self._hds, self.hd_branch = read_tree(hd_tree, header_tree_names)
                
            if self.__read_daq_status:
                ds_tree = self.combined_tree if skip_incomplete else self.full_daq_tree
                self._dss, self.ds_branch = read_tree(ds_tree, daqstatus_tree_names)
        
        else:
            raise ValueError("No vaild run directory / file found!")

        if station == 0 and run == 0: 
            self.station = self._hds['station_number'].array(entry_start=0, entry_stop=1)[0]
            self.run = self._hds['run_number'].array(entry_start=0, entry_stop=1)[0]

        else:
            self.station = station
            self.run = run

        self.data_dir = data_dir_file
        self.setEntries(0) 

        self.run_info = None
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
                

    def eventInfo(self) -> Union[Optional[mattak.Dataset.EventInfo],Sequence[Optional[mattak.Dataset.EventInfo]]]: 
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
        radiantStartWindows = self._hds['trigger_info/trigger_info.radiant_info.start_windows[24][2]'].array(**kw, library='np')
        
        infos = [] 
        info = None  # if range(0)
        for i in range(self.last-self.first): 
            
            triggerType  = "UNKNOWN"; 
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

            info = mattak.Dataset.EventInfo(eventNumber = eventNumber[i], 
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
                                            lowTrigThrs=lowTrigThrs[i] if self.__read_daq_status else None)
            
            infos.append(info)
            
        if not self.multiple:
            return info
        else:
            return infos
    
    
    def N(self) -> int: 
        return self._hds.num_entries
    
    def __unpack_cal_parameters(self, cal_file : uproot.ReadOnlyDirectory) -> numpy.ndarray:
        """
        Function that reads out 9th order polynomial parameters from the calibration file

        Returns
        ---------
        coef : np.ndarray of shape (24 * 4096, 10)
        """
        coef = numpy.stack(cal_file["coeffs_tree/coeff"].array(library = 'np'))   #stack is needed to convert ndarray of ndarrays to 'normally shaped' array, otherwise you have nested ndarrays
        return coef

    def __unpack_cal_residuals(self, cal_file : uproot.ReadOnlyDirectory) -> numpy.ndarray:
        """
        Function that reads out the residuals from the calibration file

        Returns
        ---------
        (v_residuals, residuals) : both v_residuals and residuals have shape (points, 2)

        """

        vres_dac1, vres_dac2 = cal_file["aveResid_dac1"].values(axis = 0), cal_file["aveResid_dac2"].values(axis = 0)
        residual_dac1, residual_dac2 = cal_file["aveResid_dac1"].values(axis = 1), cal_file["aveResid_dac2"].values(axis = 1)
        return numpy.stack(numpy.array([vres_dac1, vres_dac2]), axis = -1), numpy.stack(numpy.array([residual_dac1, residual_dac2]), axis = -1)
    
    def __calibrate(self, waveform_array : numpy.ndarray, param : numpy.ndarray, 
                    vres : numpy.ndarray, res : numpy.ndarray, starting_window : float | int,
                    fit_min : float = -1.3, fit_max : float = 0.7, accuracy : float = 0.01) -> numpy.ndarray:
        """
        The calibration function that transforms waveforms from ADC to voltage

        Parameters
        ------------
        waveform_array: array of one waveform, expected to have shape (24, 2048)
        param : the parameters found in the calibration file, expected to have shape (24 * 4096, 10)
        vres : the voltage points of the residuals
        res : the ADC values of the residuals
        starting_window : the sample on which the run started
        
        """

        # "discrete" inverse
        vsamples = numpy.arange(fit_min, fit_max, accuracy)
        waveform_volt = numpy.zeros((24, 2048))
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
                adcsamples = numpy.polyval(p[::-1], vsamples)
                volt = numpy.interp(adc, adcsamples, vsamples, left = fit_min, right = fit_max)
                waveform_volt[c, s] = volt
        return waveform_volt

    def wfs(self, calibrated : bool = False) -> numpy.ndarray: 
        # assert(not calibrated) # not implemented yet 

        w = None 
        if self.full or self.skip_incomplete: 
            w = self._wfs['radiant_data[24][2048]'].array(entry_start=self.first, entry_stop=self.last, library='np')
        elif not self.multiple: 
            if self.first in self.events_with_waveforms: 
                idx = self.events_with_waveforms[self.first]
                w = self._wfs['radiant_data[24][2048]'].array(entry_start=idx, entry_stop=idx+1, library='np')
        else: 
            # so ... we need to loop through and find which things we have actually have waveforms
            # start by allocating the output
            w = numpy.zeros((self.last - self.first, 24, 2048), dtype='float64' if calibrated else 'int16')

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

        # calibration
        if calibrated and self.cal_file:
            cal_param = self.__unpack_cal_parameters(self.cal_file)
            cal_residuals_v, cal_residuals_adc = self.__unpack_cal_residuals(self.cal_file)
            kw = dict(entry_start = self.first, entry_stop = self.last)
            radiantStartWindows = self._hds['trigger_info/trigger_info.radiant_info.start_windows[24][2]'].array(**kw, library='np')
            starting_window = radiantStartWindows[0, :, 0]
            w = numpy.array([self.__calibrate(ele, cal_param, cal_residuals_v, cal_residuals_adc,  starting_window) for ele in w])
        elif calibrated and not self.cal_file:
            print(f"No calibration file was found in {self.rundir}, calibration was not applied")

        if self.multiple:
            return w
        
        return None if w is None else w[0] 
    
    
    def _iterate(self, start, stop, calibrated, max_in_mem,
                 selector: Optional[Callable[[mattak.Dataset.EventInfo], bool]] = None) -> Tuple[mattak.Dataset.EventInfo, numpy.ndarray]:

        # cache current values given by setEntries(..)
        original_entry = (self.first, self.last) if self.multiple else self.entry

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
            







