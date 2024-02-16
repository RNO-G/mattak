import os 
import uproot
import errno 
import os.path
import configparser
import mattak.Dataset
import typing
from typing import Union, List, Optional, Tuple, Generator, Callable, Sequence
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

    def __init__(self, station : int, run : int, data_path : str, verbose : bool = False,
                 skip_incomplete : bool = True, read_daq_status : bool = True,
                 read_run_info : bool = True, file_preference: Optional[str] = None):
        
        self.backend = "uproot"
        self.__verbose = verbose
        self.__read_daq_status = read_daq_status
        self.__read_run_info = read_run_info

        #special case where data_path is a file

        if os.path.isfile(data_path): 
            self.data_path_is_file = True
            self.rundir = data_path 
        else: 
            self.data_path_is_file = False
            # special case where we load a directory instead of a station/run
            if station == 0 and run == 0: 
                self.rundir = data_path
            else: 
                self.rundir = "%s/station%d/run%d" % (data_path,station,run)



        if (skip_incomplete == False  and self.data_path_is_file):
            print("skip_incomplete = false is incompatible with data_path as file"); 
            skip_incomplete = True 

        self.skip_incomplete = skip_incomplete 

        # this duplicates a bunch of C++ code in mattak::Dataset
        # check for full or partial run by looking for waveforms.root
        # Only open files/trees, do not access data

        self.combined_tree = None 

        # check for a preference
        if (file_preference is not None and file_preference != "") or self.data_path_is_file: 
            preferred_file = "%s:combined" % (self.rundir) if self.data_path_is_file  else "%s/%s.root:combined" %(self.rundir, file_preference) 
            try: 
                self.combined_tree = uproot.open(preferred_file)
                if verbose: 
                    print ("Found preferred file %s" % (preferred_file))
                    print (self.combined_tree)
                self.full = False
            except Exception:
                # can't recover from this 
                if self.data_path_is_file: 
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.rundir)
                print ("Could not find preferred file %s. Falling back to normal behavior" % (preferred_file))


        

        #if we didn't load the combined_tree already , try to load full tree 
        if self.combined_tree is None:
            try:
                self.wf_file = uproot.open("%s/waveforms.root" % (self.rundir))
                if self.__verbose:
                    print ("Open waveforms.root (Found full run folder) ...") 

                self.combined_tree = None
                self.full = True

                self.wf_tree, self.wf_branch = read_tree(self.wf_file, waveform_tree_names)
                self._wfs = self.wf_tree[self.wf_branch]

                self.hd_file = uproot.open("%s/headers.root" % (self.rundir))
                self.hd_tree,self.hd_branch = read_tree(self.hd_file, header_tree_names); 
                self._hds = self.hd_tree[self.hd_branch]

                if self.__read_daq_status:
                    self.ds_file = uproot.open("%s/daqstatus.root" % (self.rundir))
                    self.ds_tree, self.ds_branch = read_tree(self.ds_file, daqstatus_tree_names)
                    self._dss = self.ds_tree[self.ds_branch]
            except Exception: 
                self.full = False

        # we haven't already loaded the full tree
        if not self.full:   
            if self.combined_tree is None: # we didn't already load our preference
                self.combined_tree = uproot.open("%s/combined.root:combined" %(self.rundir))
                if self.__verbose: 
                    print ("Found combined file")
                    print (self.combined_tree)

            self._wfs, self.wf_branch = read_tree(self.combined_tree, waveform_tree_names)

            # build an index of the waveforms we do have
            if not self.skip_incomplete: 
                wfs_included = self._wfs['event_number'].array() 
                self.events_with_waveforms = {ev: idx for idx, ev in enumerate(wfs_included)}

                self.full_head_file = uproot.open("%s/headers.root" % (self.rundir))
                self.full_head_tree,_ = read_tree(self.full_head_file, header_tree_names)
                
                if self.__read_daq_status:
                    self.full_daq_file = uproot.open("%s/daqstatus.root" % (self.rundir))
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
        if self.__read_run_info:
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

    def wfs(self, calibrated : bool =False) -> Optional[numpy.ndarray]: 
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

        # here we'd eventually handle calibration I think? 

        if self.multiple: 
            return w
        return None if w is None else w[0] 

    def _iterate(self, start : int , stop : int , calibrated: bool,  max_in_mem : int, 
                 selector: Optional[Callable[[mattak.Dataset.EventInfo],bool]] = None) -> Generator[Tuple[Optional[mattak.Dataset.EventInfo], Optional[numpy.ndarray]],None,None]:

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

