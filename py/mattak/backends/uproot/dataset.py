import os 
import uproot

import mattak.Dataset
import typing
from typing import Union,List
import numpy


waveform_tree_names = ["waveforms","wfs","wf","waveform"]
header_tree_names = ["hdr","header","hd","hds","headers"] 

class Dataset ( mattak.Dataset.AbstractDataset):

    def __init__(self, station : int, run : int, data_dir : str): 
        self.rundir = "%s/station%d/run%d" % (data_dir,station,run)

        # this duplicates a bunch of C++ code in mattak::Dataset

        # check for full or partial run by looking for waveforms.root


        try: 
            self.wf_file = uproot.open("%s/waveforms.root" % (self.rundir))
            self.combined_tree = None
            for wf_tree_name in waveform_tree_names:
                if wf_tree_name in self.wf_file:
                    self.wf_tree = self.wf_file[wf_tree_name] 
                    self.wf_branch = wf_tree_name
                    self._wfs = self.wf_tree[self.wf_branch]
                    break 
            self.hd_file = uproot.open("%s/headers.root" % (self.rundir))
            for hd_tree_name in header_tree_names:
                if hd_tree_name in self.hd_file:
                    self.hd_tree = self.hd_file[hd_tree_name] 
                    self.hd_branch = hd_tree_name
                    self._hds = self.hd_tree[self.hd_branch]
                    break 
 
            self.full = True

        except: 
            self.full = False
            self.combined_tree = uproot.open("%s/combined.root;combined" %(self.rundir))
            # get the right branch names
            for wf_branch_name in waveform_tree_names:
                if wf_branch_name in self.combined_tree:
                    self.wf_branch = wf_branch_name
                    self._wfs = self.combined_tree[self.wf_branch]
                    break 
            for hd_branch_name in header_tree_names:
                if hd_branch_name in self.combined_tree:
                    self.hd_branch = hd_branch_name
                    self._hds = self.combined_tree[self.hd_branch]
                    break 

        self.station = station
        self.run = run
        self.data_dir = data_dir
        self.setEntries(0) 

        self.run_info = None
        # try to get the run info, if we're using combined tree, try looking in tehre 
        if self.combined_Tree is not None: 
            self.run_info = uproot.open("%s/combined.root;info", %(self.rundir))
        else 
            self.run_info = uproot.open("%s/runinfo.root;info", %(self.rundir))


    def eventInfo(self) -> Union[mattak.Dataset.EventInfo,typing.Sequence[mattak.Dataset.EventInfo]]: 
        station = self._hds['station_number'].array(entry_start = self.first, entry_stop = self.last)
        run = self._hds['run_number'].array(entry_start = self.first, entry_stop = self.last)
        eventNumber = self._hds['event_number'].array(entry_start = self.first, entry_stop = self.last)
        readoutTime = self._hds['readout_time'].array(entry_start = self.first, entry_stop = self.last)
        triggerTime = self._hds['trigger_time'].array(entry_start = self.first, entry_stop = self.last)
        triggerInfo = self._hds['trigger_info'].array(entry_start = self.first, entry_stop = self.last)
        pps = self._hds['pps_num'].array(entry_start = self.first, entry_stop = self.last)
        sysclk = self._hds['sysclk'].array(entry_start = self.first, entry_stop = self.last)
        sysclk_lastpps = self._hds['sysclk_last_pps'].array(entry_start = self.first, entry_stop = self.last)
        sysclk_lastlastpps = self._hds['sysclk_last_last_pps'].array(entry_start = self.first, entry_stop = self.last)
        sampleRate = 3.2 if self.run_info is None else self.run_info.radiant_sample_rate 

        # um... yeah, that's obvious 
        radiantStartWindows = self._hds['trigger_info/trigger_info.radiant_info.start_windows[24][2]'].array(entry_start = self.first, entry_stop = self.last, library='np')
        infos : Union[None,List[mattak.Dataset.EventInfo]] = [] if self.multiple else None
        for i in range(self.last-self.first): 
            triggerType  = "UNKNOWN"; 
            if triggerInfo[i]['trigger_info.radiant_trigger']:
                which = triggerInfo[i]['trigger_info.which_radiant_trigger']
                if which == -1: 
                    which = "X" 
                triggerType = "RADIANT" + which
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
                                            radiantStartWindows = radiantStartWindows[i])

            if not self.multiple:
                return info
            infos.append(info) 

        return infos



    

    def N(self) -> int: 
        return self._hds.num_entries

    def wfs(self, calibrated : bool =False) -> numpy.ndarray: 
#        assert(not calibrated) # not implemented yet 

        w = self._wfs['radiant_data[24][2048]'].array(entry_start = self.first, entry_stop = self.last, library='np')
        if self.multiple: 
            return w
        return w[0] 
        

    def _iterate(self, start, stop, calibrated,  max_in_mem) -> typing.Tuple[mattak.Dataset.EventInfo, numpy.ndarray]:

        current_start = -1 
        current_stop = -1
        w = None
        e = None
        i = start 
        preserve_entries = self.entry
        if self.multiple:
            preserve_entries = (self.first, self.last) 
        while i < stop: 
            if i >= current_stop:
                current_start = i
                current_stop = min(i+max_in_mem, stop)
                self.setEntries((current_start, current_stop)) 
                w = self.wfs(calibrated)
                e = self.eventInfo()
                self.setEntries(preserve_entries) # hide that we're modifying these 

            rel_i = i -current_start
            i+=1
            yield (e[rel_i], w[rel_i])

 








