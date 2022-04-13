import ROOT 
import mattak.backends.pyroot.mattakloader
import mattak.Dataset
import typing
import cppyy.ll
import numpy

class Dataset(mattak.Dataset.AbstractDataset): 

    def __init__(self, station : int, run : int, data_dir : str): 
        self.ds = ROOT.mattak.Dataset(station, run, ROOT.nullptr, data_dir) 
        self.entries =0 
        self.station = station
        self.run = run 

    def N(self) -> int: 
        return self.ds.N() 

    def setEntries(self, i : typing.Union[int,typing.Sequence[int]]) -> None : 
            self.entries = i 
    
    
    def _eventInfo(self, i : int) -> mattak.Dataset.EventInfo:
        #todo: handle this in C++ code if it's too slow in Python
        self.ds.setEntry(i)
        hdr = self.ds.header() 
        if hdr is None: 
            return None

        assert( hdr.station_number == self.station)
        assert (hdr.run_number == self.run) 
        station = hdr.station_number
        run = hdr.run_number
        eventNumber = hdr.event_number
        readoutTime = hdr.readout_time
        triggerTime = hdr.trigger_time
        triggerType = "UNKNOWN"
        if (hdr.trigger_info.radiant_trigger):
            which = hdr.trigger_info.which_radiant_trigger
            if which == -1:
                which = "X" 
            triggerType = "RADIANT" + which
        elif hdr.trigger_info.lt_trigger:
            triggerType = "LT" 
        elif hdr.trigger_info.force_trigger:
            triggerType = "FORCE" 
        elif hdr.trigger_info.pps_trigger:
            triggerType = "PPS" 

        pps = hdr.pps_num
        sysclk = hdr.sysclk
        sysclkLastPPS = ( hdr.sysclk_last_pps, hdr.sysclk_last_last_pps) 
        radiantStartWindows = numpy.frombuffer( cppyy.ll.cast['uint8_t*'](hdr.trigger_info.radiant_info.start_windows), dtype='uint8', count=24*2).reshape(24,2)
        return mattak.Dataset.EventInfo(eventNumber = eventNumber, station = station, run = run, readoutTime=readoutTime, triggerTime=triggerTime, triggerType=triggerType, sysclk=sysclk, sysclkLastPPS=sysclkLastPPS, pps=pps, radiantStartWindows = radiantStartWindows)


    def eventInfo(self) -> typing.Union[mattak.Dataset.EventInfo,typing.Sequence[mattak.Dataset.EventInfo]]: 

        if isinstance(self.entries, typing.Sequence): 
            infos = [] 
            for i in self.entries: 
                infos.append(self._eventInfo(i))
            return infos 

        return self._eventInfo(self.entries)

    def _wfs(self, i : int, calibrated: bool = False): 
        self.ds.setEntry(i)
        wf = self.ds.calibrated() if calibrated else self.ds.raw()
        if wf is None: 
            return None 
        return numpy.frombuffer(cppyy.ll.cast['double*' if calibrated else 'int16_t*'](wf.radiant_data), dtype = 'float64' if calibrated else 'int16', count = 24*2048).reshape(24,2048)


    def wfs(self, calibrated : bool =False) -> numpy.ndarray: 

            
        # the simple case first
        if not isinstance(self.entries, typing.Sequence):
            return self._wfs(self.entries, calibrated) 
        
        #otherwise, let's set up the output ahead of time 

        out = numpy.zeros((len(self.entries), 24, 2048))
        for i, entry in enumerate(self.entries):
            out[i][:][:] = self._wfs(entry, calibrated)

            
        return out

        





    




