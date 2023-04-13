import ROOT
import mattak.backends.pyroot.mattakloader
import mattak.Dataset
from typing import Sequence, Union, Tuple, Optional, Callable
import numpy

try:
    import cppyy.ll

# work around a weird issue that happens on some systems
# for some reason, in some configurations it  doesn't detect free in the global namespace.
# So we'll define a different function in the global namespace with a super creative name
# that does the same thing as free then set it equal to free.
except AttributeError:
    import cppyy
    cppyy.cppdef("void freee(void *p) { free(p); }")
    cppyy.gbl.free = cppyy.gbl.freee
    import cppyy.ll


cppyy.cppdef(" bool is_nully(void *p) { return !p; }");


def isNully(p):
    return p is None or ROOT.AddressOf(p) == 0 or  cppyy.gbl.is_nully(p)

class Dataset(mattak.Dataset.AbstractDataset):

    def __init__(self, station : int, run : int, data_dir : str, verbose: bool = False, skip_incomplete: bool = True):
        
        self.backend = "pyroot"

        #special case where we load a directory instead of a station/run
        if station == 0 and run == 0:
            if verbose:
                print("Trying to load data dir!")
                
            self.ds = ROOT.mattak.Dataset()
            self.ds.setVerbose(verbose)

            self.ds.loadDir(data_dir, skip_incomplete)
            self.station = self.ds.header().station_number
            self.run = self.ds.header().run_number

            if verbose:
                print("We think we found station %d run %d" % (self.station,self.run))

        else:
            self.ds = ROOT.mattak.Dataset(station, run, ROOT.nullptr, data_dir, skip_incomplete, verbose)
            self.station = station
            self.run = run
                    
        self.data_dir = data_dir
        self.setEntries(0)

    def N(self) -> int:
        return self.ds.N()

    def _eventInfo(self, i : int) -> Optional[mattak.Dataset.EventInfo]:
        #TODO: handle this in C++ code if it's too slow in Python
        if not self.ds.setEntry(i):
            return None
        
        hdr = self.ds.header()
        
        daq_status = self.ds.status()
        radiantThrs = numpy.array(daq_status.radiant_thresholds)
        lowTrigThrs = numpy.array(daq_status.lt_trigger_thresholds)

        assert(hdr.station_number == self.station)
        assert(hdr.run_number == self.run)
        station = hdr.station_number
        run = hdr.run_number
        eventNumber = hdr.event_number
        readoutTime = hdr.readout_time
        triggerTime = hdr.trigger_time
        triggerType = "UNKNOWN"
        
        if hdr.trigger_info.radiant_trigger:
            which = hdr.trigger_info.which_radiant_trigger
            if which == -1:
                which = "X"
            triggerType = "RADIANT" + str(which)
        elif hdr.trigger_info.lt_trigger:
            triggerType = "LT"
        elif hdr.trigger_info.force_trigger:
            triggerType = "FORCE"
        elif hdr.trigger_info.pps_trigger:
            triggerType = "PPS"

        pps = hdr.pps_num
        sysclk = hdr.sysclk
        sysclkLastPPS = (hdr.sysclk_last_pps, hdr.sysclk_last_last_pps)
        radiantStartWindows = numpy.frombuffer(cppyy.ll.cast['uint8_t*'](hdr.trigger_info.radiant_info.start_windows), dtype='uint8', count=24 * 2).reshape(24, 2)
        sampleRate = 3.2  #if ( ROOT.AddressOf(self.ds.info()) ==0 or self.ds.info().radiant_sample_rate == 0)  else self.ds.info().radiant_sample_rate/1000.


        return mattak.Dataset.EventInfo(eventNumber = eventNumber,
                                        station = station,
                                        run = run,
                                        readoutTime=readoutTime,
                                        triggerTime=triggerTime,
                                        triggerType=triggerType,
                                        sysclk=sysclk,
                                        sysclkLastPPS=sysclkLastPPS,
                                        pps=pps,
                                        radiantStartWindows = radiantStartWindows,
                                        sampleRate = sampleRate,
                                        radiantThrs=radiantThrs,
                                        lowTrigThrs=lowTrigThrs)


    def eventInfo(self) -> Union[Optional[mattak.Dataset.EventInfo],Sequence[Optional[mattak.Dataset.EventInfo]]]:

        if self.multiple:
            infos = []
            for i in range(self.first, self.last):
                infos.append(self._eventInfo(i))
            
            return infos

        return self._eventInfo(self.entry)

    def _wfs(self, i : int, calibrated: bool = False):
        self.ds.setEntry(i)
        wf = self.ds.calibrated() if calibrated else self.ds.raw()
        if isNully(wf):
            return None
        
        return numpy.frombuffer(cppyy.ll.cast['double*' if calibrated else 'int16_t*'](wf.radiant_data), dtype = 'float64' if calibrated else 'int16', count=24 * 2048).reshape(24,2048)


    def wfs(self, calibrated : bool=False) -> numpy.ndarray:

        # the simple case first
        if not self.multiple:
            return self._wfs(self.entry, calibrated)

        if self.last - self.first < 0:
            return None

        out = numpy.zeros((self.last - self.first, 24, 2048), dtype='float64' if calibrated else 'int16')
        for entry in range(self.first, self.last):
            this_wfs = self._wfs(entry, calibrated)
            if this_wfs is not None:
                out[entry-self.first][:][:] = this_wfs

        return out

    # ignore max_in_mem since it doesn't save much time for us...
    def _iterate(self, start, stop, calibrated, max_in_mem,
                 selector: Optional[Callable[[mattak.Dataset.EventInfo], bool]] = None) -> Tuple[mattak.Dataset.EventInfo, numpy.ndarray]:

        if selector is not None:
            for i in range(start, stop):
                evinfo = self._eventInfo(i)
                if selector(evinfo):
                    yield evinfo, self._wfs(i, calibrated)
        else:
            for i in range(start, stop):
                yield self._eventInfo(i), self._wfs(i, calibrated)
        return










