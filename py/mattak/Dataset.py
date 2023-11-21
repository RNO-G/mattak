### Python dataset class, agnostic to backend


import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
import typing
from typing import Sequence, Union, Tuple, Optional, Generator, Callable
import numpy
import datetime


@dataclass
class EventInfo:
    """ Pure python event information. In effect duplicating the most important bits of the ROOT header"""
    eventNumber: int
    station : int
    run: int
    readoutTime : float
    triggerTime : float
    triggerType: str
    sysclk: int
    sysclkLastPPS: Tuple[int, int]  # the last 2 PPS sysclks, most recent first
    pps: int
    radiantStartWindows: numpy.ndarray
    sampleRate: float  # Sample rate, in GSa/s
    radiantThrs: numpy.ndarray
    lowTrigThrs: numpy.ndarray

##
class AbstractDataset(ABC):
    """ Abstract Base Class for accessing RNO-G data in Python , implemented either with uproot or PyROOT"""

    def setEntries(self, i : Union[int, Tuple[int, int]]):
        """
        Select entries to read out with wfs or eventInfo. Can either be a
        single entry or a tuple representating start and end entries, exclusive.
        Use (0,dataset.N()) or (0, None) to select all events. Negative indices indicate from end.
        """

        if isinstance(i, tuple):
            self.multiple = True
            self.first = i[0]
            self.last = i[1]
            if self.last is None:
                self.last = self.N()
            if self.first < 0:
                self.first += self.N()
            if self.last < 0:
                self.last += self.N()

            if self.last > self.N():
                self.last = self.N()
        else:
            self.multiple = False
            self.entry = i
            self.first = i
            self.last = i + 1


    def N(self) -> int:
        """ Return the number of events available in this dataset"""
        return 0


    @abstractmethod
    def _iterate(self, start: int , stop : int , calibrated: bool, max_entries_in_mem: int) -> Generator[Tuple[EventInfo, numpy.ndarray],None,None]:
        """ implementation-defined part of iterator"""
        pass

    def iterate(self, start : int = 0, stop : Union[int,None] = None,
                calibrated: bool = False, max_entries_in_mem : int = 256,
                selector: Optional[Callable[[EventInfo], bool]] = None) \
                -> Generator[Optional[Tuple[EventInfo, numpy.ndarray]], None, None]:
        """ Iterate over events from start to stop, holding at most max_entries_in_mem in RAM.
            Returns a tuple of EventInfo and the event waveforms (potentially calibrated).
        """
        if start < 0:
            start += self.N()
        if start < 0 or start > self.N():
            return

        if stop is None:
            stop = self.N()

        if stop < 0:
            stop += self.N()

        if stop < 0 or start > self.N():
            return

        yield from self._iterate(start, stop, calibrated, max_entries_in_mem, selector)

    @abstractmethod
    def eventInfo(self) -> Union[Optional[EventInfo], Sequence[Optional[EventInfo]]]:
        """Get selected event info(s)
           Depending on what was passed to setEntries this can return either one EventInfo or a list of them
        """
        pass

    @abstractmethod
    def wfs(self, calibrated : bool = False) -> Optional[numpy.ndarray]:
        """ Get select waveform(s).
            Depending on what was passed to setEntries, this may be a single waveform or many
        """
        pass
    
    def get_selected_wfs(self, selector: Callable[[EventInfo], bool], calibrated : bool = False) -> Optional[numpy.ndarray]:
        """ Convenience interface to use selector """
        return numpy.array([wf for _, wf in self.iterate(start=self.start, stop=self.stop, selector=selector)])


def Dataset(station : int, run : int, data_dir : str = None, backend : str= "auto", 
            verbose : bool = False, skip_incomplete : bool = True,
            read_daq_status : bool = True, read_run_info : bool = True
            preferred_file : Optional[str] = None ) -> Optional[AbstractDataset]:
   """

   This is not a class, but a factory method! This is meant to be the interface
   to rule them all for loading RNO-G data. Due to Cosmin's poor initial API
   design, it has become perhaps more complicated than it should be. 

   If data_dir is a directory and station/run are non-zero, this returns a
   dataset corresponding to the station and run using data_dir as the base
   directory (i.e. the folder hierarcy is structured something likelike
   ${data_dir}/stationX/runY/*.root). If data_dir is None, then the
   environmental variable RNO_G_DATA (or RNO_G_ROOT_DATA) will be queried and
   the run loaded from that base. data_dir can also be a URL for loading of
   files via HTTP (e.g. https://user:password@example.com/rno-g-data), though
   there may be some subtleties about escaping passwords that may differ
   betweeen different backends. 

   In the special case of setting station = 0 and run = 0, data_dir
   will be interpreted as a directory containing ROOT files, which is useful if
   you don't have the full directory hierarchy setup or want to look at data
   taken with the fakedaq.

   If data_dir is, despite its name, a file rather than directory, then that
   file will be attempted to be loaded as a combined file. Really it should
   be called data_source, changing the parameter name will break the API and
   we'll try hard not to do that. 

   The backend can be chosen explicitly ("pyroot" or "uproot") or auto will try
   to use the best one ("pyroot" if available, otherwise reverting to "uproot").

   verbose prints out things mostly useful for debugging.

   read_daq_status and read_run_info are mostly self-explanatory. Avoiding
   reading them may speed things up or work around the files not being there
   for some reason. 

   skip_incomplete affects what happens when a dataset is incomplete
   (telemetered). If True, will only index fully telemetered events. If False,
   will be indexed by full dataset, but untelemetered events will have
   waveforms of None type (if requested singly) or be all 0's (if requested
   via the bulk interface, as numpy doesn't support jagged ararys).

   preferred_file, if not None or "", and data_dir is not a file, will further
   change the loading behavior. By default,  we will try to load full waveforms
   falling back to loading combined.root. But if preferred_file is set, it will
   prefer loading ${preferred_file}.root if possible, treating it as a file in
   the same format as combined.root.  For example, you can set it to "combined"
   to load combined.root even if full waveforms are available. Or, if you have
   your own subselection in the same format (e.g. for example if you generated
   a file that is only forced triggers) this provides an arguably convenient
   way to load those. 


   """

   if data_dir is None:
       for env_var in ['RNO_G_DATA', 'RNO_G_ROOT_DATA']: 
           if env_var in os.environ: 
               data_dir = os.environ[env_var] 
               break 
       if data_dir is None:                
           print("Neither data_dir nor any relevant environmental variable (e.g. RNO_G_DATA) is defined and I don't know where else to look :(")
           return None 

   if backend == "auto":
        try:
            import ROOT
            import mattak.backends.pyroot.mattakloader
            if verbose:
                print('Using pyroot backend')
            backend = "pyroot"
        except:
            try:
                import uproot
                backend = "uproot"
                if verbose:
                    print('Using uproot backend')
            except:
                print("No backends available")
                return None

   if backend == "uproot":
        import mattak.backends.uproot.dataset
        return mattak.backends.uproot.dataset.Dataset(station, run, data_dir, verbose, skip_incomplete, read_daq_status, read_run_info, preferred_file)
   elif backend == "pyroot":
        import mattak.backends.pyroot.dataset
        return mattak.backends.pyroot.dataset.Dataset(station, run, data_dir, verbose, skip_incomplete, read_daq_status, read_run_info, preferred_file)
   else:
       print("Unknown backend (known backends are \"uproot\" and \"pyroot\")")
       return None
