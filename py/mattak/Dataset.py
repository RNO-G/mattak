### Python dataset class, agnostic to backend


import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
import typing
from typing import Sequence, Union, Tuple, Optional, Generator, Callable
import numpy
import datetime
import warnings


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
    sampleRate: Optional[float]  # Sample rate, in GSa/s
    radiantThrs: Optional[numpy.ndarray]
    lowTrigThrs: Optional[numpy.ndarray]

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
            if i < 0:
                i+= self.N() 
            if i is None: 
                i = 0 
            self.entry = i
            self.first = i
            self.last = i + 1


    def N(self) -> int:
        """ Return the number of events available in this dataset"""
        return 0


    @abstractmethod
    def _iterate(self, start: int , stop : int , calibrated: bool, max_entries_in_mem: int, selector: Optional[Callable[[EventInfo],bool]]) -> Generator[Tuple[Optional[EventInfo], Optional[numpy.ndarray]],None,None]:
        """ implementation-defined part of iterator"""
        pass

    def iterate(self, start : int = 0, stop : Union[int,None] = None,
                calibrated: bool = False, max_entries_in_mem : int = 256,
                selector: Optional[Callable[[EventInfo], bool]] = None) \
                -> Generator[Tuple[Optional[EventInfo], Optional[numpy.ndarray]], None, None]:
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
        return numpy.array([wf for _, wf in self.iterate(start=self.first, stop=self.last, selector=selector)])


def Dataset(station : int = 0, run : int = 0, data_path : Optional[str] = None, backend : str= "auto", 
            verbose : bool = False, skip_incomplete : bool = True,
            read_daq_status : bool = True, read_run_info : bool = True,
            preferred_file : Optional[str] = None, 
            *, data_dir : Optional[str] = None ) -> Optional[AbstractDataset]:
   """

   This is not a class, but a factory method! This is meant to be the interface
   to rule them all for loading RNO-G data. Due to Cosmin's poor initial API
   design, it has become perhaps more complicated than it should be. 

   In the case of setting station = 0 and run = 0, data_path will be
   interpreted as a directory containing ROOT files, which is useful if you
   don't have the full directory hierarchy setup or want to look at data taken
   with the fakedaq.

   If data_path is a file rather than directory, then that file will be
   attempted to be loaded as a combined file, ignoring the provided run and
   station numbers. 

   If data_path is a directory and station/run are non-zero, this returns a
   dataset corresponding to the station and run using data_path as the base
   directory (i.e. the folder hierarcy is structured something like
   ${data_path}/stationX/runY/*.root). 


   If data_path is None (or ""), then the environmental variable RNO_G_DATA
   (and also  RNO_G_ROOT_DATA) will be queried.  data_path can also be a URL
   for loading of files via HTTP (e.g.
   https://user:password@example.com/rno-g-data), though there may be some
   subtleties about escaping passwords that may differ betweeen different
   backends. 

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

   preferred_file, if not None or "", and data_path is not a file, will further
   change the loading behavior. By default,  we will try to load full waveforms
   falling back to loading combined.root. But if preferred_file is set, it will
   prefer loading ${preferred_file}.root if possible, treating it as a file in
   the same format as combined.root.  For example, you can set it to "combined"
   to load combined.root even if full waveforms are available. Or, if you have
   your own subselection in the same format (e.g. for example if you generated
   a file that is only forced triggers) this provides an arguably convenient
   way to load those. 

   """


   # handle deprecated name data_dir 
   if data_dir is not None: 
       warnings.warn("data_dir is deprecated, use data_path instead. This may be removed in the future, breaking your code.")
       if data_path is not None: 
           raise TypeError("Dataset received both data_path and data_dir!") 
       data_path = data_dir 


   #treat "" as an alias for None
   if data_path == "": 
       data_path = None


   if data_path is None:
       for env_var in ['RNO_G_DATA', 'RNO_G_ROOT_DATA']: 
           if env_var in os.environ: 
               data_path = os.environ[env_var] 
               break 
       if data_path is None:                
           print("Neither data_path nor any relevant environmental variable (e.g. RNO_G_DATA) is defined and I don't know where else to look :(")
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
        return mattak.backends.uproot.dataset.Dataset(station, run, data_path, verbose, skip_incomplete, read_daq_status, read_run_info, preferred_file)
   elif backend == "pyroot":
        import mattak.backends.pyroot.dataset
        return mattak.backends.pyroot.dataset.Dataset(station, run, data_path, verbose, skip_incomplete, read_daq_status, read_run_info, preferred_file)
   else:
       print("Unknown backend (known backends are \"uproot\" and \"pyroot\")")
       return None
