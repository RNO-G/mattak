### Python dataset class, agnostic to backend


import os 
from abc import ABC, abstractmethod
from dataclasses import dataclass
import typing
from typing import Sequence, Union, Tuple, Optional
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
    sysclkLastPPS: Tuple[int,int]  # the last 2 PPS sysclks, most recent first
    pps: int 
    radiantStartWindows: numpy.ndarray
    sampleRate: float  # Sample rate, in GSa/s 

##
class AbstractDataset(ABC): 
    """ Abstract Base Class for accessing RNO-G data in Python , implemented either with uproot or PyROOT"""

    def setEntries( self, i : Union[int,Tuple[int,int]]):
        """
        Select entries to read out with wfs or eventInfo. Can either be a
        single entry or a tuple representating start and end entries, exclusive.
        Use (0,dataset.N()) or (0, None) to select all events. Negative indices indicate from end. 
        """ 

        if isinstance(i, tuple):
            self.multiple  = True
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
            self.last = i+1
 

    def N(self) -> int: 
        """ Return the number of events available in this dataset"""
        return 0


    @abstractmethod
    def _iterate(self, start: int , stop : Union[int,None] , calibrated: bool, max_entries_in_mem: int)-> Tuple[EventInfo, numpy.ndarray]:
        """ implementation-defined part of iterator"""
        pass

    def iterate(self, start : int = 0, stop : Union[int,None] = None,  calibrated: bool = False, max_entries_in_mem : int = 256) -> Optional[Tuple[EventInfo, numpy.ndarray]]:
        """ Iterate over events from start to stop, holding at most max_entries_in_mem in RAM.
            Returns a tuple of EventInfo and the event waveforms (potentially calibrated). 
        """
        if start < 0: 
            start += self.N() 
        if start < 0 or start > self.N(): 
            return None

        if stop is None: 
            stop = self.N()
        
        if stop < 0: 
            stop += self.N() 

        if stop < 0 or start > self.N(): 
            return None

        return self._iterate(start,stop,calibrated,max_entries_in_mem) 



    @abstractmethod
    def eventInfo(self) -> Union[Optional[EventInfo],Optional[Sequence[EventInfo]]]:
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


def Dataset(station, run, data_dir = None, backend="auto", verbose = False): 
   """
   This is not a class, but a factory method! 
   Returns a dataset corresponding to the station and run using data_dir as the base. If data_dir is not defined,
   then the environmental variable RNO_G_DATA will be used. The backend can be chosen explicitly or auto will try to
   use the best one. 
   """

   if data_dir is None: 
       data_dir = os.environ['RNO_G_DATA'] 

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
        return mattak.backends.uproot.dataset.Dataset(station, run, data_dir)
   elif backend == "pyroot": 
        import mattak.backends.pyroot.dataset 
        return mattak.backends.pyroot.dataset.Dataset(station, run, data_dir) 
   else: 
       print("Unknown backend (known backends are \"uproot\" and \"pyroot\")")
       return None 

