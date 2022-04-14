### Python dataset class, agnostic to backend


import os 
from abc import ABC, abstractmethod
from dataclasses import dataclass
import typing
from typing import Sequence, Union, Tuple
import numpy 
import datetime 


## Pure python event information. In effect duplicating the ROOT header
## TODO: do we want this to be a struct of arrays rather than an array of structs? 
@dataclass 
class EventInfo: 
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

##
''' Abstract Base Class for accessing RNO-G data in Python , implemented either with uproot or PyROOT
  This presents an event-by-event interface. Bulk interface, uses the AbstractBulkDataset''' 

class AbstractDataset(ABC): 
    ''' Select entries to read out with wfs or eventInfo. Can either be a
    single entry or a tuple representating start and end entries, inclusive. Use (0,dataset.N()) to select all events.

    ''' 
    def setEntries( self, i : Union[int,Tuple[int,int]]):
        if isinstance(i, typing.Tuple):
            self.multiple  = True
            self.first = i[0]
            self.last = i[1]
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
        return 0

    '''Get Event Info as a dataclass'''
    @abstractmethod
    def eventInfo(self) -> Union[EventInfo,Sequence[EventInfo]]:
        pass

    @abstractmethod
    def wfs(self, calibrated : bool = False) -> numpy.ndarray:  
        pass



# This is not a class, but a factory method! 
def Dataset(station, run, data_dir = None, backend="auto"): 

   if data_dir is None: 
       data_dir = os.environ['RNO_G_DATA'] 

   if backend == "auto":
        try: 
            import ROOT 
            import mattak.backends.pyroot.mattakloader 
            print('Using pyroot backend') 
            backend = "pyroot" 
        except: 
            try: 
                import uproot
                backend = "uproot" 
                print('Using uproot backend') 
            except:
                print("No backends available")
                return None 
                
   if backend == "uproot": 
        import mattak.backends.uproot.dataset 
        return mattak.backends.uproot.dataset.Dataset(station, run, data_dir)
   else: 
        import mattak.backends.pyroot.dataset 
        return mattak.backends.pyroot.dataset.Dataset(station, run, data_dir) 

