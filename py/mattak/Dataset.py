### Python dataset class, agnostic to backend


import os 
from abc import ABC, abstractmethod
from dataclasses import dataclass
import typing
from typing import Sequence, Union
import numpy 
import datetime 


## Pure python event information. In effect duplicating the ROOT header
@dataclass 
class EventInfo: 
    eventNumber: int 
    station : int
    run: int
    readoutTime : float
    triggerTime : float
    triggerType: str
    sysclk: int
    sysclkLastPPS: typing.Tuple[int,int]  # the last 2 PPS sysclks, most recent first
    pps: int 
    radiantStartWindows: numpy.ndarray

##
''' Abstract Base Class for accessing RNO-G data in Python , implemented either with uproot or PyROOT
  This presents an event-by-event interface. Bulk interface, uses the AbstractBulkDataset''' 

class AbstractDataset(ABC): 
    ''' Select entries to read out with wfs or eventInfo. As a convenience, it's possible to just
        pass an int. 
        Returns the number of entries that are valid 
    ''' 
    @abstractmethod 
    def setEntries( self, i : Union[int,Sequence[int]]):
        pass

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
        from mattak.backends.uproot.dataset import Dataset as uprootDataset 
        return uprootDataset(station, run, data_dir)
   else: 
        from mattak.backends.pyroot.dataset import Dataset as PyROOTDataset
        return PyROOTDataset(station, run, data_dir) 

