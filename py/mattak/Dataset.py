### Python dataset class, agnostic to backend
import os
import glob
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Sequence, Union, Tuple, Optional, Generator, Callable, TypeVar
import numpy
import logging
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


class AbstractDataset(ABC):
    """
    Abstract Base Class for accessing RNO-G data in Python, implemented either with uproot or PyROOT.
    To see how to initalize a dataset object see the fuction `Dataset` defined below.
    """

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
                logging.warning("You specified a range which is larger than the amount of events stored in this dataset.")
                self.last = self.N()
        else:
            self.multiple = False
            if i >= self.N() or i < 0:
                raise ValueError(f"The entry you specified \"{i}\" does not exist.")
            self.entry = i
            self.first = i
            self.last = i + 1


    def N(self) -> int:
        """ Return the number of events available in this dataset"""
        return 0


    @abstractmethod
    def _iterate(self, start: int , stop : int , calibrated: bool, max_entries_in_mem: int,
                 selector: Optional[Callable[[EventInfo], bool]]) \
                 -> Generator[Tuple[Optional[EventInfo], Optional[numpy.ndarray]], None, None]:
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
            voltage_calibration : Optional[Union[str, TypeVar('ROOT.mattak.VoltageCalibration')]] = None,
            *, data_dir : Optional[str] = None ) -> Optional[AbstractDataset]:
    """

    This is not a class, but a factory method! This is meant to be the interface
    to rule them all for loading RNO-G data. Due to Cosmin's poor initial API
    design, it has become perhaps more complicated than it should be.

    There are a couple of arguments (`station`, `run`, `data_path`, `preferred_file`)
    which control which/how data are read. The following explains this in the order the code
    interpretes the arguments

        * If `data_path` is a file, then that file will be loaded as a "combined root file".

        * If `data_path` is a directory and `station = 0` and `run = 0` (default),
          `data_path` will be interpreted as a directory containing the ROOT files.

        * If `data_path` is a directory and `station` or `run` are non-zero, a
          dataset corresponding to the run stored in `${data_path}/stationX/runY/`
          using `data_path` as the base is returned.

        * If `data_path` is None, it will be set to either `RNO_G_DATA` or `RNO_G_ROOT_DATA`
          environmental variable if they exist on your system (otherwise an error is raised).
          This (typically) also requires that `station` and `run` are non-zero.

        * `data_path` can also be a URL for loading of files via HTTP
          (e.g. `https://user:password@example.com/rno-g-data`), though
          there may be some subtleties about escaping passwords that may differ
          betweeen different backends.

    If `data_path` is not a file, by default, mattak will always try to load the full run
    data stored in `waveforms.root`, `headers.root`, ... . Only if that fails, it will
    "fall back" and try to load `combined.root`. However, you can specify with `preferred_file`
    which root file it should load (will be interpreted as "combined root file"). For example,
    you can set it to "combined" to load combined.root even if full waveforms are available.
    Or, if you have your own subselection in the same format (e.g. for example if you generated
    a file that is only forced triggers) this provides an arguably convenient
    way to load those.

    Other Parameters
    ----------------

    backend : str (Default: "auto")
        The backend can be chosen explicitly (`"pyroot"` or `"uproot"`) or `"auto"` will try
        to use the best one (`"pyroot"` if available, otherwise reverting to `"uproot"`).

    verbose : bool
        Verbose prints out things mostly useful for debugging.

    read_daq_status : bool
        Self-explanatory. Avoiding reading them may speed things up or work around
        the files not being there for some reason.

    read_run_info : bool
        Self-explanatory. Avoiding reading them may speed things up or work around
        the files not being there for some reason.

    skip_incomplete : bool
        Affects what happens when a dataset is incomplete (telemetered - send by satellite).
        If True, will only index fully telemetered events. If False,
        will be indexed by full dataset, but untelemetered events will have
        waveforms of None type (if requested singly) or be all 0's (if requested
        via the bulk interface, as numpy doesn't support jagged ararys).

    voltage_calibration : str
        Path to a voltage calibration file. If None, check for file in run directory.
        (The pyroot backend actually allows to pass an object of type
        `ROOT.mattak.VoltageCalibration`, but this is not possible to
        implement in uproot.)

    data_dir : deprecated
        Left here for backwards compatibility
    """

    # handle deprecated name data_dir
    if data_dir is not None:
        warnings.warn("data_dir is deprecated, use data_path instead. This may be removed in the future, breaking your code.")
        if data_path is not None:
            raise TypeError("Dataset received both data_path and data_dir!")

        data_path = data_dir

    # treat "" as an alias for None
    if data_path == "":
        data_path = None

    if data_path is None:
        for env_var in ['RNO_G_DATA', 'RNO_G_ROOT_DATA']:
            if env_var in os.environ:
                data_path = os.environ[env_var]
                break

        if data_path is None:
            logging.error(
                "Neither `data_path` nor any relevant environmental variable (e.g. RNO_G_DATA) "
                "is defined and I don't know where else to look :(")
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
        return mattak.backends.uproot.dataset.Dataset(
            station, run, data_path, verbose=verbose, skip_incomplete=skip_incomplete, read_daq_status=read_daq_status,
            read_run_info=read_run_info, preferred_file=preferred_file, voltage_calibration=voltage_calibration)

    elif backend == "pyroot":
        import mattak.backends.pyroot.dataset
        return mattak.backends.pyroot.dataset.Dataset(
            station, run, data_path, verbose=verbose, skip_incomplete=skip_incomplete, read_daq_status=read_daq_status,
            read_run_info=read_run_info, preferred_file=preferred_file, voltage_calibration=voltage_calibration)

    else:
        print("Unknown backend (known backends are \"uproot\" and \"pyroot\")")
        return None

def find_voltage_calibration(rundir, station, time):
    """
    Function to find the calibration file that lays closest to given time.
    Returns None if no file was found
    The order of the search is:
        * run directory
        * under RNO_G_DATA/stationX/calibration
    Parameters
    ----------
    rundir : str
        run directory, found by each backend individually
    station : int
        station number, read from runfile to account for station = 0 case
    time: float
        time of run, read as first time in trigger times
            
    Returns
    -------
    vc_list[closest_idx] : str
        location of calibration file closest ib time to run trigger time
    None
        if no calibration file was found
    """
    # try finding a calibration file in the run directory
    vc_list = glob.glob(f"{rundir}/volCalConst*.root/")

    vc_dir = None
    if not vc_list:
        # look in VC constants directory
        for env_var in ["RNO_G_DATA", "RNO_G_ROOT_DATA"]:
            if env_var in os.environ:
                vc_dir = f"{os.environ[env_var]}/calibration/station{station}"
                vc_list = glob.glob(f"{vc_dir}/volCalConst*.root")
                break
        
        if vc_dir is None:
            logging.error(
                "Could not find a directory for the calibration files." 
                "Was RNO_G_DATA or RNO_G_ROOT_DATA defined as a system env variable?")
            return None

        if not vc_list:
            logging.error("Could not find any calibration files")
            return None
        
    # to marginally save time when there is only one file
    if len(vc_list) == 1:
        return vc_list[0]
    vc_basenames = [os.path.basename(vc) for vc in vc_list]
    # extracting bias scan start time from cal_file name
    vc_start_times = [(i, float(re.split("\W+|_", el)[3])) for i, el in enumerate(vc_basenames)]
    closest_idx = min(vc_start_times, key = lambda pair : numpy.abs(pair[1] - time))[0]
    return vc_list[closest_idx]