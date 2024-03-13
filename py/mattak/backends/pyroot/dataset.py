import ROOT
import mattak.backends.pyroot.mattakloader
import mattak.Dataset
from typing import Sequence, Union, Tuple, Optional, Callable, Generator, TypeVar
import numpy
import os.path
import warnings

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

cppyy.cppdef(" bool is_nully(void *p) { return !p; }")

def isNully(p):
    return p is None or ROOT.AddressOf(p) == 0 or cppyy.gbl.is_nully(p)


class Dataset(mattak.Dataset.AbstractDataset):

    def __init__(self, station : int, run : int, data_path : str,
                 verbose : bool = False, skip_incomplete : bool = True,
                 read_daq_status : bool = True, read_run_info : bool = True,
                 preferred_file : Optional[str] = None,
                 voltage_calibration : Optional[Union[str, TypeVar('ROOT.mattak.VoltageCalibration')]] = None):
        """
        PyROOT backend for the python interface of the mattak Dataset. See further information in
        `mattak.Dataset.Dataset`.
        """

        self.backend = "pyroot"
        self.__read_daq_status = read_daq_status
        self.__read_run_info = read_run_info

        opt = ROOT.mattak.DatasetOptions()
        self.ds = ROOT.mattak.Dataset()

        opt.partial_skip_incomplete = skip_incomplete
        opt.verbose = verbose
        if preferred_file is not None and preferred_file != "":
            opt.file_preference = preferred_file

        if data_path is not None and os.path.isfile(data_path):
            self.ds.loadCombinedFile(data_path, opt)
            self.rundir = os.path.dirname(data_path)  # Just to keep backends compatibile, not actually needed
        else:
            opt.base_data_dir = data_path
            if station == 0 and run == 0:
                self.rundir = data_path
                self.ds.loadDir(data_path, opt)
            else:
                self.rundir = f"{data_path}/station{station}/run{run}"
                self.ds.loadRun(station, run, opt)

        if self.N() < 0:
            raise IOError("Could not load run [data_path: %s, %d %d]" % (data_path, station, run))

        if self.N() == 0:
            warnings.warn("Run is empty?")
            self.station = station
            self.run = run
        else:
            self.station = self.ds.header().station_number
            self.run = self.ds.header().run_number

        if isinstance(voltage_calibration, str) or not isNully(voltage_calibration):
            # the voltage calibration has to be set as member variable. Otherwise the pointer would get deleted to early.
            self.set_calibration(voltage_calibration)
        else:
            if verbose:
                print("Looking for a calibration file")

            time = self.ds.header().trigger_time
            cal_file = mattak.Dataset.find_voltage_calibration(self.rundir, self.station, time)

            if cal_file is not None:
                if verbose:
                    print(f"Found calibration file {cal_file}")

                self.set_calibration(voltage_calibration)
            else:
                if verbose:
                    print("No calibration file found")

                self.has_calib = False

        self.data_path = data_path
        self.setEntries(0)

        if verbose:
            print("We think we found station %d run %d" % (self.station, self.run))

    def set_calibration(self, path_or_object):
        if isinstance(path_or_object, str):
            self.vc = ROOT.mattak.VoltageCalibration()
            self.vc.readFitCoeffsFromFile(path_or_object)
        else:
            self.vc = path_or_object

        self.ds.setCalibration(self.vc)
        self.has_calib = True

    def N(self) -> int:
        return self.ds.N()

    def _eventInfo(self, i : int) -> Optional[mattak.Dataset.EventInfo]:
        #TODO: handle this in C++ code if it's too slow in Python
        if not self.ds.setEntry(i):
            return None

        hdr = self.ds.header()

        if self.__read_daq_status:
            daq_status = self.ds.status()
            radiantThrs = numpy.array(daq_status.radiant_thresholds)
            lowTrigThrs = numpy.array(daq_status.lt_trigger_thresholds)
        else:
            radiantThrs = None
            lowTrigThrs = None

        if self.__read_run_info:
            runinfo = self.ds.info()
            sampleRate = runinfo.radiant_sample_rate / 1000 if not isNully(runinfo) else 3.2
        else:
            sampleRate = None

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

        return numpy.frombuffer(cppyy.ll.cast['double*' if calibrated else 'int16_t*'](wf.radiant_data),
                                dtype = 'float64' if calibrated else 'int16', count=24 * 2048).reshape(24, 2048)


    def wfs(self, calibrated : bool = False) -> Optional[numpy.ndarray]:
        if calibrated and not self.has_calib:
            raise ValueError("You requested a calibrated waveform but no calibration is available")

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

        out = numpy.asarray(out, dtype=float)

        return out


    def _iterate(self, start : int, stop : int, calibrated : bool , max_in_mem : int,
                 selector: Optional[Callable[[mattak.Dataset.EventInfo], bool]] = None) -> Generator[Tuple[Optional[mattak.Dataset.EventInfo], Optional[numpy.ndarray]],None,None]:

        if selector is not None:
            for i in range(start, stop):
                evinfo = self._eventInfo(i)
                if evinfo is not None and selector(evinfo):
                    yield evinfo, self._wfs(i, calibrated)
        else:
            for i in range(start, stop):
                yield self._eventInfo(i), self._wfs(i, calibrated)
        return
