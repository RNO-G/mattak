import ROOT
import mattak.backends.pyroot.mattakloader
import mattak.Dataset
from typing import Sequence, Union, Tuple, Optional, Callable, Generator, TypeVar
import numpy
import os.path
import warnings
import logging

logger = logging.getLogger(__name__)

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

# this is needed for newer versions of cppyy to remain backwards compatible, due to changes in uint8_t handling
cppyy.cppdef(" uint8_t* cast_uint8_t(void * x) { return (uint8_t*) x; }")
cast_uint8_t  = cppyy.gbl.cast_uint8_t


cppyy.cppdef(" int16_t* cast_int16_t(void * x) { return (int16_t*) x; }")
cast_int16_t  = cppyy.gbl.cast_int16_t

cppyy.cppdef(" uint16_t* cast_uint16_t(void * x) { return (uint16_t*) x; }")
cast_uint16_t  = cppyy.gbl.cast_uint16_t

def isNully(p):
    return p is None or ROOT.AddressOf(p) == 0 or cppyy.gbl.is_nully(p)

def _read(obj):
    try:
        return numpy.array(obj)
    except AttributeError:
        return None

# mattak.Dataset.EventInfo threshold field name -> candidate mattak::DAQStatus (C++) attribute
# names, tried in order (first one present on `daq_status` wins). More than one name covers
# fields that were renamed across DAQStatus ClassDef versions (see DAQStatus.h).
_DAQ_STATUS_THRESHOLD_FIELDS = {
    "radiantThrs": ("radiant_thresholds",),
    "lowTrigThrs": ("lt_trigger_thresholds", "lt_coinc_trigger_thresholds"),
    "lowphasedTrigThrs": ("lt_phased_trigger_thresholds",),
    "didaqCoinThrs": ("didaq_coin_thresholds",),
    "didaqPhasedTrigThrs": ("didaq_phased_trigger_thresholds",),
}


def _read_daq_status_thresholds(daq_status) -> dict:
    """
    Read the per-channel/per-beam threshold arrays off a `mattak::DAQStatus` object.

    For each entry in `_DAQ_STATUS_THRESHOLD_FIELDS`, tries each candidate C++ attribute name in
    order and reads the first one that exists on `daq_status` (via `hasattr`/`getattr`, so this
    stays robust across DAQStatus schema versions that renamed a field). A field with no matching
    attribute at all is left as `None`.

    Parameters
    ----------
    daq_status : ROOT.mattak.DAQStatus

    Returns
    -------
    dict
        Maps each `mattak.Dataset.EventInfo` threshold field name (e.g. "radiantThrs") to a
        `numpy.ndarray` (or `None` if not available). Keys match `EventInfo`'s constructor
        keywords exactly, so the result can be passed straight through as `**kwargs`.
    """
    values = {}
    for field_name, candidate_attrs in _DAQ_STATUS_THRESHOLD_FIELDS.items():
        value = None
        for attr in candidate_attrs:
            if hasattr(daq_status, attr):
                value = _read(getattr(daq_status, attr))
                break
        values[field_name] = value
    return values

def _check_digitizer_enum_in_sync():
    """ Sanity check that `mattak.Dataset.Digitizer` (python) and `mattak::Dataset::digitizer`
    (C++, src/mattak/Dataset.h) still declare the exact same members. A mismatch means the
    installed libmattak and the mattak python package come from different versions. """
    cpp_enum = ROOT.mattak.Dataset.digitizer
    cpp_members = {name: int(getattr(cpp_enum, name)) for name in dir(cpp_enum) if name[0].isupper()}
    py_members = {member.name: int(member.value) for member in mattak.Dataset.Digitizer}
    if cpp_members != py_members:
        raise RuntimeError(
            f"mattak.Dataset.Digitizer (python, {py_members}) is out of sync with "
            f"mattak::Dataset::digitizer (C++, {cpp_members}). Update mattak/Dataset.py to match "
            "src/mattak/Dataset.h.")


# Bits within `hdr.trigger_info.didaq_info.type` identifying which DiDAQ RF trigger source fired.
# Mirrors the RNO_G_TRIGGER_RF_DIDAQ_* constants of rno_g_trigger_type_t in librno-g's rno-g.h
# (the DIDAQ bit itself is not needed here since `didaq_trigger` already confirms it's set).
_DIDAQ_COINC0 = 1 << 2
_DIDAQ_COINC1 = 1 << 6
_DIDAQ_DEEP_PHASED = 1 << 3
_DIDAQ_SURF_UP = 1 << 4
_DIDAQ_SURF_DOWN = 1 << 5


def _get_trigger_type(hdr) -> str:
    """
    Classify a header's trigger source into a short descriptive string.

    Reads `hdr.trigger_info` (see `mattak::TriggerInfo` in TriggerInfo.h) and picks, in order,
    whether the event was an RF trigger from RADIANT, the low-threshold (LT/FLOWER) board, or
    DiDAQ, falling back to a forced or PPS trigger.

    Parameters
    ----------
    hdr : ROOT.mattak.Header
        The header of the event to classify.

    Returns
    -------
    str
        "RADIANT0"/"RADIANT1"/"RADIANTX" for a RADIANT RF trigger (X if which one is unknown),
        "LT" for a low-threshold board RF trigger, one of "DIDAQ_COINC0"/"DIDAQ_COINC1"/
        "DIDAQ_SURF_UP"/"DIDAQ_SURF_DOWN"/"DIDAQ_DEEP_PHASED" for a DiDAQ RF trigger (or
        "DIDAQ_X" if `didaq_trigger` is set but none of the known source bits match), "FORCE"
        for a software trigger, "PPS" for a PPS trigger, or "UNKNOWN" if nothing above matches.
    """
    ti = hdr.trigger_info

    # RADIANT/LT and DiDAQ are different, mutually exclusive digitizers -- mattak::Header
    # (Header.cc) should never set both for the same event. Catch it here rather than silently
    # picking one, since it would indicate a bug in the header trigger-type decoding.
    assert not (ti.didaq_trigger and (ti.radiant_trigger or ti.lt_trigger)), (
        f"Header for event {hdr.event_number} decoded as both a DiDAQ trigger and a "
        "RADIANT/LT trigger -- this should be impossible, check mattak::Header (Header.cc)")

    triggerType = "UNKNOWN"
    if ti.radiant_trigger:
        which = ti.which_radiant_trigger
        if which == -1:
            which = "X"
        triggerType = "RADIANT" + str(which)
    elif ti.lt_trigger:
        triggerType = "LT"
    elif ti.didaq_trigger:
        t = ti.didaq_info.type
        if t & _DIDAQ_SURF_UP:
            triggerType = "DIDAQ_SURF_UP"
        elif t & _DIDAQ_SURF_DOWN:
            triggerType = "DIDAQ_SURF_DOWN"
        elif t & _DIDAQ_DEEP_PHASED:
            triggerType = "DIDAQ_DEEP_PHASED"
        elif t & _DIDAQ_COINC0:
            triggerType = "DIDAQ_COINC0"
        elif t & _DIDAQ_COINC1:
            triggerType = "DIDAQ_COINC1"
    elif ti.force_trigger:
        triggerType = "FORCE"
    elif ti.pps_trigger:
        triggerType = "PPS"

    return triggerType


_check_digitizer_enum_in_sync()


class Dataset(mattak.Dataset.AbstractDataset):

    def __init__(self, station : int, run : int, data_path : str,
                 verbose : bool = False, skip_incomplete : bool = True,
                 read_daq_status : bool = True, read_run_info : bool = True,
                 preferred_file : Optional[str] = None,
                 voltage_calibration : Optional[Union[str, bool, TypeVar('ROOT.mattak.VoltageCalibration')]] = None,
                 cache_calibration : Optional[bool] = True):
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
        self.skip_incomplete = skip_incomplete
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

        # Look for voltage calibration if None, returns None if not found.
        # Pass voltage_calibration=False to skip searching/loading entirely.
        if voltage_calibration is None:
            voltage_calibration = mattak.Dataset.find_voltage_calibration_for_dataset(self)

        self.has_calib = False
        if voltage_calibration is not False and (isinstance(voltage_calibration, str) or not isNully(voltage_calibration)):
            # the voltage calibration has to be set as member variable. Otherwise the pointer would get deleted to early.
            self.set_calibration(voltage_calibration, cache_calibration=cache_calibration)

        self.data_path = data_path
        self.full = self.ds.isFullDataset()
        self.digitizer = mattak.Dataset.Digitizer(int(self.ds.getDigitizer()))
        self.setEntries(0)

        logger.debug("We think we found station %d run %d", self.station, self.run)

        self.run_info = None
        if isNully(self.ds.info()):
            self.__read_run_info = False
            warnings.warn("Could not read run info")
        elif self.__read_run_info:
            info = self.ds.info()
            self.run_info = mattak.Dataset.RunInfo(
                station=info.station,
                run=info.run,
                run_start_time=info.run_start_time,
                run_end_time=info.run_end_time,
                sampling_rate=info.radiant_sample_rate,
                run_config=f"{self.rundir}/cfg/acq.cfg",
                acq_start=info.acq_start_time,
                acq_stop=info.run_stop_time,
            )
        else:
            pass


    def set_calibration(self, path_or_object, cache_calibration):
        if isinstance(path_or_object, str):
            self.vc = ROOT.mattak.VoltageCalibration()
            self.vc.readFitCoeffsFromFile(path_or_object, cache_tables=cache_calibration)
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

        daq_thresholds = {field_name: None for field_name in _DAQ_STATUS_THRESHOLD_FIELDS}
        if self.__read_daq_status:
            daq_thresholds = _read_daq_status_thresholds(self.ds.status())

        # now use Dataset's faster sample rate getter
        sampleRate = self.ds.sampleRate() / 1000

        hdr = self.ds.header()
        triggerType = _get_trigger_type(hdr)

        radiantStartWindows = None
        readout_delay = None
        didaqStartOffsets = None
        didaqChannelMask = None
        didaqBeamMask = None

        if self.digitizer == mattak.Dataset.Digitizer.RADIANT:
            # The `numpy.copy(...)`` is strictly necessary. Otherwise group access via `dataset.eventInfo()`
            # results in the same `radiantStartWindows` for each event (only for the last event it is correct)
            radiantStartWindows = numpy.copy(numpy.frombuffer(
                    cast_uint8_t(hdr.trigger_info.radiant_info.start_windows),
                    dtype='uint8', count=self.NUM_CHANNELS * 2).reshape(self.NUM_CHANNELS, 2))


            readout_delay = numpy.copy(numpy.around(numpy.frombuffer(
                cppyy.ll.reinterpret_cast['float*'](self.ds.radiantReadoutDelays()),
                dtype = numpy.float32, count=self.NUM_CHANNELS)))
        elif self.digitizer == mattak.Dataset.Digitizer.DIDAQ:
            didaqStartOffsets = numpy.copy(numpy.frombuffer(
                    cast_uint16_t(hdr.trigger_info.didaq_info.start_offsets),
                    dtype='uint16', count=self.NUM_CHANNELS))

            didaqChannelMask = hdr.trigger_info.didaq_info.channel_mask
            didaqBeamMask = hdr.trigger_info.didaq_info.beam_mask

        return mattak.Dataset.EventInfo(
            eventNumber=hdr.event_number,
            station=self.station,
            run=self.run,
            readoutTime=hdr.readout_time,
            triggerTime=hdr.trigger_time,
            triggerType=triggerType,
            sysclk=hdr.sysclk,
            sysclkLastPPS=(hdr.sysclk_last_pps, hdr.sysclk_last_last_pps),
            pps=hdr.pps_num,
            radiantStartWindows=radiantStartWindows,
            sampleRate=sampleRate,
            hasWaveforms=self.ds.rawAvailable(),
            readoutDelay=readout_delay,
            didaqStartOffsets=didaqStartOffsets,
            didaqChannelMask=didaqChannelMask,
            didaqBeamMask=didaqBeamMask,
            **daq_thresholds,
        )


    def eventInfo(self) -> Union[Optional[mattak.Dataset.EventInfo],Sequence[Optional[mattak.Dataset.EventInfo]]]:
        if self.multiple:
            return [self._eventInfo(idx) for idx in range(self.first, self.last)]

        return self._eventInfo(self.entry)

    def _wfs(self, i : int, calibrated: bool = False):
        self.ds.setEntry(i)
        wf = self.ds.calibrated() if calibrated else self.ds.raw()
        if isNully(wf):
            return None

        if self.digitizer == mattak.Dataset.Digitizer.DIDAQ:
            wfs = numpy.frombuffer(cast_uint8_t(wf.didaq_data), dtype="uint8",
                count=self.NUM_CHANNELS * wf.buffer_length).reshape(
                    self.NUM_CHANNELS, wf.buffer_length)
        else:

            if calibrated:
                wfs = numpy.frombuffer(cppyy.ll.cast['double*'](wf.radiant_data), dtype="float64",
                    count=self.NUM_CHANNELS * wf.buffer_length).reshape(
                        self.NUM_CHANNELS, wf.buffer_length)
            else:
                # FS: I think a np.copy is not necessary here because we do it in wfs()
                wfs = numpy.frombuffer(cast_int16_t(wf.radiant_data), dtype="int16",
                    count=self.NUM_CHANNELS * wf.buffer_length).reshape(
                        self.NUM_CHANNELS, wf.buffer_length)
        return wfs


    def wfs(self, calibrated : bool = False) -> Optional[numpy.ndarray]:
        if calibrated and not self.has_calib:
            raise ValueError("You requested a calibrated waveform but no calibration is available")

        # the simple case first
        if not self.multiple:
            # here a copy is needed to avoid overwriting the waveform in memory
            return numpy.copy(self._wfs(self.entry, calibrated))

        if self.last - self.first < 0:
            return None

        out = None
        for entry in range(self.first, self.last):
            this_wfs = self._wfs(entry, calibrated)
            if this_wfs is not None:
                if out is None:
                    out = numpy.zeros((self.last - self.first, *this_wfs.shape), dtype=this_wfs.dtype)
                out[entry-self.first][:][:] = this_wfs

        return numpy.asarray(out, dtype=float)


    def _iterate(
            self, start : int, stop : int, calibrated : bool , max_in_mem : int,
            selectors: Optional[Union[Callable[[mattak.Dataset.EventInfo], bool], Sequence[Callable[[mattak.Dataset.EventInfo], bool]]]] = None,
            override_skip_incomplete : Optional[bool] = None,
            copy : bool = True) -> Generator[Tuple[Optional[mattak.Dataset.EventInfo], Optional[numpy.ndarray]], None, None]:
        """
        PyROOT implementation of the iterator. See `mattak.Dataset.AbstractDataset.iterate`
        for the meaning of the arguments.

        If `copy` is True (default), each yielded waveform array is a copy (cast to float).
        If False, the yielded array references the memory of the underlying ROOT object,
        which is overwritten when the next event is read: the yielded array is only valid
        until the next iteration step (and keeps the raw dtype, i.e., int16 for uncalibrated
        waveforms).
        """

        skip_incomplete = override_skip_incomplete or self.ds.getOpt().partial_skip_incomplete

        def copy_wfs(wfs):
            # numpy.array (unlike numpy.asarray) guarantees a copy also when
            # the dtype already matches (i.e., for calibrated waveforms)
            if wfs is None or not copy:
                return wfs
            return numpy.array(wfs, dtype=float)

        if selectors is not None:
            if not isinstance(selectors, (list, numpy.ndarray)):
                selectors = [selectors]

            for i in range(start, stop):
                evinfo = self._eventInfo(i)
                wfs = self._wfs(i, calibrated)
                if skip_incomplete and wfs is None:
                    continue
                if evinfo is not None and numpy.all([selector(evinfo) for selector in selectors]):
                    yield evinfo, copy_wfs(wfs)
        else:
            for i in range(start, stop):
                wfs = self._wfs(i, calibrated)
                if skip_incomplete and wfs is None:
                    continue
                yield self._eventInfo(i), copy_wfs(wfs)
