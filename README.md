# Mattak 

Mattak is a Greenlandic meal consisting of whale skin and blubber. Blubber is a
store of energy. That's close enough to data storage... 

Mattak is (will eventually be) a multilingual package containing readers and
helpers for RNO-G in C++, Python (via both uproot and PyROOT) and maybe even JS (via
rootjs). Writers are provided in C++ ( and Python via PyROOT)
These are all in the same package so that they can be kept in concert with each other. 

# Quickstart
The Python interface for mattak can be installed by running:

    pip install git+https://github.com/RNO-G/mattak.git

In order to use the (faster) PyROOT backend, you need a working installation of [ROOT](https://root.cern/install/). If no ROOT installation can be found, only the uproot backend will be installed.

# Logging

Mattak has two independent logging systems, one per language. They do not share
configuration: suppressing one does not affect the other.

## C++ (ROOT message system)

The C++ classes log through ROOT's message system (`TError.h`), using the
global `Info` / `Warning` / `Error` functions with the originating function as
the location, e.g.

    Error in <mattak::Dataset::loadDir>: Failed to load headers.root in /path/to/run

Severities are used as follows: `Error` for operations that failed (e.g. a run
that could not be loaded), `Warning` for recoverable problems (e.g. missing
optional files like pedestals or runinfo), and `Info` for debug/progress
messages, which are only emitted when the `verbose` option is set (e.g.
`DatasetOptions::verbose` or `Dataset::setVerbose`).

Verbosity is controlled globally through ROOT:

  * C++: `gErrorIgnoreLevel = kError;` (suppresses Info and Warning)
  * PyROOT: `ROOT.gErrorIgnoreLevel = ROOT.kError`
  * `.rootrc`: `Root.ErrorIgnoreLevel: Error`
  * Messages can be redirected/captured with `SetErrorHandler()`.

Exceptions: interactive progress output (e.g. the per-channel dots printed
while `VoltageCalibration` recalculates fits) goes directly to stdout, and the
code paths that are also compiled without ROOT (`applyVoltageCalibration` and
the `mattak_noroot` pybind module) write errors to stderr since `TError.h` is
unavailable there.

## Python (standard `logging` module)

The Python package logs through the named logger `"mattak"` (the backends use
child loggers of it), with `error` / `warning` / `debug` used analogously to
the C++ side. It can be configured with the standard mechanisms, e.g.

    import logging
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("mattak").setLevel(logging.WARNING)

## Connecting the two (PyROOT backend)

With the PyROOT backend you get messages from *both* systems (the C++ library
logs through ROOT, the Python wrapper through `logging`). The convenience
function `mattak.Dataset.set_log_level` controls both with a single call: it
sets the level of the `"mattak"` logger and translates it to
`gErrorIgnoreLevel` (process-global, see above):

    import logging
    import mattak.Dataset
    mattak.Dataset.set_log_level(logging.ERROR)   # only show errors, from both systems

The `verbose` argument of `mattak.Dataset.Dataset` is coupled to the log
level: if not given, it defaults to True when the `"mattak"` logger is
enabled for DEBUG (so `set_log_level(logging.DEBUG)` also enables the C++
debug messages of subsequently created datasets); passing `verbose=True`
explicitly conversely lowers the `"mattak"` logger to DEBUG.

One limitation is inherent to ROOT's design: the C++ side cannot be filtered
to mattak only, so e.g. `set_log_level(logging.ERROR)` also silences ROOT
warnings from everything else in the process.

# Data directory hierarchy 

All data is assumed to reside in a directory named ${RNO_G_DATA}, which in the C++ and Python readers is often accessed via an environmental variable. This may be a mounted directory, or could be an http directory (assuming the HTTP server it resides on supports byte ranges). 


The general schema for ``run-associated" data (i.e., data taken synchronously with a run)  is: 

    ${RNO_G_DATA}/station${STATION_NUMBER}/run${RUN_NUMER}/${PACKET_TYPE}.root 

Config files for each run live in: 

    ${RNO_G_DATA}/station${STATION_NUMBER}/run${RUN_NUMBER}/cfg 

And any auxilliary files we might generate will be under: 

    ${RNO_G_DATA}/station${STATION_NUMBER}/run${RUNNUMBER}/aux 

Note that this might result in a lot of files in a directory. 


# Run-associated Data packets

## Header

Run-associated. This is all metadata associated with a single event, not including waveforms. This is kept separated because some operations need only the metadata and also because we can transfer all metadata but only a portion of waveforms 

## Waveforms

These are all waveforms associated with each event.

## DAQStatus

Asynchronous DAQ-related data (e.g. thresholds, scalers, etc.) recorded every once in a while asynchronous to events but still as part of a run. 

# Non-run associated data packets  (TODO) 

Eventually we might have sensor data, calibration data, GNSS data and who knows what else in the same file tree. That will be defined later... 


