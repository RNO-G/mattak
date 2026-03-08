"""
This script computes monitoring information RNO-G run data. It extracts event-level and run-level information
and stores them in a dedicated monitoring.root file format.

To run the script, provide the path to the run directory containing a waveforms.root and headers.root as an argument:
    python mon_write.py <run_directory>

To read the monitoring.root file, use the following provided python script:
    python mon_read.py <monitoring_file_path>
"""
import ROOT

import mattak.Dataset
from NuRadioReco.utilities import fft, units

from collections import defaultdict
import numpy as np
import sys
import pathlib
import os

from NuRadioReco.modules.RNO_G.channelBlockOffsetFitter import fit_block_offsets
from NuRadioReco.modules.RNO_G.channelGlitchDetector import diff_sq, unscramble

from NuRadioReco.utilities import logging as nu_logging
nu_logging.set_general_log_level(nu_logging.ERROR)  # suppress warnings from NuRadio

NR_CHANNELS = 24
NR_SAMPESRATES = 2048
OFFSET_BLOCK_SIZE = 128


def calculate_glitch_test_statistic(wf):
    """ Calculate a test statistic for glitch detection based on the second derivative of the waveform """
    wf_us = unscramble(wf)
    return (diff_sq(wf_us) - diff_sq(wf_us)) / np.var(wf_us)


def assign_numpy_array_to_cpp_vector(cpp_vector, np_array):
    """ Assigns a 1D numpy array to a ROOT std::vector """
    if np_array.ndim != 1:
        raise ValueError("Only 1D numpy arrays are supported.")

    cpp_vector.clear()
    cpp_vector.reserve(len(np_array))
    for val in np_array:
        cpp_vector.push_back(val)


def get_run_summary(dataset):
    """ Create and return RunSummary object with run/station number and number of events"""
    run_summary = ROOT.mattak.RunSummary()
    run_summary.run_number = dataset.run
    run_summary.station_number = dataset.station
    return run_summary


def write_event_summary(event_summary, event_info, wfs):
    """ Fill the EventSummary object with information from the header (event_info) and waveforms """
    event_summary.event_number = event_info.eventNumber
    event_summary.block_offset.clear()

    rms = np.std(wfs, axis=1).astype(np.float32)
    assign_numpy_array_to_cpp_vector(event_summary.rms, rms)

    amax = np.max(np.abs(wfs), axis=1).astype(np.float32)
    assign_numpy_array_to_cpp_vector(event_summary.max_abs_amplitude, amax)

    glitching_test_statitic = []
    block_offsets = []

    for wf in wfs:

        glitching_test_statitic.append(calculate_glitch_test_statistic(wf))

        offsets = fit_block_offsets(
            wf, block_size=OFFSET_BLOCK_SIZE, sampling_rate=event_info.sampleRate,
            max_frequency=50*units.MHz, mode='auto', return_trace=False,
            maxiter=5, tol=1e-6)
        block_offsets.append(offsets)

        offset_vec = ROOT.std.vector("float")()
        assign_numpy_array_to_cpp_vector(offset_vec, offsets)
        event_summary.block_offset.push_back(offset_vec)

    glitching_test_statitic = np.asarray(glitching_test_statitic, dtype=np.float32)
    assign_numpy_array_to_cpp_vector(event_summary.glitching_test_statitic, glitching_test_statitic)


if len(sys.argv) != 2:
    print("Usage: python mon_write.py <run_directory>")
    sys.exit(1)

run_dir = pathlib.Path(sys.argv[1])

if not run_dir.is_dir():
    print(f"Error: {run_dir} is not a valid directory.")
    sys.exit(1)

if not (run_dir / "waveforms.root").is_file() and not (run_dir / "combined.root").is_file():
    print(f"Error: {run_dir} does not contain the required root files.")
    sys.exit(1)

if not os.access(run_dir, os.W_OK):
    print(f"Error: {run_dir} is not writable.")
    sys.exit(1)

monitoring_file_path = run_dir / "monitoring.root"

dataset = mattak.Dataset.Dataset(data_path=sys.argv[1], backend='pyroot')

# We allow to update an existing monitoring file. This is needed
# as the rno-g-autoconverter script may be run multiple times on the
# same run directory as new data arrives, and we want to avoid losing
# previously computed monitoring information.
if monitoring_file_path.exists():
    f = ROOT.TFile(str(monitoring_file_path), "UPDATE")
    t = f.Get("events")
    run_summary = f.Get("RunSummary")

    if (run_summary.run_number != dataset.run or
        run_summary.station_number != dataset.station):
        print(f"Error: Existing monitoring file has run number {run_summary.run_number} "
              f"and station number {run_summary.station_number}, which do not match "
              f"the current dataset with run number {dataset.run} and station number "
              f"{dataset.station}.")
        sys.exit(1)

    print(f"Monitoring file for run {run_summary.run_number} already exists "
        f"containing {run_summary.n_events} events, will update with new "
        "events if needed.")

    event_summary = ROOT.mattak.EventSummary()
    t.SetBranchAddress("EventSummary", event_summary)

    event_ids = ROOT.RDataFrame(t).AsNumpy(["event_number"])["event_number"]
    update_file = True
else:
    f = ROOT.TFile(str(monitoring_file_path), "CREATE")
    t = ROOT.TTree("events", "Event Summarty Tree")

    run_summary = get_run_summary(dataset)
    event_summary = ROOT.mattak.EventSummary()
    t.Branch("EventSummary", event_summary)

    update_file = False



event_counts = defaultdict(int)
avg_spectra = defaultdict(lambda:
    np.zeros((NR_CHANNELS, NR_SAMPESRATES // 2 + 1), dtype=np.float32))

update_needed = False
rms = []
for ev, wfs in dataset.iterate():

    # Event already exists in file
    if update_file and ev.eventNumber in event_ids:
        continue

    update_needed = True

    write_event_summary(event_summary, ev, wfs)
    t.Fill()

    event_counts["total"] += 1
    event_counts[ev.triggerType] += 1
    frequencies = fft.freqs(wfs.shape[1], ev.sampleRate)  # function caches...

    specs = np.abs(fft.time2freq(wfs, ev.sampleRate))**2
    avg_spectra["total"] += specs
    avg_spectra[ev.triggerType] += specs


if not update_needed:
    print("No new events to add to monitoring file, exiting.")
    f.Close()
    sys.exit(0)

for trigger_type in avg_spectra:
    avg_spectra[trigger_type] /= event_counts.get(trigger_type, 1)

assign_numpy_array_to_cpp_vector(run_summary.frequencies, frequencies)

def fill_spectra(run_summary_obj, trigger_type, prev_event_number):
    """
    Helper function to write 2d numpy array spectra to the corresponding
    field in the run summary, taking into account existing data if the file
    is being updated.
    """
    for i in range(NR_CHANNELS):
        if update_file and event_counts.get(trigger_type, 0):
            w1 = prev_event_number / (prev_event_number + event_counts[trigger_type])
            w2 = event_counts[trigger_type] / (prev_event_number + event_counts[trigger_type])

            avg_spectra[trigger_type][i] = (
                avg_spectra[trigger_type][i] * w2 + np.array(run_summary_obj[i]) * w1)

        vec = ROOT.std.vector("float")()
        assign_numpy_array_to_cpp_vector(vec, avg_spectra[trigger_type][i])

        if not update_file:
            run_summary_obj.push_back(vec)
        else:
            run_summary_obj[i] = vec

fill_spectra(run_summary.avg_spectrum, "total", run_summary.n_events)
fill_spectra(run_summary.avg_spectrum_force, "FORCE", run_summary.n_forced_triggers)
fill_spectra(run_summary.avg_spectrum_rf0, "RADIANT0", run_summary.n_rf0_triggers)
fill_spectra(run_summary.avg_spectrum_rf1, "RADIANT1", run_summary.n_rf1_triggers)
fill_spectra(run_summary.avg_spectrum_lt, "LT", run_summary.n_lt_triggers)

# Adding event counts to run summary
run_summary.n_events += event_counts.get("total", 0)
run_summary.n_forced_triggers += event_counts.get("FORCE", 0)
run_summary.n_rf0_triggers += event_counts.get("RADIANT0", 0)
run_summary.n_rf1_triggers += event_counts.get("RADIANT1", 0)
run_summary.n_lt_triggers += event_counts.get("LT", 0)

t.Write("", ROOT.TObject.kOverwrite)
run_summary.Write("RunSummary", ROOT.TObject.kOverwrite)

f.Close()
