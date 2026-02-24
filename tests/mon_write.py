import ROOT

import mattak.Dataset
from NuRadioReco.utilities import fft

from collections import defaultdict
import numpy as np
import sys
import pathlib
import os


def assign_numpy_array_to_cpp_vector(cpp_vector, np_array):
    """ Assigns a 1D numpy array to a ROOT std::vector """
    if np_array.ndim != 1:
        raise ValueError("Only 1D numpy arrays are supported.")

    cpp_vector.clear()
    cpp_vector.reserve(len(np_array))
    for val in np_array:
        cpp_vector.push_back(val)



def get_run_summary(dataset):
    run_summary = ROOT.mattak.RunSummary()
    run_summary.frun_number = dataset.run
    run_summary.fstation_number = dataset.station
    run_summary.fevent_count = dataset.N()
    return run_summary

def write_event_summary(event_summary, event_info, wfs):

    rms = np.std(wfs, axis=1).astype(np.float32)
    assign_numpy_array_to_cpp_vector(event_summary.rms, rms)

    amax = np.max(np.abs(wfs), axis=1).astype(np.float32)
    assign_numpy_array_to_cpp_vector(event_summary.max_abs_amplitude, amax)

    event_summary.event_number = event_info.eventNumber


run_dir = pathlib.Path(sys.argv[1])

if not run_dir.is_dir():
    print(f"Error: {run_dir} is not a valid directory.")
    sys.exit(1)

if not (run_dir / "waveforms.root").is_file():
    print(f"Error: {run_dir} does not contain waveforms.root.")
    sys.exit(1)

if not os.access(run_dir, os.W_OK):
    print(f"Error: {run_dir} is not writable.")
    sys.exit(1)

monitoring_file_path = run_dir / "monitoring.root"

dataset = mattak.Dataset.Dataset(data_path=sys.argv[1], backend='pyroot')

f = ROOT.TFile(str(monitoring_file_path), "RECREATE")
t = ROOT.TTree("events", "Event Summarty Tree")

event_summary = ROOT.mattak.EventSummary()
run_summary = get_run_summary(dataset)

t.Branch("EventSummary", event_summary)

event_counts = defaultdict(int)

avg_spectra = defaultdict(lambda: np.zeros((24, 1025), dtype=np.float32))


rms = []
for ev, wfs in dataset.iterate():

    write_event_summary(event_summary, ev, wfs)
    t.Fill()

    event_counts["total"] += 1
    event_counts[ev.triggerType] += 1
    frequencies = fft.freqs(wfs.shape[1], ev.sampleRate)  # function caches...

    specs = np.abs(fft.time2freq(wfs, ev.sampleRate))**2
    avg_spectra["total"] += specs
    avg_spectra[ev.triggerType] += specs


for trigger_type in avg_spectra:
    avg_spectra[trigger_type] /= event_counts.get(trigger_type, 1)

assign_numpy_array_to_cpp_vector(run_summary.frequencies, frequencies)

def fill_spectra(run_summary_obj, trigger_type):
    for i in range(24):
        vec = ROOT.std.vector("float")()
        assign_numpy_array_to_cpp_vector(vec, avg_spectra[trigger_type][i])
        run_summary_obj.push_back(vec)

fill_spectra(run_summary.avg_spectrum, "total")
fill_spectra(run_summary.avg_spectrum_force, "FORCE")
fill_spectra(run_summary.avg_spectrum_rf0, "RADIANT0")
fill_spectra(run_summary.avg_spectrum_rf1, "RADIANT1")
fill_spectra(run_summary.avg_spectrum_lt, "LT")

run_summary.n_events = event_counts.get("total", 0)
run_summary.n_force_triggers = event_counts.get("FORCE", 0)
run_summary.n_rf0_triggers = event_counts.get("RADIANT0", 0)
run_summary.n_rf1_triggers = event_counts.get("RADIANT1", 0)
run_summary.n_lt_triggers = event_counts.get("LT", 0)

f.WriteObject(run_summary, "run")

f.Write()
f.Close()
