import ROOT

import mattak.Dataset
from NuRadioReco.utilities import fft,units

from collections import defaultdict
import numpy as np
import sys
import pathlib
import os

from NuRadioReco.modules.RNO_G.channelBlockOffsetFitter import fit_block_offsets
SAMPLING_BLOCK_SIZE = 64
READOUT_SIZE = 2048
OFFSET_BLOCK_SIZE = int(2*SAMPLING_BLOCK_SIZE)
OFFSET_BLOCK_COUNT = int(READOUT_SIZE/OFFSET_BLOCK_SIZE)
## glitch detection methods lifted from NuRadioMC glitchDetector class to remove reader dependency
## https://github.com/nu-radio/NuRadioMC/blob/develop/NuRadioReco/modules/RNO_G/channelGlitchDetector.py
def diff_sq(eventdata,lab4d_sampling_blocksize=SAMPLING_BLOCK_SIZE):
    """
    Returns sum of squared differences of samples across seams of 128-sample chunks.

    `eventdata`: channel waveform
    """
    block_size = lab4d_sampling_blocksize
    twice_block_size = 2 * block_size

    runsum = 0.0
    for chunk in range(len(eventdata) // twice_block_size - 1):
        runsum += (eventdata[chunk * twice_block_size + block_size - 1] - eventdata[chunk * twice_block_size + block_size]) ** 2
    return np.sum(runsum)
def unscramble(trace,lab4d_readout_size=READOUT_SIZE,lab4d_sampling_blocksize=SAMPLING_BLOCK_SIZE):
    """
    Applies an unscrambling operation to the passed `trace`.
    Note: the first and last sampling block are unusable and hence replaced by zeros in the returned waveform.

    Parameters
    ----------
    `trace`: channel waveform
    """

    readout_size = lab4d_readout_size
    block_size = lab4d_sampling_blocksize
    twice_block_size = 2 * block_size

    new_trace = np.zeros_like(trace)

    for i_section in range(len(trace) // block_size):
        section_start = i_section * block_size
        section_end = i_section * block_size + block_size
        if i_section % 2 == 0:
            new_trace[(section_start + twice_block_size) % readout_size :\
                        (section_end + twice_block_size) % readout_size] = trace[section_start:section_end]
        elif i_section > 1:
            new_trace[(section_start - twice_block_size) % readout_size :\
                        (section_end - twice_block_size) % readout_size] = trace[section_start:section_end]
            new_trace[0:block_size] = 0

    return new_trace

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
    run_summary.frun_number = dataset.run
    run_summary.fstation_number = dataset.station
    run_summary.fevent_count = dataset.N()
    return run_summary

def write_event_summary(event_summary, event_info, wfs):
    event_data = defaultdict(np.array)

    rms = np.std(wfs, axis=1).astype(np.float32)
    assign_numpy_array_to_cpp_vector(event_summary.rms, rms)

    amax = np.max(np.abs(wfs), axis=1).astype(np.float32)
    assign_numpy_array_to_cpp_vector(event_summary.max_abs_amplitude, amax)
    
    glitching_test_statitic = []
    block_offsets = []
    vector_of_floats = ROOT.std.vector('float')
    vector_of_offsets = ROOT.std.vector(vector_of_floats)()
    for wf in wfs:
        
        wf_us = unscramble(wf)
        glitch_ts = (diff_sq(wf_us) - diff_sq(wf_us)) / np.var(wf_us)
        glitching_test_statitic = np.append(glitching_test_statitic,glitch_ts)

        offsets = fit_block_offsets(
            wf, block_size=OFFSET_BLOCK_SIZE, sampling_rate=3.2*units.GHz,
            max_frequency=50*units.MHz, mode='auto', return_trace=False,
            maxiter=5, tol=1e-6)
        block_offsets.append(offsets)
        
        offset_vec = ROOT.std.vector("float")()
        assign_numpy_array_to_cpp_vector(offset_vec, offsets)
        vector_of_offsets.push_back(offset_vec)
        
    assign_numpy_array_to_cpp_vector(event_summary.glitching_test_statitic, glitching_test_statitic.astype(np.float32))
    event_summary.block_offset = vector_of_offsets
    
    event_summary.event_number = event_info.eventNumber

    event_data['rms'] = rms
    event_data['maximum_amplitude'] = amax
    event_data['glitching_test_statitic'] = glitching_test_statitic
    event_data['block_offset'] = np.array(block_offsets)

    return event_data



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

avg_block_offset = np.zeros((24,OFFSET_BLOCK_COUNT),dtype=np.float32)

glitch_counts = np.zeros(24,dtype=np.int32)
glitch_ts_cut = 0.0

rms = []
for ev, wfs in dataset.iterate():

    event_data = write_event_summary(event_summary, ev, wfs)
    t.Fill()

    event_counts["total"] += 1
    event_counts[ev.triggerType] += 1
    frequencies = fft.freqs(wfs.shape[1], ev.sampleRate)  # function caches...

    specs = np.abs(fft.time2freq(wfs, ev.sampleRate))**2
    avg_spectra["total"] += specs
    avg_spectra[ev.triggerType] += specs
    avg_block_offset += event_data['block_offset']
    for ch, ts in enumerate(event_data['glitching_test_statitic']):
        if ts > glitch_ts_cut:
            glitch_counts[ch] += 1



for trigger_type in avg_spectra:
    avg_spectra[trigger_type] /= event_counts.get(trigger_type, 1)

avg_block_offset /= event_counts.get(trigger_type, 1)


assign_numpy_array_to_cpp_vector(run_summary.frequencies, frequencies)

def fill_spectra(run_summary_obj, trigger_type):
    for i in range(24):
        vec = ROOT.std.vector("float")()
        assign_numpy_array_to_cpp_vector(vec, avg_spectra[trigger_type][i])
        run_summary_obj.push_back(vec)
def fill_vec_per_channel(run_summary_obj,avg_block_offset):
    for i in range(24):
        vec = ROOT.std.vector("float")()
        assign_numpy_array_to_cpp_vector(vec, avg_block_offset[i])
        run_summary_obj.push_back(vec)
fill_spectra(run_summary.avg_spectrum, "total")
fill_spectra(run_summary.avg_spectrum_force, "FORCE")
fill_spectra(run_summary.avg_spectrum_rf0, "RADIANT0")
fill_spectra(run_summary.avg_spectrum_rf1, "RADIANT1")
fill_spectra(run_summary.avg_spectrum_lt, "LT")
fill_vec_per_channel(run_summary.avg_block_offset,avg_block_offset)


run_summary.n_events = event_counts.get("total", 0)
run_summary.n_force_triggers = event_counts.get("FORCE", 0)
run_summary.n_rf0_triggers = event_counts.get("RADIANT0", 0)
run_summary.n_rf1_triggers = event_counts.get("RADIANT1", 0)
run_summary.n_lt_triggers = event_counts.get("LT", 0)

run_summary.glitch_counts = glitch_counts

f.WriteObject(run_summary, "run")

f.Write()
f.Close()
