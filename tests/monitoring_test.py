import mattak
import mattak.Dataset
import time
import numpy as np
import argparse
import os
import ROOT
# from ROOT.mattak import monitoring
import monitoring
import mattak.backends.pyroot.mattakloader
from NuRadioReco.utilities import units
# ROOT.gInterpreter.Declare('#include "mattak/Monitoring.h"')
# ROOT.gSystem.Load("libmattak.so")

from NuRadioReco.utilities.trace_utilities import get_split_trace_noise_RMS,get_signal_to_noise_ratio

RNO_G_DATA = os.environ["RNO_G_DATA"]
HERE = os.path.dirname(os.path.abspath(__file__))
try:
    daqstatus = ROOT.mattak.DAQStatus()
    print("Successfully created DAQStatus object:", daqstatus)
except Exception as e:
    print("Error creating DAQStatus object:", e)
try:
    monitoring_obj = ROOT.mattak.Monitoring(run,station)
    print("Successfully created Monitoring object:", monitoring_obj)
except Exception as e:
    print("Error creating Monitoring object:", e)

HERE = os.path.dirname(os.path.abspath(__file__))

## Run MonitoringAnalyzer to generate monitoring.root
def extract_all_parameters(self,event):
    """Example processor to extract all parameters."""
    print("Extracting all parameters")

from NuRadioReco.utilities.trace_utilities import get_split_trace_noise_RMS,get_signal_to_noise_ratio 
def calculate_vrms(self,event):
    traces,times = event.get_waveforms()
    noises = [get_split_trace_noise_RMS(traces[i], segments=4,lowest=2) for i in range(len(traces))] 
    snrs = [get_signal_to_noise_ratio(traces[i],noises[i],window_size=3) for i in range(len(traces))]
    # sampling_rate = event.get_sampling_rate()
    sampling_rate = (1./(times[0][1]-times[0][0]))
    print("event",event.get_id(),"traces shape",np.array(traces).shape)
    print(f"sampling rate {sampling_rate:.2f}")

    self.metadata["Vrms"] = [noises] if "Vrms" not in self.metadata else self.metadata["Vrms"] + [noises]
    self.metadata["snr"] = [snrs] if "snr" not in self.metadata else self.metadata["snr"] + [snrs]

    
from NuRadioReco.modules.RNO_G.channelBlockOffsetFitter import channelBlockOffsets, fit_block_offsets, _calculate_block_offsets
block_fitter = channelBlockOffsets()
def detect_block_offset(self,event):
    traces,times = event.get_waveforms()
    #(n_events, n_channels, n_cuncks)
    block_offsets = _calculate_block_offsets(np.array(traces),block_size=128)
    self.metadata["block_offsets"] = [block_offsets] if "block_offsets" not in self.metadata else self.metadata["block_offsets"] + [block_offsets]

from NuRadioReco.modules.RNO_G import channelGlitchDetector
glitch_detector = channelGlitchDetector.channelGlitchDetector()
glitch_detector.begin()
def detect_channel_glitches(self,event):
    """Glitch detection processor."""
    station = event.get_station()
    glitch_detector.run(event,station)
    has_glitch = channelGlitchDetector.has_glitch(event)
    self.metadata["has_glitch"] = [has_glitch] if "has_glitch" not in self.metadata else self.metadata["has_glitch"] + [has_glitch]

def summarize_metadata(metadata):
    """Summarize the collected metadata."""
    summary = {}
    for key, values in metadata.items():
        if key == "block_offsets":
            # For block offsets, summarize the distribution of offsets
            all_offsets = np.array(values)
            summary[key] = {
                "name": "Offsets by channel and block (ADC counts)",
                "mean": np.mean(all_offsets,axis=0),
                "std": np.std(all_offsets,axis=0),
                "min": np.min(all_offsets,axis=0),
                "max": np.max(all_offsets,axis=0)
            }
        if key == "Vrms":
            all_vrms = np.array(values)
            summary[key] = {
                "name": "Vrms by channel (mV)",
                "mean": np.mean(all_vrms,axis=0),
                "std": np.std(all_vrms,axis=0),
                "min": np.min(all_vrms,axis=0),
                "max": np.max(all_vrms,axis=0)
            }
        if key == "snr":
            all_snrs = np.array(values)
            summary[key] = {
                "name": "SNR by channel",
                "mean": np.mean(all_snrs,axis=0),
                "std": np.std(all_snrs,axis=0),
                "min": np.min(all_snrs,axis=0),
                "max": np.max(all_snrs,axis=0)
            }
        if key == "has_glitch":
            all_glitches = np.array(values)
            summary[key] = {
                "name": "Glitch presence by event",
                "glitches": all_glitches,
                "num_glitches": np.sum(all_glitches),
                "total_events": len(all_glitches),
                "glitch_rate": np.sum(all_glitches)/len(all_glitches) if len(all_glitches) > 0 else 0.
            }
    return summary
def print_summary(summary):
    """Print the summary of metadata."""
    for key, stats in summary.items():
        if key == "block_offsets":
            continue
        print(f"\n{stats['name']}:")
        if key == "has_glitch":
            print(f"  Total events: {stats['total_events']}")
            print(f"  Events with glitches: {stats['num_glitches']}")
            print(f"  Glitch rate: {stats['glitch_rate']:.2%}")
        else:
            print(f"  Mean: {stats['mean']}")
            print(f"  Std: {stats['std']}")
            print(f"  Min: {stats['min']}")
            print(f"  Max: {stats['max']}")
station = 23
run = 3400
analyzer = monitoring.MonitoringAnalyzer(directory=RNO_G_DATA,output_dir=HERE,backend="pyroot",debug=True)
analyzer.add_processor(monitoring.default_processor)
# analyzer.add_processor(extract_all_parameters)
analyzer.add_processor(calculate_vrms)
analyzer.add_processor(detect_block_offset)
analyzer.add_processor(detect_channel_glitches)
analyzer.run(station=[station], run=[run], output_file="test_monitoring.root")  
meta = analyzer.metadata
print("Metadata collected during processing:")
summary = summarize_metadata(meta)
print_summary(summary)
### Open and Read Monitoring.root file for verification
# f = ROOT.TFile.Open(test_monitoring, "READ")
f = ROOT.TFile.Open(os.path.join(HERE,"test_monitoring.root"), "READ")
print("Reading Monitoring.root file:",f.GetName())
obj = f.Get("Monitoring")
print("Retrieved Monitoring object from file:")
print("station number",obj.station_number)
print("run number",obj.run_number)
print("readoutDelay",obj.eventParameters['readoutDelay'])

#########
analyzer.end()
f.Close()
