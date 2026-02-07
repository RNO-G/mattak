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

station = 23
run = 999
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

    self.metadata["event_info"]["Vrms"] = [noises] if "Vrms" not in self.metadata["event_info"] else self.metadata["event_info"]["Vrms"] + [noises]
    self.metadata["event_info"]["snr"] = [snrs] if "snr" not in self.metadata["event_info"] else self.metadata["event_info"]["snr"] + [snrs]

    
from NuRadioReco.modules.RNO_G import channelBlockOffsetFitter
block_fitter = channelBlockOffsetFitter.channelBlockOffsets()
def fit_block_offsets(self,event):
    pass

from NuRadioReco.modules.RNO_G import channelGlitchDetector
glitch_detector = channelGlitchDetector.channelGlitchDetector()
glitch_detector.begin()
def detect_channel_glitches(self,event):
    """Glitch detection processor."""
    station = event.get_station()
    glitch_detector.run(event,station)
    has_glitch = channelGlitchDetector.has_glitch(event)
    self.metadata["event_info"]["has_glitch"] = [has_glitch] if "has_glitch" not in self.metadata["event_info"] else self.metadata["event_info"]["has_glitch"] + [has_glitch]


analyzer = monitoring.MonitoringAnalyzer(directory=RNO_G_DATA,output_dir=HERE,backend="pyroot")
analyzer.add_processor(monitoring.default_processor)
# analyzer.add_processor(extract_all_parameters)
analyzer.add_processor(calculate_vrms)
analyzer.add_processor(fit_block_offsets)
analyzer.add_processor(detect_channel_glitches)
analyzer.run(station=[station], run=[run], output_file="test_monitoring.root")  

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
f.Close()
