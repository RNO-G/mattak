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
def extract_all_parameters(self):
    """Example processor to extract all parameters."""
    print("Extracting all parameters")
    self.monitoringData.station_number = self.data.station
    self.monitoringData.run_number = self.data.run
    self.monitoringData.num_events = self.data.N()
    print("Number of events:", self.monitoringData.num_events)
    ie = 0
    for ev in self.data.iterate():
        info = ev[0]
        info_dict = info.__dict__
        self.metadata[ie]= info_dict
        # self.metadata[ie]['wfs'] = ev[1]
        ie += 1
def calculate_vrms(self):
    num_events = self.data.N()
    print("Calculating vrms for", num_events, "events")
    vrms_list = [[0 for ch in range(24)] for ie in range(num_events)]
    ie = 0
    for ev in self.data.iterate():
        wfs = ev[1]
        for ch,wf in enumerate(wfs):
            noise = get_split_trace_noise_RMS(wf, segments=4, lowest=2)
            vrms = get_signal_to_noise_ratio(wf,noise, window_size=3)
            vrms_list[ie][ch] = vrms
        ie += 1
    self.metadata['vrms'] = vrms_list


analyzer = monitoring.MonitoringAnalyzer(directory=RNO_G_DATA,output_dir=HERE,backend="pyroot")
# analyzer.add_processor(monitoring.default_processor)
analyzer.add_processor(extract_all_parameters)
analyzer.add_processor(calculate_vrms)
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
