import mattak
import mattak.Dataset
import time
import numpy as np
import argparse
import os
import ROOT
# from ROOT.mattak import monitoring
import mattak.backends.pyroot.mattakloader
# ROOT.gInterpreter.Declare('#include "mattak/Monitoring.h"')
# ROOT.gSystem.Load("libmattak.so")

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
    monitoring = ROOT.mattak.Monitoring(run,station)
    print("Successfully created Monitoring object:", monitoring)
except Exception as e:
    print("Error creating Monitoring object:", e)

# quit())
HERE = os.path.dirname(os.path.abspath(__file__))

preferred_file= "combined"
combined_file = os.path.join(RNO_G_DATA, f"station{station:d}/run{run}/{preferred_file}.root")
file_dirname = os.path.dirname(combined_file)
backend = "pyroot"

### Read combined.root file
print(f"Load datasets with station = {station}, run = {run}, data_dir = {RNO_G_DATA}, preferred_file = {preferred_file}")
d = mattak.Dataset.Dataset(station, run, data_path=RNO_G_DATA, backend=backend, preferred_file="combined")
print("Number of events:",d.N())
print(d.eventInfo())
# print(d.wfs())

mean = 0
start = time.time()
i = 0
print("Iterate over events:")
for ev in d.iterate():
    # print("i:",i,ev)
    # print(ev[1]) ## ev[0]: eventInfo, ev[1]: wfs
    mean += np.average(ev[1])
    i += 1
end = time.time()

print("mean",mean / d.N())
print("time:{:.3f} seconds".format(end - start))
## Create monitoring obeject
monitoring_obj = ROOT.mattak.Monitoring()
print("Monitoring object:")
print(monitoring_obj)
print("run_number:",monitoring_obj.run_number)
print("station_number:",monitoring_obj.station_number)
monitoring_obj.runParameters["test_float"] = 3.14
my_array = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
monitoring_obj.eventParameters["test_array"] = my_array
print("extra runParameters:",monitoring_obj.runParameters["test_float"])
print("extra eventParameters:",monitoring_obj.eventParameters["test_array"])

### Write Monitoring.root fil
# with ROOT.TFile.Open(test_monitoring, "RECREATE") as m:
#     monitoring_tree = ROOT.TTree("Monitoring")
#     monitoring_obj = ROOT.mattak.Monitoring() ##
#     t.Branch("mon", monitoring_obj)
#     for event in d.iterate():
#         # Fill monitoring_obj with relevant data from event
#         monitoring_obj.run_number = event[0]['run']
#         monitoring_obj.radiant_voltage = d.station
#         monitoring_tree.Fill()

### Open and Read Monitoring.root file for verification
# f = ROOT.TFile.Open(test_monitoring, "READ")
f = ROOT.TFile.Open(os.path.join(RNO_G_DATA, "monitoring/test_monitoring.root"), "READ")
print("Reading Monitoring.root file:",f.GetName())
obj = f.Get("Monitoring")
print("Retrieved Monitoring object from file:")
print("station number",obj.station_number)
print("run number",obj.run_number)

#########
f.Close()
