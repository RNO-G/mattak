import mattak
import mattak.Dataset
import time
import numpy
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


# test_object = ROOT.mattak.Dataset()
# print("Successfully created Dataset object:", test_object)
# try:
#     daqstatus = ROOT.mattak.DAQStatus()
#     print("Successfully created DAQStatus object:", daqstatus)
# except Exception as e:
#     print("Error creating DAQStatus object:", e)
try:
    monitoring = ROOT.mattak.Monitoring(run,station,RNO_G_DATA)
    print("Successfully created Monitoring object:", monitoring)
except Exception as e:
    print("Error creating Monitoring object:", e)
# monitoring = ROOT.mattak.Monitoring()
quit()
HERE = os.path.dirname(os.path.abspath(__file__))
test_monitoring = os.path.join(HERE, "test_monitoring.root")
preferred_file= "header"
combinedFile = os.path.join(RNO_G_DATA, f"station{station:d}/run{run}/{preferred_file}.root")
file_dirname = os.path.dirname(combinedFile)
backend = "pyroot"

# import ROOT
# f = ROOT.TFile.Open(combinedFile, "READ")

### Read combined.root file
print(f"Load datasets with station = {station}, run = {run}, data_dir = {RNO_G_DATA}, preferred_file = {preferred_file}")
d = mattak.Dataset.Dataset(station, run, data_path=RNO_G_DATA, backend=backend, preferred_file=preferred_file)
print("Number of events:",d.N())
print(d.eventInfo())
# print(d.wfs())

# d.setEntries((1,2))
# print(d.eventInfo())
# print(d.wfs())

mean = 0
start = time.time()
i = 0
print("Iterate over events:")
for ev in d.iterate():
    # print("i:",i,ev)
    # print(ev[1]) ## ev[0]: eventInfo, ev[1]: wfs
    mean += numpy.average(ev[1])
    i += 1
end = time.time()

print(mean / d.N())
print("time:{:.3f} seconds".format(end - start))
## Create monitoring obeject
monitoring_obj = ROOT.mattak.Monitoring()
print(monitoring_obj)
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

# for key in f.GetListOfKeys():
#     print(key.GetName())
#     obj = key.ReadObj()
#     print(obj)


#########
# f.Close()
