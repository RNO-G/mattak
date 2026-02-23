import ROOT
import sys

import mattak.backends.pyroot.mattakloader  # ensure the library is loaded before importing the class
import mattak.backends.pyroot.dataset

import numpy as np
import cppyy


f = ROOT.TFile(sys.argv[1], "READ")

event_tree = f.Get("events")
print(event_tree)
run_summary = f.Get("run")

event_summary = ROOT.mattak.EventSummary()



event_tree.SetBranchAddress("EventSummary", event_summary)
# t.GetEntry(0)  # We typically only store one entry (one run) in the monitoring tree, so we read the first entry.

print("Station number:", run_summary.fstation_number)
print("Number of events:", run_summary.fevent_count)

# f.close()

# df = ROOT.RDataFrame("events", sys.argv[1])
# print(df)

rmss = []

for i in range(event_tree.GetEntries()):
    event_tree.GetEntry(i)

    rms = np.frombuffer(cppyy.ll.cast['float*'](event_summary.rms_per_channel), dtype=np.float32, count=24)
    # rms = np.array(event_summary.rms_per_channel, dtype=np.float32)
    print(rms)
    rmss.append(rms)


rmss = np.array(rmss).flatten()

from matplotlib import pyplot as plt
fig, ax = plt.subplots()
print(rmss)
ax.hist(rmss, bins=20)
ax.set_xlabel("RMS per channel / ADC counts")
plt.show()