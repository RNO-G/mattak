import ROOT
import sys

import mattak.backends.pyroot.mattakloader  # ensure the library is loaded before importing the class
import mattak.backends.pyroot.dataset

import numpy as np
import cppyy


f = ROOT.TFile(sys.argv[1], "READ")

event_tree = f.Get("events")
run_summary = f.Get("run")

print("Station number:", run_summary.fstation_number)
print("Number of events:", run_summary.fevent_count)


rmss = []

for entry in event_tree:
    event_summary = entry.EventSummary

    # rms = np.frombuffer(cppyy.ll.cast['float*'](event_summary.rms_per_channel), dtype=np.float32, count=24)

    rms = np.array(event_summary.rms_per_channel, dtype=np.float32)
    rmss.append(rms)


rmss = np.array(rmss).flatten()

from matplotlib import pyplot as plt
fig, ax = plt.subplots()
print(rmss)
ax.hist(rmss, bins=np.linspace(0., 50))
ax.set_xlabel("RMS per channel / ADC counts")
plt.show()