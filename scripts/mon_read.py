"""
Example script to read monitoring information from a ROOT file created by mon_write.py.
It demonstrates how to access event-level and run-level information stored in the monitoring.root file format.
"""
import ROOT
import sys
import re

import mattak.backends.pyroot.mattakloader  # ensure the library is loaded before importing the class
import mattak.backends.pyroot.dataset

import matplotlib.pyplot as plt
import numpy as np


def convert_rdf_to_numpy(rdf):
    """ Returns a dictionary of column names to numpy arrays from a ROOT RDataFrame. """
    # Filter out columns that are not part of the EventSummary struct or are ROOT internal fields
    # We also remove the "EventSummary." prefix and any array size annotations like "[24]"
    # for cleaner keys in the resulting dictionary

    column_names = [
        key for key in rdf.GetColumnNames()
        if key.startswith("EventSummary.") and key.split(".")[-1] not in ["fBits", "TObject", "fUniqueID"]
    ]

    data = rf.AsNumpy(column_names)
    data = {re.sub(r"\[\d+\]", "", k.replace("EventSummary.", "")): v for k, v in data.items()}
    return data


f = ROOT.TFile(sys.argv[1], "READ")

event_tree = f.Get("events")
run_summary = f.Get("RunSummary")

print("Station number:", run_summary.station_number)
print("Number of events:", run_summary.n_events)

# Access option 1: directly from the tree
for entry in event_tree:
    event_summary = entry.EventSummary

    event_id = event_summary.event_number
    rms = np.array(event_summary.rms, dtype=np.float32)


# Access option 2: via ROOT's RDataFrame
rf = ROOT.RDataFrame(event_tree)
data = convert_rdf_to_numpy(rf)

fig, ax = plt.subplots(1, 1)

ax.plot(np.array(run_summary.frequencies), np.array(run_summary.avg_spectrum_force[0]))
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("RMS")
ax.legend()
plt.show()