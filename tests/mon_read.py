import ROOT
import sys

import mattak.backends.pyroot.mattakloader  # ensure the library is loaded before importing the class
import mattak.backends.pyroot.dataset

import numpy as np
import cppyy

def convert_rdf_to_numpy(rdf):
    """ Returns a dictionary of column names to numpy arrays from a ROOT RDataFrame. """
    column_names = [
        key for key in rdf.GetColumnNames()
        if key.startswith("EventSummary.") and key.split(".")[-1] not in ["fBits", "TObject", "fUniqueID"]
    ]

    data = rf.AsNumpy(column_names)
    data = {k.replace("EventSummary.", "").replace("[24]", ""): v for k, v in data.items()}
    return data


f = ROOT.TFile(sys.argv[1], "READ")

event_tree = f.Get("events")
run_summary = f.Get("run")

print("Station number:", run_summary.station_number)
print("Number of events:", run_summary.n_events)

# Access option 1: directly from the tree

for entry in event_tree:
    event_summary = entry.EventSummary

    event_id = event_summary.event_number
    rms = np.array(event_summary.rms_per_channel, dtype=np.float32)

# Access option 2: via ROOT's RDataFrame

rf = ROOT.RDataFrame(event_tree)
data = convert_rdf_to_numpy(rf)
print(data)
