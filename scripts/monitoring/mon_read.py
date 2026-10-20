"""
Example script to read monitoring information from a ROOT file created by mon_write.py.
It demonstrates how to access event-level and run-level information stored in the monitoring.root file format.

Handles both RADIANT and DiDAQ runs, telling the two apart by the trigger types stored in the file.
"""
import ROOT
import cppyy
import sys
import re


import mattak.backends.pyroot.mattakloader  # ensure the library is loaded before importing the class
import mattak.backends.pyroot.dataset

import matplotlib.pyplot as plt
import numpy as np

# (plot label, RunSummary field) per trigger type. mon_write.py fills only the spectra of the
# digitizer the run was taken with, FORCE being the one trigger type common to both.
SPECTRA_RADIANT = [
    ("Force", "avg_spectrum_force"), ("LT", "avg_spectrum_lt"),
    ("RF0", "avg_spectrum_rf0"), ("RF1", "avg_spectrum_rf1")]
SPECTRA_DIDAQ = [
    ("Force", "avg_spectrum_force"), ("Coinc0", "avg_spectrum_didaq_coinc0"),
    ("Coinc1", "avg_spectrum_didaq_coinc1"), ("Surf up", "avg_spectrum_didaq_surf_up"),
    ("Surf down", "avg_spectrum_didaq_surf_down"), ("Deep phased", "avg_spectrum_didaq_deep_phased")]


def convert_rdf_to_numpy(rdf):
    """ Returns a dictionary of column names to numpy arrays from a ROOT RDataFrame. """
    # Filter out columns that are not part of the EventSummary struct or are ROOT internal fields
    # We also remove the "EventSummary." prefix and any array size annotations like "[24]"
    # for cleaner keys in the resulting dictionary

    column_names = [
        key for key in rdf.GetColumnNames()
        if not key.startswith("EventSummary.") and key.split(".")[-1] not in ["fBits", "TObject", "fUniqueID"]
    ]

    data = rf.AsNumpy(column_names)
    data = {re.sub(r"\[\d+\]", "", k.replace("EventSummary.", "")): v for k, v in data.items()}
    return data


f = ROOT.TFile(sys.argv[1], "READ")

event_tree = f.Get("events")
run_summary = f.Get("RunSummary")

# The DiDAQ spectra only exist (are only filled) for a DiDAQ run, and vice versa
is_didaq = len(run_summary.avg_spectrum_didaq_deep_phased) > 0
is_radiant = len(run_summary.avg_spectrum_lt) > 0

if not is_radiant and not is_didaq:
    raise ValueError("Could not identify digitizer, now phased array triggers found ...")

print("Station number:", run_summary.station_number)
print("Number of events:", run_summary.n_events)
print("Digitizer:", "DIDAQ" if is_didaq else "RADIANT")

block_offsets = []
# Access option 1: directly from the tree
for entry in event_tree:
    event_summary = entry.EventSummary

    event_id = event_summary.event_number
    rms = np.array(event_summary.rms, dtype=np.float32)
    block_offsets.append(np.array(event_summary.block_offset))

# Block offsets are RADIANT-only, mon_write.py leaves them empty for DiDAQ events
if not is_didaq:
    block_offsets = np.array(block_offsets, dtype=np.uint16)
    fig, ax = plt.subplots(1, 1)
    ax.hist(block_offsets[:, 0], label="Channel 0")
    ax.set_xlabel("max. abs. block offset")
    ax.legend()

# Access option 2: via ROOT's RDataFrame - CURRENTLY NOT RECOMMANDED AS INT TYPES ARE NOT PROPERLY CONVERTED TO PYTHON INT
rf = ROOT.RDataFrame(event_tree)
data = convert_rdf_to_numpy(rf)


# Access RunSummary information
fig, ax = plt.subplots(1, 1)

frequencies = np.array(run_summary.frequencies)
for label, field in (SPECTRA_DIDAQ if is_didaq else SPECTRA_RADIANT):
    spectrum = getattr(run_summary, field)
    if len(spectrum):  # empty if this trigger type did not fire during the run
        ax.plot(frequencies, np.array(spectrum[0]), lw=1, label=label)

ax.set_xlabel("Frequency (GHz)")
ax.set_ylabel("Average spectrum")
ax.set_yscale("log")
ax.legend(title="Channel 0")

plt.show()