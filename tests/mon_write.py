import ROOT

import mattak.Dataset

import numpy as np
import sys
import pathlib
import os

def get_run_summary(dataset):
    run_summary = ROOT.mattak.RunSummary()
    run_summary.frun_number = dataset.run
    run_summary.fstation_number = dataset.station
    run_summary.fevent_count = dataset.N()
    return run_summary

def write_event_summary(event_summary, event_info, wfs):

    rms = np.std(wfs, axis=1).astype(np.float32).tolist()
    for i in range(24):
        event_summary.rms_per_channel[i] = int(rms[i])

    return rms

run_dir = pathlib.Path(sys.argv[1])

if not run_dir.is_dir():
    print(f"Error: {run_dir} is not a valid directory.")
    sys.exit(1)

if not (run_dir / "waveforms.root").is_file():
    print(f"Error: {run_dir} does not contain waveforms.root.")
    sys.exit(1)

if not os.access(run_dir, os.W_OK):
    print(f"Error: {run_dir} is not writable.")
    sys.exit(1)

monitoring_file_path = run_dir / "monitoring.root"

dataset = mattak.Dataset.Dataset(data_path=sys.argv[1], backend='pyroot')

f = ROOT.TFile(str(monitoring_file_path), "RECREATE")
t = ROOT.TTree("events", "Event Summarty Tree")

event_summary = ROOT.mattak.EventSummary()
run_summary = get_run_summary(dataset)

t.Branch("EventSummary", event_summary)

rms = []
for ev, wfs in dataset.iterate():

    # rms.append(write_event_summary(event_summary, ev, wfs))
    rms = np.std(wfs, axis=1).astype(np.float32).tolist()
    for i in range(24):
        event_summary.rms_per_channel[i] = int(rms[i])

    t.Fill()
    # print(np.array(event_summary.rms_per_channel))


f.WriteObject(run_summary, "run")

f.Write()
f.Close()


f = ROOT.TFile(str(monitoring_file_path), "READ")
event_tree = f.Get("events")

event_tree.Print()
for entry in event_tree:
    evt = entry.EventSummary
    print(np.array(evt.rms_per_channel))

# event_summary = ROOT.mattak.EventSummary()
# t.SetBranchAddress("EventSummary", event_summary)

# for i in range(event_tree.GetEntries()):
#     event_tree.GetEntry(i)
#     print(np.array(event_summary.rms_per_channel))

# del run_summary
# del event_summary

# from matplotlib import pyplot as plt
# fig, ax = plt.subplots()
# ax.hist(np.array(rms).flatten(), bins=np.linspace(0., 50))
# ax.set_xlabel("RMS per channel / ADC counts")
# plt.show()