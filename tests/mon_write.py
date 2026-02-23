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

    rms = np.std(wfs, axis=1).astype(np.float32)
    for i in range(24):
        event_summary.rms_per_channel[i] = rms[i]

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

    write_event_summary(event_summary, ev, wfs)
    t.Fill()


f.WriteObject(run_summary, "run")

f.Write()
f.Close()
