#/usr/bin/env python3

import mattak.Dataset

# from NuRadioReco.modules.RNO_G.channelGlitchDetector import glitch_test_statistic

import awkward as ak
import numpy as np
import argparse
import uproot
import json
import os


def _write_root_file(f, dic, prev_key=""):
    for key, value in dic.items():
        if isinstance(value, dict):

            _write_root_file(f, value, prev_key=prev_key + key + "/")
        else:

            if isinstance(value, (int, float)):
                f[prev_key + key] = ak.Array([value])
            elif isinstance(value, list):
                f[prev_key + key] = ak.Array(value)
            elif isinstance(value, str):
                f[prev_key + key] = value
            else:
                raise TypeError(f"Unsupported type ({type(value)}) to be stored in root file for {prev_key + key}.")



def write_root_file(monitoring_data, fname="monitoring.root"):
    f = uproot.recreate(fname)
    _write_root_file(f, monitoring_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create monitoring data from data run")

    parser.add_argument(
        "rundir",
        type=str,
        help="Path to the run directory"
    )

    args = parser.parse_args()

    assert os.path.isdir(args.rundir), f"Run directory {args.rundir} does not exist (or is not a directory)"

    required_files = ["headers.root", "waveforms.root"]
    for f in required_files:
        assert os.path.isfile(os.path.join(args.rundir, f)), f"Required file {f} not found in run directory {args.rundir}"

    dataset = mattak.Dataset.Dataset(data_path=args.rundir)

    monitoring_data = {
        "header": {
            "run": int(dataset.run),
            "station": int(dataset.station),
            "duration": dataset.duration(),
        },
        "trigger_rate": {
            "total": dataset.trigger_rate(),
            "LT": dataset.trigger_rate("LT"),
            "RADIANT0": dataset.trigger_rate("RADIANT0"),
            "RADIANT1": dataset.trigger_rate("RADIANT1"),
            "RADIANTX": dataset.trigger_rate("RADIANTX"),
            "FORCE": dataset.trigger_rate("FORCE"),
        }
    }


    dataset.setEntries((0, dataset.N()))


    monitoring_data["events"] = {}
    for ev, wfs in dataset.iterate():

        # glitch_ts = [
        #     glitch_test_statistic(wf, 2048, 64) for wf in wfs
        # ]

        vrms = np.std(wfs, axis=-1)


        monitoring_data["events"][str(ev.eventNumber)] = {
            "trigger_type": ev.triggerType,
            # "glitch_ts": glitch_ts,
            "vrms": vrms.tolist(),
        }

    with open("output.json", "w") as f:
        json.dump(monitoring_data, f, indent=4)


    write_root_file(monitoring_data)