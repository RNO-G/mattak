import sys
import json
import uproot
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from mattak.backends.uproot.dataset import read_tree, daqstatus_tree_names


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Read daqstatus.root files')
    parser.add_argument('paths', type=str, nargs="+", help='Path to the daqstatus.root files')
    parser.add_argument('--output_file', type=str, default="daqstatus.json", help='Output file')

    args = parser.parse_args()

    keys = ['lt_trigger_thresholds[4]', 'readout_time_lt']

    data = defaultdict(list)

    for path in sys.argv[1:]:

        try:
            ds_file = uproot.open(path)
            ds_tree, ds_branch = read_tree(ds_file, daqstatus_tree_names)

            ds = ds_tree[ds_branch]

            for key in keys:
                data[key].extend(np.array(ds[key]).tolist())
        except Exception as e:
            print(f"Error reading {path}: {e}")

    with open(args.output_file, 'w') as f:
        json.dump(data, f)