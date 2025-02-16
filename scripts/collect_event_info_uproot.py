import mattak.backends.uproot.dataset as uproot_dataset
import collect_runinfo
import os
import argparse
import uproot
from collections import defaultdict
import numpy as np
import time


class Dataset(uproot_dataset.Dataset):
    def __init__(self, path, read_daq_status=True):

        self.__read_daq_status = read_daq_status
        self._radiantThrs = None
        self._lowTrigThrs = None
        self.run_info = None
        self.skip_incomplete = True

        # Necessary for the eventInfo method
        self._wfs = None
        self.events_with_waveforms = {}

        self.rundir = path
        if not os.path.isdir(path):
            raise ValueError(f"Pass path ({path}) is not a directory")

        file_names = ["headers.root", "daqstatus.root"]
        for file_name in file_names:
            if not os.path.exists(os.path.join(path, file_name)):
                raise ValueError(f"File {file_name} does not exist in {path}")

        self.hd_file = uproot.open("%s/headers.root" % (self.rundir))
        self.hd_tree, self.hd_branch = uproot_dataset.read_tree(self.hd_file, uproot_dataset.header_tree_names)
        self._hds = self.hd_tree[self.hd_branch]

        if self.__read_daq_status:
            self.ds_file = uproot.open("%s/daqstatus.root" % (self.rundir))
            self.ds_tree, self.ds_branch = uproot_dataset.read_tree(self.ds_file, uproot_dataset.daqstatus_tree_names)
            self._dss = self.ds_tree[self.ds_branch]

        self.setEntries((0, self.N()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Collect event info from RNO-G root files')
    parser.add_argument('paths', type=str, nargs="+", help='Path to the directory containing the uproot files')
    parser.add_argument('--read_daq_status', type=bool, default=True, help='Read daq status file')
    parser.add_argument('--output_file', type=str, default="event_info.npz", help='Output file name')
    args = parser.parse_args()

    keys = ["triggerTime", "triggerType", "lowTrigThrs"]

    data = defaultdict(list)

    t0 = time.time()
    for path in args.paths:
        dataset = Dataset(path, args.read_daq_status)
        event_infos = dataset.eventInfo()
        data["run_number"].append(getattr(event_infos[0], "run"))
        data["number_of_events"].append(dataset.N())
        for key in keys:
            data[key].extend([getattr(event_info, key) for event_info in event_infos])

        flower_gain_codes = collect_runinfo.read_flower_gain_code(
            os.path.join(path, "aux/runinfo.txt"))

        data["flower_gain_codes"].append(flower_gain_codes["flower_gain_codes"])


    data = {key: np.array(data[key]) for key in data}
    print(f"Time taken: {time.time() - t0:.2f} s")
    for key in data:
        print(key, data[key].shape)

    np.savez(args.output_file, **data)