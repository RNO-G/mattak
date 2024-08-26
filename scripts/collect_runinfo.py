import argparse
# import uproot
import os
import glob
import collections
import configparser
import json
import numpy as np


def read_runinfo_txt(file: str) -> dict:
    # we'll abuse configparser to read the runinfo, but we have to add a dummy
    # section to properly abuse it
    config = configparser.ConfigParser()
    with open(file, 'r') as fruninfo:
        config.read_string('[dummy]\n' + fruninfo.read())

    return {k: config["dummy"][k] for k in config["dummy"]}


def read_flower_gain_code(file: str) -> dict:
    runinfo_dir = os.path.dirname(file)
    flower_gain_code_files = glob.glob(f"{runinfo_dir}/flower_gain_codes.*.txt")

    if len(flower_gain_code_files) == 0:
        return {"flower_gain_codes": None}

    assert len(flower_gain_code_files) == 1, f"Found {len(flower_gain_code_files)} flower_gain_code files in {runinfo_dir}."
    flower_gain_code_file = flower_gain_code_files[0]

    with open(flower_gain_code_file, 'r') as f:
        flower_gain_codes_str = f.readlines()

    flower_gain_codes = flower_gain_codes_str[1].strip("\n")
    return {"flower_gain_codes": flower_gain_codes}


def read_comment(file: str) -> dict:
    runinfo_dir = os.path.dirname(file)
    comment_files = glob.glob(f"{runinfo_dir}/comment.txt")

    if len(comment_files) == 0:
        return {"comment": None}

    assert len(comment_files) == 1, f"Found {len(comment_files)} comment files in {runinfo_dir}."
    comment_file = comment_files[0]

    with open(comment_file, 'r') as f:
        comment = f.read()

    return {"comment": comment}



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Collect run information from a list of RNO-G runinfo.txt files. Store the collected information in a single json file.')

    parser.add_argument('run_dir_files', type=str, nargs="+",
                        help='If a single path is provided interpret it as root directory to search for `runinfo.txt` files in any subdirectory. '
                        'If multiple paths are provided interpret them as a list of `runinfo.txt` files.')

    parser.add_argument('--output_file', type=str, default=None, help='The file to write the collected information to. '
                        'If None (default), the name will be determine from the provided input files.')

    args = parser.parse_args()


    filelist = args.run_dir_files
    if len(filelist) == 1 and os.path.isdir(filelist[0]):
        filelist = glob.glob(filelist[0] + "**/*runinfo.txt", recursive=True)

    print(f"Founds {len(filelist)} runinfo files.")

    data = collections.defaultdict(list)

    run_info_keys = [
        'station', 'run', 'run-start-time', 'librno-g-git-hash', 'rno-g-ice-software-git-hash', 'free-space-mb-output-partition',
        'free-space-mb-runfile-partition', 'radiant-fwver', 'radiant-fwdate', 'radiant-bm-fwver', 'radiant-bm-fwdate', 'radiant-samplerate',
        'flower-fwver', 'flower-fwdate', 'run-end-time', 'flower_gain_codes', 'comment', 'total-number-of-events-written']

    for file in filelist:
        run_info = read_runinfo_txt(file)

        for key in run_info_keys:
            value = run_info.get(key, None)
            data[key].append(value)

        flower_gain_codes = read_flower_gain_code(file)
        data["flower_gain_codes"].append(flower_gain_codes["flower_gain_codes"])
        comment = read_comment(file)
        data["comment"].append(comment["comment"])

    for key in run_info:
        if key not in run_info_keys:
            print(f"Found unknown key {key} in runinfo files.... skipping.")

    if args.output_file is None:
        runs = sorted([int(run) for run in data["run"]])
        station = np.unique(data["station"])
        assert len(station) == 1, f"Found multiple stations in the runinfo files: {station}"
        output_file = f"runinfo_summary_st{station[0]}_{runs[0]}-{runs[-1]}_{len(runs)}.json"
    else:
        output_file = args.output_file

    print(f"Writing collected run information to {output_file}.")
    json.dump(data, open(output_file, 'w'), indent=4)