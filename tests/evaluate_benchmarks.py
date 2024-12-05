import argparse
import json
import numpy as np
import sys

argparser = argparse.ArgumentParser()
argparser.add_argument('--benchmark-file', default='benchmark.json', help='Path to json file with benchmark data')
argparser.add_argument('--threshold', type=float, default=1.5, help='Maximum increase over median previous benchmark; will raise an error if exceeded.')
argparser.add_argument('tag', type=str, help='Tag of benchmark to test against other benchmarks')

args = argparser.parse_args()

with open(args.benchmark_file, 'r') as f:
    benchmarks = json.load(f)

test_benchmark = benchmarks.pop(args.tag)
exit_code = 0
for key, test_value in test_benchmark.items():
    old_values = []
    for run_tag, run in benchmarks.items():
        if key in run:
            old_values.append(run[key])

    if not len(old_values):
        print(f'Skipping key: {key} with no previous benchmark data.')
        continue

    reference = np.median(old_values)
    print(f'{key:20s} : {test_value*1e3:-7.3f} ms / {reference*1e3:-7.3f} ({test_value/reference*100:-4.0f} %)')
    exit_code += test_value/reference > args.threshold

if exit_code:
    print(f"!!! {exit_code} benchmark tests have failed !!!")

sys.exit(exit_code)
