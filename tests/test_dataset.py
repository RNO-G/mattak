import mattak.Dataset
import time
import numpy
import argparse

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Test mattak.")
    argparser.add_argument('--station', type=int, default=0)
    argparser.add_argument('--run', type=int, default=0)
    argparser.add_argument('--data_dir', type=str, default=None)
    argparser.add_argument('--voltage_calibration', '-vc', type=str, default=None)
    argparser.add_argument('--calibrate', action="store_true")
    argparser.add_argument('--backend', nargs='*', help="Which backend(s) to test", default=['pyroot', 'uproot'])
    argparser.add_argument(
        '--no-skip-incomplete', action='store_true',
        help='If True, test the `skip_incomplete=False` option which tries to read (header)data also for events for which no waveforms are present.' )
    argparser.add_argument('--benchmark', type=str, default=None, help="If provided, store benchmark results in this file.")
    args = argparser.parse_args()

    calibrated = args.calibrate
    skip_incomplete = not args.no_skip_incomplete
    print(f"Test settings: voltage_calibration: {calibrated} / skip_incomplete: {skip_incomplete}")

    for backend in args.backend:
        print(f">----- Testing backend: {backend} -----<")
        d = mattak.Dataset.Dataset(
            args.station, args.run, data_path=args.data_dir,
            backend=backend, verbose=True,
            voltage_calibration=args.voltage_calibration,
            skip_incomplete=skip_incomplete)

        print(d.N())
        print(d.eventInfo())
        if skip_incomplete:
            print(d.wfs(calibrated=calibrated))

        d.setEntries((1, 2))
        print(d.eventInfo())
        if skip_incomplete:
            print(d.wfs(calibrated=calibrated))

        mean = 0
        start = time.time()
        for idx, ev in enumerate(d.iterate(calibrated=calibrated)):
            if skip_incomplete:
                mean += numpy.average(ev[1])

        end = time.time()
        time_per_event = (end - start) / (idx + 1)
        print(mean / d.N())
        print("Total time:", end - start)
        print("Time per event:", time_per_event)
        if args.benchmark is not None:
            import os
            import json
            benchmark_file = os.path.join(os.path.dirname(__file__), 'benchmark.json')
            if not os.path.exists(benchmark_file):
                benchmarks = {}
            else:
                with open(benchmark_file, 'r') as f:
                    benchmarks = json.load(f)

            tag = args.benchmark
            # store benchmark results separately per backend / calibration setting / skip_incomplete
            subtag = f'{backend}' + ['', '/cal'][calibrated] + ['/inc', ''][skip_incomplete]
            if not tag in benchmarks:
                benchmarks[tag] = {}
            benchmarks[tag][subtag] = time_per_event

            with open(benchmark_file, 'w') as f:
                json.dump(benchmarks, f, indent=4)

