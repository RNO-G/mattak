import mattak.Dataset
import time
import numpy
import argparse


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Test mattak.")
    argparser.add_argument('--station', type=int, default=21)
    argparser.add_argument('--run', type=int, default=300)
    argparser.add_argument('--data_dir', type=str, default=None)
    argparser.add_argument('--voltage_calibration', '-vc', type=str, default=None)
    argparser.add_argument('--backend', nargs='*', help="Which backend(s) to test", default=['pyroot', 'uproot'])
    args = argparser.parse_args()

    calibrated = args.voltage_calibration is not None

    for backend in args.backend:
        print(f">----- Testing backend: {backend} -----<")
        d = mattak.Dataset.Dataset(
            args.station, args.run, data_path=args.data_dir,
            backend=backend, verbose=True,
            voltage_calibration=args.voltage_calibration)

        print(d.N())
        print(d.eventInfo())
        print(d.wfs(calibrated=calibrated))

        d.setEntries((1, 2))
        print(d.eventInfo())
        print(d.wfs(calibrated=calibrated))

        mean = 0
        start = time.time()
        for idx, ev in enumerate(d.iterate(calibrated=calibrated)):
            mean += numpy.average(ev[1])

        end = time.time()
        print(mean / d.N())
        print("Total time:", end - start)
        print("Time per event:", (end - start) / (idx + 1))