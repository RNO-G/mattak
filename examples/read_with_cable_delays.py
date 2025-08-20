import mattak.Dataset
import argparse

from NuRadioReco.detector import detector


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="")

    argparser.add_argument('data_dir', type=str, default=None)
    argparser.add_argument('--backend', type=str, help="Which backend to use", default='uproot')

    args = argparser.parse_args()

    d = mattak.Dataset.Dataset(0, 0, data_path=args.data_dir, backend=args.backend, verbose=True)

    d.setEntries((1, 20))

    det = detector.Detector(source="rnog_mongo", always_query_entire_description=False)

    wfs = d.wfs(calibrated=False)
    trace_start_times = d.tracesStartTime(det=det)

    print(trace_start_times)