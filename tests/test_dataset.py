import mattak.Dataset
import time
import numpy
import argparse

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Test mattak.")
    argparser.add_argument('--station', type=int, default=21)
    argparser.add_argument('--run', type=int, default=300)
    argparser.add_argument('--data_dir', type=str, default=None)
    argparser.add_argument('--backend', nargs='*', help="Which backend(s) to test", default=['pyroot', 'uproot'])
    args = argparser.parse_args()

    for backend in args.backend:
        print(f">----- Testing backend: {backend} -----<")
        d = mattak.Dataset.Dataset(
            args.station, args.run, data_dir=args.data_dir,
            backend=backend, verbose=True)
        print(d.N())
        print(d.eventInfo())
        print(d.wfs())

        d.setEntries((1,2))
        print(d.eventInfo())
        print(d.wfs())

        mean = 0
        start = time.time()
        for ev in d.iterate():
        #    print ( ev[0].eventNumber,numpy.average(ev[1]))
            mean += numpy.average(ev[1])
        end = time.time()
        print (mean/d.N())
        print ("time:", end-start)

