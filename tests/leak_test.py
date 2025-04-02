import gc
import mattak.Dataset
import resource
import time
import sys
import logging
import argparse 
logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument('--backend', choices=['uproot','pyroot'], default='pyroot')
parser.add_argument('--station', type=int, default=21)
parser.add_argument('--start_run', type=int, default=0)
parser.add_argument('--end_run', type=int, default=2000)
args = parser.parse_args()

max_N = 0
n_runs = 0

for run in range(args.start_run,args.end_run):
    try:
        d = mattak.Dataset.Dataset(
            station=args.station, run=run,
            skip_incomplete=False,
            read_daq_status=False, read_run_info=False,
            backend=args.backend
        )
        d.setEntries((0, d.N()))
        einfo = d.eventInfo()
        if d.N() > max_N: 
            max_N= d.N()
        print (args.station,run, d.N())
        print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, max_N)


        n_runs += 1
    except (IOError, ReferenceError, KeyError,ValueError) as e:
        logger.warning('run error', exc_info=e)
        continue

