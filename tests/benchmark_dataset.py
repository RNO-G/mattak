import numpy as np
import uproot
import mattak.Dataset
import time
import os

station=11
run=1101
mattakreader = mattak.Dataset.Dataset(
    station, run, 
    # data_dir ='/mnt/pnfs/inbox/',
    data_dir=None,
    backend='pyroot', 
    read_daq_status = False,
    read_run_info = False
)
t0 = time.time()
mattakreader.setEntries((0, None))
mattakreader.eventInfo()
print(f'mattak.Dataset (pyroot) took {time.time()-t0:.3f} s')

mattakreader = mattak.Dataset.Dataset(
    station, run, 
    # data_dir ='/mnt/pnfs/inbox/',
    data_dir=None,
    backend='uproot', 
    read_daq_status = False, 
    read_run_info = False
)

t0=time.time()
mattakreader.setEntries((0, None))
mattakreader.eventInfo()
print(f'mattak.Dataset (uproot backend) took {time.time()-t0:.3f} s')

t0=time.time()
f = uproot.open("%s/station%d/run%d/combined.root" % (os.environ['RNO_G_DATA'], station,run))
hdr = f['combined']['header']
_ = hdr['event_number'].array(library='np')
print(f'uproot took {time.time()-t0:.3f} s')
