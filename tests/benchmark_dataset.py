import numpy as np
import uproot
import ROOT
import mattak.Dataset
import time
import os
from os import path

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
print ("nentries: " +str( mattakreader.N()))
t0 = time.time()
mattakreader.setEntries((0, None))
mattakreader.eventInfo()
print(f'mattak.Dataset (pyroot backend) took {time.time()-t0:.3f} s')

t0 = time.time()
d =ROOT.mattak.Dataset(station,run) 
for i in range(d.N()): 
    d.setEntry(i)
    d.header() 
print(f'mattak::Dataset (pyroot) took {time.time()-t0:.3f} s')

t0 = time.time()
ROOT.gInterpreter.ProcessLine("mattak::Dataset d (%d,%d); for (int i = 0; i < d.N(); i++) { d.header(); }" % (station,run))
print(f'mattak::Dataset (cling) took {time.time()-t0:.3f} s')

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
fname = "%s/station%d/run%d/combined.root" % (os.environ['RNO_G_DATA'], station,run); 
f = None
tree_name = "combined"
if os.path.exists(fname):
    f = uproot.open(fname)
else:
    fname = "%s/station%d/run%d/headers.root" % (os.environ['RNO_G_DATA'], station,run); 
    tree_name = "header" 
    f = uproot.open(fname)

hdr = f[tree_name]['header']
_ = hdr['event_number'].array()
_ = hdr['station_number'].array()
_ = hdr['run_number'].array()
_ = hdr['trigger_time'].array()
_ = hdr['trigger_info'].array()
_ = hdr['pps_num'].array()
_ = hdr['sysclk'].array()
_ = hdr['sysclk_last_pps'].array()
_ = hdr['sysclk_last_last_pps'].array()
_ = hdr['trigger_info/trigger_info.radiant_info.start_windows[24][2]'].array()

print(f'uproot took {time.time()-t0:.3f} s')
