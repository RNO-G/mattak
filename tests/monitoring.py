

import numpy as np
import ROOT
from array import array
import re
import logging
import argparse
import os
from glob import glob 

import time
import functools
def timer(func):
    """ Timer decorator """
    @functools.wraps(func) ## preserve doc string
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        duration = end - start
        print(f"Finished {func.__name__} in {duration:.4f} seconds")
        # logging.info(f"{func.__name__} took {end - start:.4f} seconds")
        return result
    return wrapper
@timer
def load_data(station,run,path):
    """ Load data from specific station, run in path"""
    return mattak.Dataset.Dataset(int(station), int(run), data_path=path, backend='pyroot')
@timer
def make_monitoring(data,station,run,output_file=None):
    """ Parse through data and make monitoring.root"""
    ################ name and create ################ 
    if output_file is None:
        output_file = "monitoring.root"
    output =  ROOT.TFile(output_file, "RECREATE")
    ################ create tree/branches ################ 
    mon = ROOT.mattak.Monitoring(int(station),int(run))
    event_tree = ROOT.TTree("eventSummary", "Event-level monitoring")
    run_tree = ROOT.TTree("runSummary", "Run-level monitoring")
    run_tree.Branch("Monitoring",mon)

    buffer_vrms = array('f',[0.0]*24)
    event_tree.Branch("Vrms",buffer_vrms,"Vrms[24]/F")

    ################ fill values ################ 
    num_events= 0
    all_vrms = []
    ## Main Event Loop
    for ev in data.iterate():
        waveform = np.array(ev[1]) ## shape[24][2048]
        vrms = np.array([ np.mean(wf**2)**0.5 for wf in waveform])
        # vrms = np.mean(waveform**2,axis=-1)**0.5  ## vrms[24] simple vrms calculation for now
        all_vrms.append(vrms)
        buffer_vrms[:] = array('f',vrms) ## put values into buffer

        event_tree.Fill()
        num_events+=1
    mon.num_events = num_events
    mon.rms_per_channel = array('f',np.mean(np.array(all_vrms)**2,axis=0)**0.5)

    run_tree.Fill()
    ################ write and close ################ 
    run_tree.Write()
    event_tree.Write()
    output.Close()
    print(f"Total: {num_events} events")
@timer
def check_file(file_path):
    """ check file structure"""
    f = ROOT.TFile(file_path,"READ")
    list_keys = []
    for key in f.GetListOfKeys():
        list_keys.append(key.GetName())
        print("key's name", key.GetName())
    expect_key = ['runSummary','eventSummary']
    for ek in expect_key:
        if ek not in list_keys:
            raise ValueError(f"{ek} not in List of Keys")
    runSummary = f['runSummary']
    eventSummary = f['eventSummary']

    all_rms = []
    for ie,event in enumerate(eventSummary):
        vrms = np.array(event.Vrms)
        all_rms.append(vrms)
        #print(f"{ie} Vrms\n:",vrms[:4])
    list_branches = [key.GetName() for key in eventSummary.GetListOfBranches()]
    # runSummary.Scan()
    runSummary.GetEntry(0)
    mon = runSummary.Monitoring
    print(f"Station: {mon.station_number}, Run: {mon.run_number}")
    print(f"Events: {mon.num_events}")
    rms_per_channel = np.array(mon.rms_per_channel,dtype=float)
    print('rms:',rms_per_channel[:4])
    print('avg:',(np.mean(np.array(all_rms)**2,axis=0)**0.5)[:4])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--station",type=int,required=True,help="Station number :int")
    parser.add_argument("--run",type=int,required=False,help="Run number :int (Optional)")
    parser.add_argument("--output_dir",type=str,default=None,help="Output directory :string (optional)")
    args = parser.parse_args()

    HERE = os.path.dirname(__file__)
    data_path = os.environ["RNO_G_DATA"] ## path contains station*/run*/*.root
    print("RNO_G_DATA:",data_path)

    glob_path = os.path.join(data_path,f"station{args.station}",f"run{'*' if (args.run == None) else args.run }/")
    paths = glob(glob_path)

    import mattak.backends.pyroot.mattakloader
    import mattak.Dataset

    for path in paths:
        m = re.search(r"station(\d+)/run(\d+)", path)
        st,r = int(m.groups()[0]),int(m.groups()[1])

        data = load_data(st,r,data_path)

        ## create monitoring.root in the run directory
        # output_file = os.path.join(path,"monitoring.root")

        ## create monitoring.root in the monitoring.py's directory
        output_file = os.path.join(HERE,"monitoring.root")

        make_monitoring(data,st,r,output_file=output_file)
        
    ################################################################
    ## Reading objects from file
    check_file(output_file)