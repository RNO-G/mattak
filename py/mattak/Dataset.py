### Python dataset class, agnostic to backend


import os 

from mattak.backends.uproot import Dataset as uprootDataset 
from mattak.backends.pyroot import Dataset as PyROOTDataset

df Dataset(station, run, data_dir = None, backend="uproot"): 

   if data_dir is None: 
       data_dir = os.environ['RNO_G_DATA'] 

    if backend == "uproot": 
        return uprootDataset(station, run, data_dir)
    else 
        return PyROOTDataset(station, run, data_dir) 




