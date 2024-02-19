import mattak.Dataset
import time
import numpy
import os

stat = 23
r = 325

stations = [stat, stat, stat, 0, 0]
runs = [r, r, r, 0, 0]
data_dirs = [None, os.environ["RNO_G_DATA"], None,
             f'{os.environ["RNO_G_DATA"]}/station{stat}/run{r}',
             f'{os.environ["RNO_G_DATA"]}/station{stat}/run{r}/combined.root']
preferred_files = [None, None, "combined", None, None]

for station, run, data_dir, preferred_file in zip(stations, runs, data_dirs, preferred_files):
    print(f"Load datasets with station = {station}, run = {run}, data_dir = {data_dir}, preferred_file = {preferred_file}")
    for backend in ("uproot", "pyroot"):
        print(backend)
        d = mattak.Dataset.Dataset(station, run, data_path=data_dir, backend=backend, preferred_file=preferred_file)
        print(d.N())
        print(d.eventInfo())
        print(d.wfs())

        d.setEntries((1,2))
        print(d.eventInfo())
        print(d.wfs())

        mean = 0
        start = time.time()
        for ev in d.iterate():
            mean += numpy.average(ev[1])
        end = time.time()

        print(mean / d.N())
        print("time:", end - start)
