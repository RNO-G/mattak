import mattak.Dataset
import time
import numpy


for backend in ("pyroot", "uproot"):
    print(backend)
    d = mattak.Dataset.Dataset(21, 300, None, backend)
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
