import mattak.Dataset
import time

d = mattak.Dataset.Dataset(23, 999, None, verbose=True)
print(d.N())

tigger_filter = lambda event_info: event_info.triggerType == "FORCE"
time_filter = lambda event_info: 1663968535 < event_info.triggerTime < 1663968971

start = time.time()
d.setEntries((0, d.N()))
for ev in d.iterate(selector=tigger_filter):
    ei, wfs = ev
    print(ei.triggerType)

d.setEntries((0, d.N()))
for ev in d.iterate(selector=time_filter):
    ei, wfs = ev
    print(ei.triggerTime)

end = time.time()
print ("time:", end - start)



