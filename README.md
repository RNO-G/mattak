# Mattak 

Mattak is a Greenlandic meal consisting of whale skin and blubber. Blubber is a
store of energy. That's close enough to data storage... 

Mattak is (will eventually be) a multilingual package containing readers and
helpers for RNO-G in C++, Python (via both uproot and PyROOT) and maybe even JS (via
rootjs). Writers are provided in C++ ( and Python via PyROOT)
These are all in the same package so that they can be kept in concert with each other. 


# Data directory hierarchy 

All data is assumed to reside in a directory named ${RNO_G_DATA}, which in the C++ and Python readers is often accessed via an environmental variable. This may be a mounted directory, or could be an http directory (assuming the HTTP server it resides on supports byte ranges). 


The general schema for ``run-associated" data (i.e., data taken synchronously with a run)  is: 

    ${RNO_G_DATA}/station${STATION_NUMBER}/run${RUN_NUMER}/${PACKET_TYPE}.root 

Config files for each run live in: 

    ${RNO_G_DATA}/station${STATION_NUMBER}/run${RUN_NUMBER}/cfg 

And any auxilliary files we might generate will be under: 

    ${RNO_G_DATA}/station${STATION_NUMBER}/run${RUNNUMBER}/aux 

Note that this might result in a lot of files in a directory. 


# Run-associated Data packets

## Header

Run-associated. This is all metadata associated with a single event, not including waveforms. This is kept separated because some operations need only the metadata and also because we can transfer all metadata but only a portion of waveforms 

## Waveforms

These are all waveforms associated with each event.

## DAQStatus

Asynchronous DAQ-related data (e.g. thresholds, scalers, etc.) recorded every once in a while asynchronous to events but still as part of a run. 

# Non-run associated data packets  (TODO) 


Sensor data (not implemented yet) is stored based on the time collected, as is GNSS:

   ${RNO_G_DATA}/station${STATION_NUMBER}/sensor/${YEAR}/${MONTH}/${DAY} 
   ${RNO_G_DATA}/station${STATION_NUMBER}/gnss/${YEAR}/${MONTH}/${DAY} 

Calibration data may eventually be stored under ${RNO_G_DATA}/station_${STATION_NUMBER}/calib but this is subject to change

There may also be auxiliary files globally in ${RNO_G_DATA} and ${RNO_G_DATA}/station_${STATION_NUMBER} (things like mappings, who knows?)


