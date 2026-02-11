import mattak
import mattak.Dataset
import time
import numpy as np
import argparse
import os
import ROOT
# from ROOT.mattak import monitoring
import monitoring
import mattak.backends.pyroot.mattakloader
from NuRadioReco.utilities import units
# ROOT.gInterpreter.Declare('#include "mattak/Monitoring.h"')
# ROOT.gSystem.Load("libmattak.so")

from NuRadioReco.utilities.trace_utilities import get_split_trace_noise_RMS,get_signal_to_noise_ratio

RNO_G_DATA = os.environ["RNO_G_DATA"]
HERE = os.path.dirname(os.path.abspath(__file__))
# try:
#     daqstatus = ROOT.mattak.DAQStatus()
#     print("Successfully created DAQStatus object:", daqstatus)
# except Exception as e:
#     print("Error creating DAQStatus object:", e)
# try:
#     monitoring_obj = ROOT.mattak.Monitoring(run,station)
#     print("Successfully created Monitoring object:", monitoring_obj)
# except Exception as e:
#     print("Error creating Monitoring object:", e)

## Run MonitoringAnalyzer to generate monitoring.root
from NuRadioReco.utilities.trace_utilities import get_split_trace_noise_RMS,get_signal_to_noise_ratio 
def calculate_vrms(self,event):
    times,traces = event.get_waveforms() ## n-channels, n-samples
    noises = [get_split_trace_noise_RMS(traces[i], segments=4,lowest=2) for i in range(len(traces))] 
    snrs = [get_signal_to_noise_ratio(traces[i],noises[i],window_size=3) for i in range(len(traces))]

    sampling_rate = (1./(times[0][1]-times[0][0]))
    print("event",event.get_id())
    print(f"sampling rate {sampling_rate:.2f}")

    self.metadata["Vrms"] = [noises] if "Vrms" not in self.metadata else self.metadata["Vrms"] + [noises]
    self.metadata["snr"] = [snrs] if "snr" not in self.metadata else self.metadata["snr"] + [snrs]

    
from NuRadioReco.modules.RNO_G.channelBlockOffsetFitter import channelBlockOffsets, fit_block_offsets, _calculate_block_offsets
block_fitter = channelBlockOffsets()
def detect_block_offset(self,event):
    times,traces = event.get_waveforms()
    #(n_events=1, n_channels, n_cuncks)
    block_offsets = _calculate_block_offsets(np.array(traces),block_size=128)
    self.metadata["block_offsets"] = [block_offsets] if "block_offsets" not in self.metadata else self.metadata["block_offsets"] + [block_offsets]

from NuRadioReco.modules.RNO_G import channelGlitchDetector
glitch_detector = channelGlitchDetector.channelGlitchDetector()
glitch_detector.begin()
def detect_channel_glitches(self,event):
    """Glitch detection processor."""
    station = event.get_station()
    glitch_detector.run(event,station)
    has_glitch = channelGlitchDetector.has_glitch(event)
    self.metadata["has_glitch"] = [has_glitch] if "has_glitch" not in self.metadata else self.metadata["has_glitch"] + [has_glitch]
def summarize_metadata(metadata):
    """Summarize the collected metadata."""
    summary = {}
    for key, values in metadata.items():
        if key == "event_info":
            trigger_types = np.array([values[eid]['triggerType'] for eid in values],dtype=str)
            trigger_time = np.array([values[eid]['triggerTime'] for eid in values],dtype=np.float64)
            readout_time = np.array([values[eid]['readoutTime'] for eid in values],dtype=np.float64)
            lowTrigThrs = np.array([values[eid]['lowTrigThrs'] for eid in values],dtype=np.float64)
            # lowTrigThrs1,lowTrigThrs2,lowTrigThrs3,lowTrigThrs4 = lowTrigThrs.T
            types = {"LT","RADIANT","FORCE"}
            summary[key] = {
                "name": "Event information",
                "num_events": len(values),
                "event_ids": list(values.keys()),
                "trigger_types": trigger_types,
                "trigger_time": trigger_time,
                "readout_time": readout_time,
                "lowTrigThrs":lowTrigThrs
            }
            summary["low_trig_thresholds"] = {
                "name": "LT trigger thresholds by channel (0,3)",
                "mean": np.mean(lowTrigThrs,axis=0),
                "std": np.std(lowTrigThrs,axis=0),
                "min": np.min(lowTrigThrs,axis=0),
                "max": np.max(lowTrigThrs,axis=0)
            }
            summary["trigger_type_distribution"] = {
                "name": "Trigger type distribution",
                "types": types,
                "counts": {t: np.sum([t in tt for tt in trigger_types]) for t in types}
            }
        if key == "block_offsets":
            # For block offsets, summarize the distribution of offsets
            all_offsets = np.array(values)
            summary[key] = {
                "name": "Offsets by channel and block (ADC counts)",
                "mean": np.mean(all_offsets,axis=0),
                "std": np.std(all_offsets,axis=0),
                "min": np.min(all_offsets,axis=0),
                "max": np.max(all_offsets,axis=0)
            }
        if key == "has_glitch":
            all_glitches = np.array(values)
            summary[key] = {
                "name": "Glitch presence by event",
                "glitches": all_glitches,
                "num_glitches": np.sum(all_glitches),
                "total_events": len(all_glitches),
                "glitch_rate": np.sum(all_glitches)/len(all_glitches) if len(all_glitches) > 0 else 0.
            }
        if key == "Vrms":
            all_vrms = np.array(values)
            print("all_vrms shape",all_vrms.shape)
            summary[key] = {
                "name": "Vrms by channel",
                "mean": np.mean(all_vrms,axis=0),
                "std": np.std(all_vrms,axis=0),
                "min": np.min(all_vrms,axis=0),
                "max": np.max(all_vrms,axis=0)
            }
        if key == "snr":
            all_snrs = np.array(values)
            summary[key] = {
                "name": "SNR by channel",
                "mean": np.mean(all_snrs,axis=0),
                "std": np.std(all_snrs,axis=0),
                "min": np.min(all_snrs,axis=0),
                "max": np.max(all_snrs,axis=0)
            }
    return summary
def print_summary(summary):
    """Print the summary of metadata."""
    out_string = "=== Metadata Summary ===\n"
    for key, stats in summary.items():
        if key == "block_offsets":
            continue
        if key == "event_info":
            # print("number of events:", stats["num_events"])
            out_string += f"{stats['name']}:\n"
            out_string += f"  Number of events: {stats['num_events']}\n"
        elif key == "trigger_type_distribution":
            out_string += f"  Trigger type distribution:\n"
            for t in stats["types"]:
                out_string += f"    {t}: {stats['counts'][t]:d} events\n"
        elif key == "low_trig_thresholds":
            out_string += f"{stats['name']}:\n"
            for ch in range(stats['mean'].shape[0]):
                out_string += f"  Channel {ch}: mean={stats['mean'][ch]:.2f}, std={stats['std'][ch]:.2f}, min={stats['min'][ch]:.2f}, max={stats['max'][ch]:.2f}\n"

        elif key == "has_glitch":
            out_string += f"\n{stats['name']}:\n"
            # print(f"  Total events: {stats['total_events']}")
            # print(f"  Events with glitches: {stats['num_glitches']}")
            # print(f"  Glitch rate: {stats['glitch_rate']:.2%}")
            out_string += f"  Total events: {stats['total_events']}\n"
            out_string += f"  Events with glitches: {stats['num_glitches']}\n"
            out_string += f"  Glitch rate: {stats['glitch_rate']:.2%}\n"
        else:
            for k,v in stats.items():
                out_string += f"  {k.capitalize()}: {v}\n"
    print(out_string)
    return out_string
def dump_metadata(self,event):
    """Dump metadata to a Monitoring object"""
    obj = self.monitoringData
    metadata = self.metadata
    # if summary is None:
    #     summary = summarize_metadata(metadata)
    ## set single-value metadata fields
    obj.run_number = run
    obj.station_number = station
    ## set event parameter arrays
    eventParameters = {}
    if "Vrms" in metadata:
        vrms = np.array(metadata["Vrms"],dtype=np.float64)
        for ch in range(vrms.shape[1]):
            eventParameters[f"Vrms_{ch}"] = vrms[:,ch]
    if "snr" in metadata:
        snr = np.array(metadata["snr"],dtype=np.float64)
        for ch in range(snr.shape[1]):
            eventParameters[f"snr_{ch}"] = snr[:,ch]
    if "block_offsets" in metadata:
        block_offsets = np.array(metadata["block_offsets"],dtype=np.int32)
        for ch in range(block_offsets.shape[1]):
            for blk in range(block_offsets.shape[2]):
                eventParameters[f"block_offset_ch{ch}_blk{blk}"] = block_offsets[:,ch,blk]
    if "has_glitch" in metadata:
        has_glitch = np.array(metadata["has_glitch"],dtype=bool)
        eventParameters["has_glitch"] = has_glitch
    if "event_info" in metadata:
        trigger_types = np.array([metadata["event_info"][eid]['triggerType'] for eid in metadata["event_info"]],dtype=str)
        trigger_time = np.array([metadata["event_info"][eid]['triggerTime'] for eid in metadata["event_info"]],dtype=np.float64)
        readout_time = np.array([metadata["event_info"][eid]['readoutTime'] for eid in metadata["event_info"]],dtype=np.float64)
        lowTrigThrs = np.array([metadata["event_info"][eid]['lowTrigThrs'] for eid in metadata["event_info"]],dtype=np.float64)
        #
        # print("assigning trigger types",trigger_types[:10])
        strings = ROOT.vector('string')()
        for t in trigger_types:
            strings.push_back(t)
        obj.trigger_type = strings
        # eventParameters["trigger_time"] = trigger_time
        # eventParameters["readout_time"] = readout_time
        for ch in range(lowTrigThrs.shape[1]):
            eventParameters[f"lowTrigThrs_ch{ch}"] = lowTrigThrs[:,ch]
    print("assigning event parameters to Monitoring object")
    obj.eventParameters = eventParameters
    
station = 23
run = 3400
analyzer = monitoring.MonitoringAnalyzer(directory=RNO_G_DATA,output_dir=HERE,backend="pyroot",debug=True)
analyzer.add_processor(monitoring.default_processor)
analyzer.add_processor(calculate_vrms)
analyzer.add_processor(detect_block_offset)
analyzer.add_processor(detect_channel_glitches)
analyzer.add_processor(dump_metadata)

analyzer.run(station=[station], run=[run], output_file="test_monitoring.root")  
# meta = analyzer.metadata
# summary = summarize_metadata(meta)
# summary_text = print_summary(summary)
# dump_metadata(meta, analyzer.monitoringData,summary=summary)
analyzer.monitoringData.run_number = run
print("Test run number",analyzer.monitoringData.run_number)

analyzer.end()


### Open and Read Monitoring.root file for verification
# f = ROOT.TFile.Open(test_monitoring, "READ")
f = ROOT.TFile.Open(os.path.join(HERE,"test_monitoring.root"), "READ")
print("Reading Monitoring.root file:",f.GetName())
obj = f.Get("Monitoring")
print("Retrieved Monitoring object from file:")

print("station number",obj.station_number)
print("run number",obj.run_number)
obj.run_number = 3411
print("run number",obj.run_number)
print("Event parameters:")  # print first 10 values of lowTrigThrs for channel 0
keys = ['Vrms_0','Vrms_1','Vrms_2','Vrms_3','snr_0','snr_1','snr_2','snr_3']
for key in keys:
    print(key,(np.array(obj.eventParameters[key],dtype=np.float64))[:5])
print("run parameters:",obj.runParameters.size()) 
print(obj.runParameters[0].name)
print(" mean",np.array(obj.runParameters[0].mean)[:5])
print(" std",np.array(obj.runParameters[0].std)[:5])
print(" max",np.array(obj.runParameters[0].max)[:5])
print(" min",np.array(obj.runParameters[0].min)[:5])  
#########
f.Close()
