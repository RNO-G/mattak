"""
monitoring.py

This module provides the MonitoringAnalyzer class for reading and processing
combined.root files to generate monitoring.root outputs. It supports
PyROOT (and uproot backends in the future) for ROOT file handling and allows for the addition
of custom data processing functions.

Usage:
    - Initialize the analyzer with a directory containing combined.root files.
    - Add any desired processing functions.
    - Call the run method to process the files and generate outputs.
"""

import numpy as np
import ROOT
import mattak.backends.pyroot.mattakloader
from pathlib import Path
import re
import mattak.Dataset
from NuRadioReco.modules.io.RNO_G.readRNOGDataMattak import readRNOGData
from NuRadioReco.utilities import units

class MonitoringAnalyzer:
    """Analyzer for reading combined.root files and creating monitoring.root output."""
    
    def __init__(self, directory, 
                 output_dir=None,
                 backend=None,
                 debug=False):
        """
        Initialize the analyzer with input directory and output directory.
        
        Args:
            directory: Path to directory containing combined.root files
            output_dir: Path to output directory for monitoring.root files
            backend: ROOT backend to use ("pyroot" or "uproot")
        """
        self.directory = Path(directory)
        if output_dir is None:
            self.output_dir = self.directory / "monitoring"
        else:
            self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.root_files = []
        self.current_file = None
        ## self.data = {} ## mattak.Dataset.Dataset object
        self.readerRNOG = None
        self.metadata = {}
        self.monitoringData = ROOT.mattak.Monitoring()
        self.processors = []
        self.backend = self._detect_backend(backend)
        self.debug = debug
        if self.backend == "uproot":
            print("uproot backend is under development. expect limited functionality.")
        print(f"Using backend: {self.backend}")
        self.station_run_info = {}
    
    def _detect_backend(self,backend):
        """Detect available ROOT backend (pyroot or uproot)."""
        available_backends = []
        try:
            import ROOT
            available_backends.append("pyroot")
        except ImportError:
            pass
        if backend in available_backends:
            return backend
        try:
            import uproot
            available_backends.append("uproot")
        except ImportError:
            raise RuntimeError("Neither PyROOT nor uproot is available")
        print("Available backends:", available_backends)
        return available_backends[0]
    def _detect_station_run_info(self,station=None,run=None):
        """Extract station and run info from directory structure."""
        try:
            # List subdirectories to find station and run patterns
            subdirs = [d.name for d in self.directory.iterdir() if d.is_dir()]
            print(f"Subdirectories found: {subdirs}")
            
            station_dirs = [d for d in subdirs if d.startswith("station")]
            if station_dirs:
                for station_dir in station_dirs:
                    station_num = float(station_dir.replace("station", ""))
                    if station_num not in self.station_run_info and (station is None or station_num in station):
                        self.station_run_info[station_num] = {}
                    
                    # Check for run subdirectories within station
                    for station_dir in station_dirs:
                        station_path = self.directory / station_dir
                        run_dirs = [d.name for d in station_path.iterdir() if d.is_dir() and d.name.startswith("run")]
                        if run_dirs:
                            for run_dir in run_dirs:
                                run_path = station_path / run_dir
                                run_num = float(run_dir.replace("run", ""))
                                if run_num not in self.station_run_info[station_num] and (run is None or run_num in run):
                                    self.station_run_info[station_num][run_num] = {}
                                    self.station_run_info[station_num][run_num]["path"] = str(run_path)
            else:
                # Try to extract from station{st_num}_run{run_num}.root files
                root_files = list(self.directory.glob("*station*_run*.root"))
                if root_files:
                    for root_file in root_files:
                        match = re.search(r'station(\d+)_run(\d+)', root_file.name)
                        if match:
                            station,run = match.groups()
                            print(f"Found station {station}, run {run} in file {root_file.name}")

        except Exception as e:
            print(f"Could not detect station/run info: {e}")
        
        print("Station/Run info:")
        for st in self.station_run_info:
            print(f"Station: {st}")
            for run in self.station_run_info[st]:
                print(f"  Run {run}: \t{self.station_run_info[st][run]['path']}")
        return self.station_run_info
    def run(self,files=None,station=None,run=None,output_file=None):
        """Loop through all combined.root files in directory and process them."""
        if files:
            ## accept user-specified input files
            if not isinstance(files, list):
                files = [files]
            self.root_files = files
        else:
            ## auto-detect input files from directory structure
            self._detect_station_run_info(station=station,run=run)
            ## collect combined.root files
            for st in self.station_run_info:
                for r in self.station_run_info[st]:
                    file_path = self.station_run_info[st][r]["path"]+"/combined.root"
                    ## check if edited combined.root exists, if not skip
                    if not Path(file_path).exists():
                        print(f"File {file_path} does not exist. Skipping station {st}, run {r}.")
                        continue
                    else: 
                        self.root_files.append(file_path)
        print(f"Found {len(self.root_files)} combined.root files to process.")
        if not self.root_files:
            print("No files to process. Exiting.")
            return
        ## initialize rnog reader
        self.readerRNOG = readRNOGData(run_table_path=None, load_run_table=False)

        for file in self.root_files: 
            self.current_file = file  
            print(f"Reading file with backend {self.backend}")

            self.current_file = self.station_run_info[st][r]["path"]+"/combined.root"
            try:
                # self.data = mattak.Dataset.Dataset(int(st), int(r), data_path=str(self.directory), backend=self.backend, preferred_file="combined")
                self.readerRNOG.begin(self.current_file, 
                                      convert_to_voltage=False,
                                      read_calibrated_data=False,
                                      apply_baseline_correction='none',
                                      run_types=["physics"],
                                      run_time_range=None,
                                      max_trigger_rate=0 * units.Hz,
                                      mattak_kwargs={"backend":"pyroot"},
                                      overwrite_sampling_rate=3.2*units.GHz,
                                      max_in_mem=256,use_fallback_time=True
                                      )
                print(f"Loaded data for station {int(st)}, run {int(r)}")

            except:
                print(f"Could not load data for station {int(st)}, run {int(r)}")
            if self.readerRNOG:
                for ie,event in enumerate(self.readerRNOG.run()):
                    self.process_data(event)
                    if self.debug and ie >= 1:
                        print("Debug mode: processed 10 events, stopping.")
                        break  
                if output_file is None:
                    output_file = self.output_dir / "test_monitoring.root"
                else:
                    output_file = self.output_dir / output_file
                self.write_monitoring_root(str(output_file))
    def end(self):
        self.close()
    
    def add_processor(self, processor):
        """Add a processing function to the pipeline."""
        self.processors.append(processor)
    
    def process_data(self, event):
        """Process data through registered processors sequentially."""
        for processor in self.processors:
            processor(self,event)
            print(f"Applied event processor: {processor.__name__}")
    
    def write_monitoring_root(self, output_file):
        """Write processed data to monitoring.root."""

        try:
            if self.backend == "pyroot":
                return self._write_pyroot(output_file)
            else:
                return self._write_uproot(output_file)
        except Exception as e:
            print(f"Error writing file: {e}")
            return False
    
    def _write_pyroot(self, output_file):
        """Write using PyROOT."""
        import ROOT
        
        # Write TObject (monitoringData) if it exists
        if self.monitoringData:
            output = ROOT.TFile(output_file, "RECREATE")
            output.WriteObject(self.monitoringData, "Monitoring")        
            output.Close()
            print(f"Successfully wrote to {output_file}")
            return True
        else:
            print("No monitoring data to write.")
            return False
    
    def _write_uproot(self, output_file):
        """Write using uproot."""
        import uproot
        with uproot.writing.identify.to_TFile(output_file) as file:
            for key, arrays in self.data.items():
                file[key] = uproot.writing.identify.to_TTree(arrays)
        print(f"Successfully wrote to {output_file}")
        return True
    def close(self):
        """Close any open files or resources."""
        if self.readerRNOG:
            self.readerRNOG.end()
            self.readerRNOG = None
        self.current_file = None
        self.metadata = {}
        self.monitoringData = None
    
def default_processor(self,event):
    """Basic processor example."""
    print("Running default processor")
    keys = ["triggerType", "triggerTime", "readoutTime", "radiantThrs", "lowTrigThrs"]
    event_info = self.readerRNOG.get_events_information(keys=keys)
    self.metadata["event_info"] = event_info
    print("Extracted event information:")
    for key in keys:
        try:
            print(f"{key}: {event_info[key]}")
        except KeyError:
            print(f"{key} not found in event information.")

if __name__ == "__main__":
    import os
    data_path = os.environ["RNO_G_DATA"]
    print("RNO_G_DATA:",data_path)
    analyzer = MonitoringAnalyzer(directory=data_path)
    analyzer.add_processor(default_processor)
    analyzer.run(station=[23], run=[999])
