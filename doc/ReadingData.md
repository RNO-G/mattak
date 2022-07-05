# Reading data with `mattak` 

This is a quick guide to reading data with `mattak.` If you have data in ROOT format, see one of the first sections. If you want to look at raw data, look at the last section.

## Building `mattak` 

If you plan on using C++ or the `PyROOT` backend, this is necessary for now. In
the future the `setup.py` will automatically do this so you only need to worry
about it if you're developing or working in C++. 

Requirements include `ROOT` and `cmake`, and potentially others (to be
documented) Typically, just typing `make` will work. The `Makefile` just wraps
a `cmake` build, so you can also use `cmake` directly. 

There are special requirements if you need to read/convert raw data (see all the way below) 


## Reading data with ROOT/C++, the easy way 

Once `mattak` is built, you can use the `mattak::Dataset` class to facilitate
reading data. You'll need to load the shared library (`libmattak.so` on Linux,
`libmattak.dylib` on MacOS) before doing anything. The `rootlogon.C` in the
`mattak` root directory will automatically load it on Linux.  

If you have the "standard" directory structure, you can do something like: 


    gSystem->Load("path/to/libmattak.so"); // if not already loaded! 
     
    mattak::Dataset d(21 /*station number*/, 
                     1000 /*run number*/, 
                     nullptr  /*Voltage Calibration ,not used yet */, 
                     data_dir /* root of data directory structure, or pass nullptr to use RNO_G_DATA env variable, in which case you can just pass the first two arguments since the latter two default to nullptr*/ 
    d.setEntry(10); /* load the 10th entry */
    d.header(); /*pointer to mattak::Header */ 
    d.raw(); /* raw waveforms, eventually d.calibrated() will work too */ 
    d.raw()->drawWaveforms(); /* pretty plot waveforms */ 
    d.raw()->radiant_data[0] /* access first channels sample values */

Note that this will automatically prefer the full dataset, if available, over the partial dataset. 

If you don't have the full directory structure, but just files from a run, you can do something like: 

    mattak::Dataset d; 
    d.loadDirectory(dir) ; 


## Reading data with ROOT/C++, the manual way

`mattak::Dataset` is just a wrapper over opening files, loading `TTree`s, and setting `TTree` branch pointers accordingly. You of course can always do something like: 
 
    mattak::Waveforms * wf = 0; 
    TFile f("path/to/waveforms.root"); 
    TTree * t = (TTree*) f.Get("waveforms") // note that the name of the tree may vary, sorry!
    t->SetBranchAddress("waveforms",&wf); //note that the name of the branch may vary, sorry! 
    t->GetEntry(10); 


## Reading data with PyROOT by wrapping `mattak::Dataset`

All of the `ROOT` code can be transcribed to `PyROOT`. For example, abbreviated: 

    import ROOT
    ROOT.gSystem.Load("path/to/libmattak.so") 
    d = ROOT.mattak.Dataset(21,1000,0, "path/to/data/dir")
    d.setEntry(10)
    d.raw().drawWaveforms() 


## Reading data with the "stable" Python API 

Here ``stable" is in quotes because it hasn't been finalized yet. 

The plan is to provide a stable Python API that uses either `PyROOT`or `uproot`
backends. Eventually `setup.py` will work (either compiling `mattak` if `ROOT`
is available, otherwise just compiling some bits that can be wrapped with
`uproot`). 

Since calibration is not yet implemented though, you can already use this
interface just by adding the `py` directory under the `mattak` root directory to your `PYTHONPATH`

Then you can do something like:  

    import mattak.Dataset
    d = mattak.Dataset.Dataset(station=21, run=1000, data_dir=datadir) # if you don't pass data_dir, it will use env var RNO_G_DATA, you can also pass a backend explicitly insitead of using auto
    d.setEntries(10) # load the 10th event
    d.eventInfo() # event info 
    d.wfs() # waveform data
    d.setEntries((20,30)) # bulk interface, load the 20th through 29th events
    d.eventInfo() # list of 10 event infos 
    d.wfs() # 10 events worth of data 

Both "uproot" and "pyroot" backends are supported. If you don't pass a backend
to `Dataset` it will use `PyROOT` if it can load `libmattak.so` (WARNING: code
probably doesn't work on Macs yet, someone who has a Mac should add support)
or `uproot` otherwise. The interface is identical. 

Eventually it will be possible to access calibrated data this way. 

# Converting raw data

If you have raw data you must convert (you probably don't unless you're at
Summit or working on the DAQ), you must build `mattak` with `librno-g` support.
The easiest way to do this is to clone and build (just type `make`) `librno-g`
(https://github.com/rno-g/librno-g) in the parent directory of `mattak`. Then
you can `make configure` in `mattak` (which just run `ccmake ../` in `build`)
and enable `LIBRNO_G_SUPPORT`. It will automatically look in the parent directory for it. 

Then when you build, it will add support for constructing `mattak` objects from raw data and for `rno-g-convert` 

