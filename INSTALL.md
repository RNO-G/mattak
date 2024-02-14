# Notes

This guide is only relevant if you want to manually install mattak, for example
because you want to use the C++ version or contribute to development. If you only want to use mattak from Python, you can just run

    pip install git+https://github.com/RNO-G/mattak.git

For users just interesting in reading data you may want to also read doc/ReadingData.md.

mattak is a C++ library that also has a python API. 

The python API supports either a ROOT or uproot backend for reading data. Depending on which backend you use, the requirements may differ. 

# Build from source

## Requirements

c++:
 - make 
 - cmake3+
 - ROOT>=6.22 (For C++ usage or for Python usage using PyROOT backend) 
 - pybind11 (Eventually, for Python usage using the uproot backend that calls C++ code for calibration) 
 - A compiler supporting at least C++11  (possibly C++14 in the future) 

Python 3.6+:
 - numpy
 - dataclasses
 - uproot  (For Python usage using uproot backend) 
 - awkward==1.\* (for Python usage using uproot backend) 


## C++ Compilation

This builds the C++ library with ROOT support. ROOTless support not yet mature so not documented yet. 

mattak uses cmake as a build system, though there is a Makefile to steer cmake. 
So you can:
```
  make # or make configure, if you want to use the TUI to e.g. support converting raw data 
```
To install somewhere: 
```
  make install # this will try to install to RNO_G_INSTALL_DIR if it exists 
```

Alternatively, you can do this manually with cmake 

``` 
mkdir build # or another directory
cd build
cmake ../ # or ccmake if you want to configure using the TUI,  
make
make install # this will try to install to RNO_G_INSTALL_DIR

```

## Setting paths


```
export PYTHONPATH=$PYTHONPATH:PATH/TO/SOURCE/py # for python bindings, until we get setup.py working 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RNO_G_INSTALL_DIR/lib
export RNO_G_DATA="PATH/TO/DATA"  # this can be via http if desired
```

Depending on how you have ROOT installed, and if you are using PyROOT, it might be necessary to also add the root libdir to the PYTHONPATH..

# Make sure something works: 
```
cd tests
python test_dataset.py 
```

