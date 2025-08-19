# Installation manual

`mattak` is a C++ library that also has a python API. This manual includes a description of how to "Build mattak from source" as well as how to install it via pip. The former is relevant for all those who like to contribute to / develop the C++ version, while the latter is sufficient when using mattak in python. The python API supports either a ROOT or uproot backend for reading data. Depending on which backend you use, the requirements may differ. For users interesting in reading data you may want to also read doc/ReadingData.md.

## Installing with pip

If you only want to use mattak from Python, you can just run
```
  pip install git+https://github.com/RNO-G/mattak.git
```
Another option is to clone the repository and pip install it with the `-e` (`--editable`) flag:
```
  git clone git@github.com:RNO-G/mattak.git
  cd mattak
  pip install -e .
```
This installs mattak but allows you to make changes to the python API for your installtion.


## Build mattak from source

### Requirements

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


### C++ Compilation
This builds the C++ library with ROOT support. ROOTless support not yet mature so not documented yet.

Clone the git repository, create an install directory, declare install directory via enviornmental variable:
```
  git clone git@github.com:RNO-G/mattak.git
  cd mattak
  mkdir install
  export RNO_G_INSTALL_DIR=$PWD/install
```
mattak uses cmake as a build system, so you can do:
```
  mkdir build
  cd build
  cmake ..
  make
  make install
```
It can be that the pybind11 directory is not automatically found. In that case you can tell cmake where to look for it:
```
  cmake .. -Dpybind11_DIR=PATH/TO/YOUR/PYTHON/LIBS/python3.9/site-packages/pybind11/share/cmake/pybind11/
```
Finding the correct path to the python library can be tricky, in particluar if you are using the virtual enviornment because your enviornment might use pybind11 from the system installed python.

There is also a a Makefile to steer cmake. So instead of the above cmake instructions you can do (in the mattak base directory):
```
  make # or make configure, if you want to use the TUI to e.g. support converting raw data
```
To install somewhere:
```
  make install # this will try to install to RNO_G_INSTALL_DIR if it exists
```

### Install with librno-g support
To install mattak with librno-g support (this is for example necessary when converting RNO-G raw data), you have to set the cmake flag `-DLIBRNO_G_SUPPORT=ON`. (It might be necessary to also specify the location of librno-g using `-DLIBRNO_G_PATH`, by default it looks in ../librno-g - relative to the mattak directory)

### Setting paths

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
