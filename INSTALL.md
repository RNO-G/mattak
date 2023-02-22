# Build from source

## Requirements

c++:
 - cmake3+
 - ROOT>=6.22

Python:
 - dataclasses
 - uproot
 - awkward==1.*


## Make

```

mkdir build
mkdir install
export RNO_G_INSTALL_DIR=PATH/TO/SOURCE/install
cd build
make
make install

```

## Run tests


```
export PYTHONPATH=$PYTHONPATH:PATH/TO/SOURCE/py
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:PATH/TO/SOURCE/install/lib
export RNO_G_DATA="PATH/TO/DATA"
```

It might be necessary to also add the root libdir to the PYTHONPATH..
```

cd tests
python test_dataset.py 
```

