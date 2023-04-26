
# "Easy" Python install for mattak 

# If we detect ROOT, we build mattak and copy
# libmattak.so to the destination (and then the pyroot backend does some hacks
# to find libmattak.so, first trying to import it in case some version has been
# added to LD_LIBRARY_PATH, otherwise checking the same directory).

# We also build a C++ module using pybind11 that doesn't depend on ROOT for use with the uproot backend. 
# So really you need a C++ compiler and cmake installed, not to mention pybind11... 

import os
import sys
import pathlib
import subprocess
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import logging
logger = logging.getLogger('setup.py')

# raise NotImplementedError("setup.py implementation not done") 

try:
    import ROOT 
except ImportError: 
    ROOT = None


# class MattakExtension(Extension): 
#     def __init__(self, name): 
#         Extension.__init__(*self,name, sources=[])


# class MattakBuild(build_ext): 

#     def run(self): 
#         print('Starting Mattak build') 

#         #check if we have cmake available, and a C++ compiler, if not complain 
#         try: 
#             subprocess.check_output(['c++', '--version'])
#         except OSError: 
#             raise RuntimeEror("can't find a C++ compiler. did you really think this was going to work? ") 

#         try: 
#             subprocess.check_output(['cmake', '--version'])
#         except OSError: 
#             raise RuntimeEror("can't find cmake. all hope is lost.") 

#         try: 
#             subprocess.check_output(['make', '--version'])
#         except OSError: 
#             raise RuntimeEror("can't find make. how do you have cmake but no make? Are you a ninja? Haha, get it?") 


#         base_dir = pathlib.Path(__file__).parent.resolve() 
#         dest_dir = self.get_ext_fullpath(ext.name)
   
#         print(dest_dir) 

#         # this just builds mattak, and installs it into the right place! 
#         if ROOT is not None: 
#             print("---ROOT has been detected, building full mattak---") 

#             ROOT_build_dir = self.build_temp + "/mattak-build" 
#             os.makedirs(ROOT_build_dir, exist_ok = True) 
#             subprocess.check_call(['cmake', base_dir ], cwd=ROOT_build_dir)
#             subprocess.check_call(['cmake','--build','.'], cwd=ROOT_build_dir)

#             # now copy libmattak.so to the destination, crazy right? 
#             os.makedirs(dest_dir + "/backends/pyroot", exist_ok = True)
#             self.copy_file(ROOT_build_dir+"/libmattak.so", dest_dir  + "/backends/pyroot")

#         else: 
#             print("---ROOT has NOT been detected. Skipping full mattak---") 
#             print("   (note you may need to uninstall and reinstall if you later install ROOT and want to use the ROOT backend") 


#         #now let's compile the module just used for pybind... 
#         #this lives in py/backends/uproot/_cxx 
#         print("---Building ROOTless Python bindings---") 

#         pybind_build_dir = self.build_temp + "/pybind-build"
#         os.makedirs(pybind_build_dir, exist_ok = True) 
#         subprocess.check_call(['cmake -D ROOTLESS=yes', base_dir ], cwd=pybind_build_dir)
#         subprocess.check_call(['cmake','--build','.'], cwd=pybind_build_dir)
#         os.makedirs(dest_dir + "/backends/uproot", exist_ok = True)
#         print("---Finding a home for them---") 
#         self.copy_file(pybind_build_dir+"/mattak-noroot.so", dest_dir  + "/backends/pyroot/_cxx")

logger.warning(
    '\n'.join([
        '\n',78 * '-',
        'Currently, only the mattak uproot backend can be installed using setup.py!',
        'For better performance, we recommend installing ROOT and installing mattak manually.',
        'See https://github.com/RNO-G/mattak/blob/main/INSTALL.md for instructions.',
        78*'-'
        ])
)
setup (name='mattak', 
       version='0.1',
       url='https://github.com/RNO-G/mattak',
       packages=find_packages('./py'),
       package_dir={'':'py'}
) 


