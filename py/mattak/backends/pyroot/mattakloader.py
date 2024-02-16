import ROOT 
import sys
import platform
import os
from mattak import __path__ as mattak_path

# this is where pip puts the compiled files
mattak_path = os.path.join(mattak_path[0], 'build/lib/')

loaded = False
loaded_path = None
libmattakName = None
currentPlatform = platform.platform()

def silent_load(what): 
  current = ROOT.gErrorIgnoreLevel
  ROOT.gErrorIgnoreLevel = ROOT.kFatal
  ret = ROOT.gSystem.Load(what)
  ROOT.gErrorIgnoreLevel = current
  return ret
 
if 'macOS' in currentPlatform:
    # print('macOS detected...')
    libmattakName = 'libmattak.dylib'
else: 
    libmattakName = 'libmattak.so'

if not loaded: 
    if not silent_load(libmattakName):
        loaded_path = "LD_LIBRARY_PATH"
        print('Successfully found ' + libmattakName + ' in LD_LIBRARY_PATH')
        loaded = True
    elif not silent_load('build/'+libmattakName):
        print('Successsfully found ' + libmattakName + ' in build')
        loaded_path = "build"
        loaded = True
    elif not silent_load(os.path.join(mattak_path, libmattakName)):
        print('Successsfully found ' + libmattakName + ' in ' + mattak_path)
        loaded_path = mattak_path
        loaded = True
    else:
        for path in sys.path:
            if not silent_load(path + '/mattak/backends/pyroot/'+libmattakName):
                print('Successsfully found ' + libmattakName + ' in ',path)
                loaded_path = path
                loaded = True
                break 

if not loaded: 
    raise Exception('Could not load '+ libmattakName) 