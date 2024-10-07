import ROOT
import sys
import platform
import os
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level = logging.DEBUG)
from mattak import __path__ as mattak_path

# this is where pip puts the compiled files
mattak_include_path = os.path.join(mattak_path[0], 'build/include/')
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

try:
    import ROOT
    ds = ROOT.mattak.Dataset()
    logger.debug("Found libmattak to be already loaded.")
    loaded = True
except:
    pass


if not loaded:
    if not silent_load(libmattakName):
        loaded_path = "LD_LIBRARY_PATH"
        logger.debug('Successfully found ' + libmattakName + ' in LD_LIBRARY_PATH')
        loaded = True
    elif not silent_load('build/'+libmattakName):
        logger.debug('Successsfully found ' + libmattakName + ' in build')
        loaded_path = "build"
        loaded = True
    else:
        ROOT.gInterpreter.AddIncludePath(mattak_include_path)
        if not silent_load(os.path.join(mattak_path, libmattakName)):
            logger.debug('Successsfully found ' + libmattakName + ' in ' + mattak_path)
            loaded_path = mattak_path
            loaded = True
        else:
          for path in sys.path:
             if not silent_load(path + '/mattak/backends/pyroot/'+libmattakName):
                 logger.debug('Successsfully found ' + libmattakName + ' in ', path)
                 loaded_path = path
                 loaded = True
                 break

if not loaded:
    raise Exception('Could not load '+ libmattakName)
