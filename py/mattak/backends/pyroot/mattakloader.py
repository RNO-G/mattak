import ROOT
import sys
import platform
import os
import logging
logger = logging.getLogger(__name__)
from mattak import __path__ as mattak_path

loaded = False
loaded_path = None
currentPlatform = platform.platform()

def silent_load(what):
    current = ROOT.gErrorIgnoreLevel
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    ret = ROOT.gSystem.Load(what)
    ROOT.gErrorIgnoreLevel = current
    return ret

if 'macOS' in currentPlatform:
    libmattakName = 'libmattak.dylib'
else:
    libmattakName = 'libmattak.so'

try:
    if ROOT.mattak is not None:
        logger.debug("Found libmattak to be already loaded.")
        loaded = True
except:
    pass

# Candidates to try, in order of decreasing specificity: (library path, matching
# include dir or None). The include dir is only registered with cling for the
# candidate we actually end up loading.
candidates = []

# Install prefix used by CMake (-DRNO_G_INSTALL_DIR / environment variable)
install_dir = os.environ.get('RNO_G_INSTALL_DIR')
if install_dir:
    candidates.append((os.path.join(install_dir, 'lib', libmattakName),
                       os.path.join(install_dir, 'include')))

# Bare name: resolved via ROOT's dynamic path (LD_LIBRARY_PATH, DYLD_LIBRARY_PATH, ...)
candidates.append((libmattakName, None))

# A build directory relative to the current working directory
candidates.append((os.path.join('build', libmattakName), None))

# This is where pip puts the compiled files
candidates.append((os.path.join(mattak_path[0], 'build/lib', libmattakName),
                   os.path.join(mattak_path[0], 'build/include')))

# Anywhere on sys.path
candidates += [(os.path.join(p, 'mattak/backends/pyroot', libmattakName), None)
               for p in sys.path]

if not loaded:
    for lib, include_dir in candidates:
        # Load() returns 0 if it loaded the library, 1 if it was already loaded
        if silent_load(lib) in (0, 1):
            if include_dir is not None:
                ROOT.gInterpreter.AddIncludePath(include_dir)
            loaded_path = lib
            loaded = True
            logger.debug('Successfully loaded %s from %s', libmattakName, lib)
            break

if not loaded:
    raise Exception('Could not load ' + libmattakName + ', tried: '
                    + ', '.join(lib for lib, _ in candidates))
