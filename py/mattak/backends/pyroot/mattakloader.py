import ROOT 
import sys

loaded = False

def silent_load(what): 
  current = ROOT.gErrorIgnoreLevel
  ROOT.gErrorIgnoreLevel = ROOT.kFatal
  ret = ROOT.gSystem.Load(what)
  ROOT.gErrorIgnoreLevel = current
  return ret

if not loaded: 
    if not silent_load("libmattak.so"):
        print('successfully found libmattak.so in LD_LIBRARY_PATH')
        loaded = True
    elif not silent_load("build/libmattak.so"):
        print('successsfully found libmattak.so in build')
        loaded = True
    else:
        for path in sys.path:
            if not silent_load(path + '/mattak/backends/pyroot/libmattak.so'):
                print('successsfully found libmattak.so in ',path)
                loaded = True
                break 

if not loaded: 
    raise('Could not load libmattak.so') 


