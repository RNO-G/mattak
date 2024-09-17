import argparse
import os
import logging
import ROOT

#Loads mattak
ROOT.gSystem.Load('/home/ruben/Documents/software/mattak/install/libmattak.so')

def get_vref(run_folder:str):
    """
    Function that finds the pedestal used for the given run in the pedestal.root file.
    Parameters
    ----------
    run_folder: string
    
    Returns
    -------
    vref: float
        pedestal used in the run, if no pedestal found defaults to 1.5
    """
    if "pedestal.root" in os.listdir(run_folder):
        logging.debug("Opening pedestal file")
        try:
            pedestalFile = ROOT.TFile.Open(f"{run_folder}/pedestal.root", "READ")
            pedestalTree = pedestalFile.pedestals
            pedestalLeaf = pedestalTree.GetLeaf("vbias")
            pedestalTree.GetEntry(0)
            pedestal = pedestalLeaf.GetValue()
        except:
            logging.warning("Unable to open pedestal file / find pedestal value, using default value of 1.5")
            pedestal = 1.5
    
    else:
        logging.warning(f"No pedestal.root file found in folder {run_folder}, using default value of 1.5")
        pedestal = 1.5
    
    return pedestal



if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = '%(prog)s',
                                    description = "Script to convert bias scans (in root format) to voltage calibration constants")
    parser.add_argument("bias_scan",
                        help = "path to bias scan file")
    parser.add_argument("--destination_folder", required = False,
                        default = None,
                        help = "Optional path to folder where to store the resulting volCalConst file \
                                default path is in the same run folder of the given bias scan")
    parser.add_argument("--vref", required = False,
                        type = float,
                        default = None,
                        help = "Reference of the bias scan, the fit will cross the bias scan at the vref.\
                                by default the script looks for the pedestal in the pedestal.root file of the run.")
    args = parser.parse_args()

    if args.bias_scan.split(".")[-1] != "root":
        logging.debug("Detected non-root format of bias scan, loading librno-g")
        # Loads librno-g
        ROOT.gSystem.Load('/home/ruben/Documents/software/librno-g/build/librno-g.so')

    bias_scan_path = os.path.abspath(args.bias_scan)
    bias_scan_directory =  os.path.dirname(bias_scan_path)

    if args.vref is None:
        if "run" in bias_scan_directory:
            vref = get_vref(bias_scan_directory)
        else:
            logging.warning("Provides bias scan was found not to be in a run folder.\
                            An associated pedestal.root file can hence not be found. \
                            Using default value of vref = 1.5")
            vref = 1.5
    else:
        vref = args.vref

    VC = ROOT.mattak.VoltageCalibration(args.bias_scan, vref)
    if args.destination_folder is not None:
        os.chdir(args.destination_folder)
    else:
        os.chdir(bias_scan_directory)                                   # save in same directory as bias scan
    VC.saveFitCoeffsInFile()
    del VC                                                          # probably not necessary but there have been pyroot memory leak issues in the past