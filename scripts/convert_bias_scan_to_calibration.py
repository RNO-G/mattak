#! /usr/bin/env python3
import argparse
import os
import logging
import ROOT

#Loads mattak, make sure env variable LD_LIBRARY_PATH includes path to your mattak install
ROOT.gSystem.Load('libmattak.so')
ROOT.gSystem.Load('librno-g.so')

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
    pedestal = 1.5

    for pedestalname in ["pedestal.root", "pedestals.root"]:
        try:
            pedestalFile = ROOT.TFile.Open(f"{run_folder}/{pedestalname}", "READ")
            pedestalTree = pedestalFile.pedestals
            pedestalLeaf = pedestalTree.GetLeaf("vbias")
            pedestalTree.GetEntry(0)
            pedestal = pedestalLeaf.GetValue()
        except:
            continue
        else:
            break

    if pedestal == 1.5:
        logging.warning("Unable to open pedestal(s) file / find pedestal value, using default value of 1.5")

    return pedestal



if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = '%(prog)s',
                                    description = "Script to convert bias scans (in root format) to voltage calibration constants,\
                                                   using PYROOT")
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

    bias_scan_path = os.path.abspath(args.bias_scan)
    bias_scan_directory =  os.path.dirname(bias_scan_path)

    if args.vref is None:
        if ( os.path.exists(f"{bias_scan_directory}/pedestal.root") or os.path.exists(f"{bias_scan_directory}/pedestals.root") ):
            vref = get_vref(bias_scan_directory)
        else:
            logging.warning("Could not find \"pedestal(s).root\" in the same folder. Using default value for reference voltage of 1.5V")
            vref = 1.5
    else:
        vref = args.vref

    vc = ROOT.mattak.VoltageCalibration(args.bias_scan, vref)

    # Jump into the directory in which to store the calibration
    if args.destination_folder is not None:
        os.chdir(args.destination_folder)
    else:
        # save in same directory as bias scan
        os.chdir(bias_scan_directory)

    vc.saveFitCoeffsInFile()
    # probably not necessary but there have been pyroot memory leak issues in the past
    del vc
