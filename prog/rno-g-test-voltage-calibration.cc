#include "mattak/VoltageCalibration.h"


int main(int nargs, char ** args)
{
  const char * filename = nargs > 1 ? args[1]  : "bias_scan.dat.gz";

  mattak::VoltageCalibration vc(filename);
  vc.saveFitCoeffsInFile();

}
