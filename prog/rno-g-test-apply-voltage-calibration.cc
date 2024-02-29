#include "mattak/VoltageCalibration.h" 
#include "mattak/Dataset.h" 
#include <stdio.h>
#include "TCanvas.h"
#include <stdlib.h>



int main (int nargs, char ** args) 
{
  if (nargs < 4) 
  {
    fprintf(stderr,"usage: rno-g-test-apply-voltage-calibration station run vcfile\n"); 
    return 1; 
  }

  mattak::DatasetOptions opt; 
  opt.verbose = true; 
  mattak::VoltageCalibration vc(args[3]); 
  opt.calib = &vc; 


  mattak::Dataset d(atoi(args[1]), atoi(args[2]), opt); 

  TCanvas c("c","c",1920,1080); 
  d.calibrated(true)->drawWaveforms(mattak::WaveformPlotOptions(), &c); 
  c.SaveAs("vc0.gif+100"); 
  c.Clear();
  d.calibrated(true)->drawWaveforms(mattak::WaveformPlotOptions(), &c); 
  c.SaveAs("vc.gif+100"); 
  d.setEntry(1); 
  c.Clear();
  d.calibrated(true)->drawWaveforms(mattak::WaveformPlotOptions(), &c); 
  c.SaveAs("vc.gif+100"); 
  c.Clear();
  d.calibrated(true)->drawWaveforms(mattak::WaveformPlotOptions(), &c); 
  c.SaveAs("vc.gif+100"); 


}
