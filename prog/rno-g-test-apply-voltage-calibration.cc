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

  int station = atoi(args[1]); 
  int run = atoi(args[2]); 

  mattak::Dataset d(station,run, opt); 

  TCanvas c("c","c",1920,1080); 
  for (int i = 0; i <  (d.N() < 10 ? d.N() : 10) ; i++) 
  {
    d.setEntry(i); 
    printf("%d\n", i); 
    //do it twice, just in case
    for (int j = 0; j < 2; j++) 
    {
      d.calibrated(j > 0)->drawWaveforms(mattak::WaveformPlotOptions(), &c); 
      c.SaveAs(Form("vc_%d_%d.gif+50", station, run)); 
      mattak::WaveformPlotOptions optraw; 
      optraw.color = kRed +2; 
      d.raw(j > 0)->drawWaveforms(optraw, &c); 
      c.SaveAs(Form("vc_%d_%d.gif+50", station, run)); 
    }
  }


}
