#include "mattak/VoltageCalibration.h" 
#include "mattak/Dataset.h" 
#include <stdio.h>
#include "TCanvas.h"
#include <stdlib.h>



int main (int nargs, char ** args) 
{
  if (nargs < 4) 
  {
    fprintf(stderr,"usage: rno-g-test-apply-voltage-calibration station run vcfile max=10 output=1\n"); 
    return 1; 
  }

  mattak::DatasetOptions opt; 
  opt.verbose = true; 
  mattak::VoltageCalibration vc(args[3]); 
  opt.calib = &vc; 

  int station = atoi(args[1]); 
  int run = atoi(args[2]); 
  int max = nargs >= 4 ? atoi(args[4]) : 10; 
  int output = nargs >=5 ? atoi(args[5]) : 1; 

  mattak::Dataset d(station,run, opt); 

  if (max <=0 || max > d.N()) max = d.N(); 
  TCanvas c("c","c",1920,1080); 
  for (int i = 0; i <  max ; i++) 
  {
    d.setEntry(i); 
    printf("%d\n", i); 
    //do it twice, just in case
    for (int j = 0; j < 2; j++) 
    {
      auto cal = d.calibrated(j > 0); 
      auto raw = d.raw(j > 0); 

      if (output) 
      {
        cal->drawWaveforms(mattak::WaveformPlotOptions(), &c); 
        c.SaveAs(Form("vc_%d_%d.gif+50", station, run)); 
        mattak::WaveformPlotOptions optraw; 
        optraw.color = kRed +2; 
        raw->drawWaveforms(optraw, &c); 
        c.SaveAs(Form("vc_%d_%d.gif+50", station, run)); 
      }
    }
  }

  return 0; 

}
