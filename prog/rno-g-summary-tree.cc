#include <iostream> 
#include "mattak/Waveforms.h" 
#include "mattak/Header.h" 
#include "mattak/DAQStatus.h" 
#include "TFile.h" 
#include "TTree.h" 
#include "TMath.h" 


int main(int nargs, char ** args) 
{
  if (nargs < 4)
  {
    std::cout << " Usage: rno-g-summary-tree output.root waveforms.root headers.root" <<std::endl; 
    return 1; 
  }

  TFile of(args[1],"RECREATE"); 

  TTree * sum = new TTree("sum","Summary tree"); 
  
  double Vpp[24]; 
  double rms[24]; 
  double dt = -1; 
  mattak::Waveforms * wf = 0;
  mattak::Header * hd = 0;

  sum->Branch("Vpp", Vpp,"Vpp[24]/D"); 
  sum->Branch("rms", rms,"rms[24]/D"); 
  sum->Branch("dt", &dt); 

  TFile wf_f(args[2]); 
  TTree * wfs = (TTree*) wf_f.Get("wf"); 
  wfs->SetBranchAddress("wf",&wf); 

  TFile hd_f(args[3]); 
  TTree * hds = (TTree*) hd_f.Get("hdr"); 
  hds->SetBranchAddress("hdr",&hd); 

  sum->Branch("hdr",&hd); 

  double last_t = 0; 
  for (int i = 0; i < hds->GetEntries(); i++) 
  {

    hds->GetEntry(i); 
    wfs->GetEntry(i); 

    if (i > 1) 
    {
      dt = hd->trigger_time - last_t; 
    }
    last_t= hd->trigger_time; 

    for (int ichan = 0; ichan <24; ichan++) 
    {
      double vmax = 0; 
      double vmin = 0; 

      //use more accurate way of computing RMS 
      double oldM1 = 0; 
      double M1 = 0, S1 = 0; 
      for (int isamp = 0; isamp  < wf->buffer_length; isamp++) 
      {
        double val = wf->radiant_data[ichan][isamp]; 
        M1 = (val-oldM1) / (wf->buffer_length+1); 
        S1 += (val - oldM1) * (val - M1); 
        oldM1 = M1; 
        if (val > vmax) vmax = val; 
        if (val < vmin) vmin = val; 
      }

      rms[ichan] = sqrt(S1 / (wf->buffer_length-1)); 
      Vpp[ichan] = vmax-vmin; 
    }

    of.cd(); 
    sum->Fill(); 
  }

  of.Write(); 

  return 0; 

}
