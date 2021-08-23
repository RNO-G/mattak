#include <iostream> 
#include "mattak/Waveforms.h" 
#include "mattak/Header.h" 
#include "mattak/Pedestals.h" 
#include "mattak/DAQStatus.h" 
#include <algorithm>
#include <random> 
#include <vector> 

#include <unistd.h>
#include <stdlib.h>
#include "TFile.h" 
#include "TTree.h" 
#include "TRandom.h" 
#include "TRandom3.h" 



int main (int nargs, char ** args) 
{
  if (nargs < 5) 
  {
    std::cout << " Usage: rno-g-combine output.root waveforms.root headers.root daqstatus.root [fraction=1]" << std::endl; 
    return 1; 
  }

  double frac = 1; 
  if (nargs > 5) frac = atof(args[5]); 
  if (frac > 1 ) frac = 1; 
  if (frac <= 0) 
  {
    std::cerr << "Keep fraction 0 or smaller, no output" << std::endl; 
    return 1; 
  }

  TFile of(args[1],"RECREATE","CombinedFile",ROOT::CompressionSettings(ROOT::kLZMA,5)); 
  mattak::Waveforms * wf = 0;
  mattak::Header * hd = 0;
  mattak::DAQStatus * ds = 0;

  TTree * out = new TTree("combined","Combined RNO-G Tree"); 
  out->Branch("waveforms",&wf); 
  out->Branch("header",&hd); 
  out->Branch("daqstatus",&ds); 


  TFile wf_f(args[2]); 
  TTree * wfs = (TTree*) wf_f.Get("wf"); 
  wfs->SetBranchAddress("wf",&wf); 

  int nevents = wfs->GetEntries(); 

  TFile hd_f(args[3]); 
  TTree * hds = (TTree*) hd_f.Get("hdr"); 
  hds->SetBranchAddress("hdr",&hd); 
  hds->GetEntry(0); 

  TFile ds_f(args[4]); 
  TTree * dss = (TTree*) ds_f.Get("ds"); 
  dss->SetBranchAddress("ds",&ds); 
  dss->BuildIndex("readout_time_radiant"); 


  std::vector<int> entries; 
  if (frac < 1) 
  {
    TRandom3 r(hd->station_number * 1e8+hd->run_number); 
    entries.reserve(frac*nevents + 3*sqrt(frac*nevents)); 
    int i = 0;
    double invfrac = -log(1-frac); 
    while ( i < nevents) 
    {
      int I = 1 + floor(log(r.Rndm()) * invfrac); 
      i = (i+I); 
      if (i < nevents) 
      {
        entries.push_back(i); 
      }
      else if (!entries.size()) // if we got nothing, try again . This breaks reproducibility but only for small numbers of events, I think! 
      {
        i = 0; 
      } 
    }
  }


  int N = frac < 1 ? entries.size() : nevents; 

  for (int i = 0; i < N; i++) 
  {
    int entry = frac < 1 ? entries[i] : i; 

    wfs->GetEntry(entry); 
    hds->GetEntry(entry); 
    int status_entry = dss->GetEntryNumberWithBestIndex(hd->readout_time); 
    if (status_entry >= 0) 
    {
      dss->GetEntry(status_entry); 
    }
    else
    {
      std::cerr << "Warning, entry " << status_entry << ": events predate daqstatus" <<std::endl ;  
    }
    of.cd(); 
    out->Fill(); 
  }


  of.Write(); 
  return 0; 
}
