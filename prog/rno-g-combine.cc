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
  int N = frac < 1 ? nevents * frac : nevents; 
  std::vector<int> entries; 

  //TODO: std::sample is probalby better here
  if (frac < 1) 
  {
    std::random_device rd; 
    std::mt19937 g(rd()); 

    entries.reserve(nevents); 
    for (int i = 0; i < nevents; i++) 
    {
      entries.push_back(i); 
    }
    std::shuffle(entries.begin(), entries.end(),g); 
    entries.resize(N); 
    std::sort(entries.begin(), entries.end()); 
  }

  TFile hd_f(args[3]); 
  TTree * hds = (TTree*) hd_f.Get("hdr"); 
  hds->SetBranchAddress("hdr",&hd); 

  TFile ds_f(args[4]); 
  TTree * dss = (TTree*) ds_f.Get("ds"); 
  dss->SetBranchAddress("ds",&ds); 
  dss->BuildIndex("readout_time_radiant"); 


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
