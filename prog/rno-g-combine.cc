#include <iostream> 
#include "mattak/Waveforms.h" 
#include "mattak/Header.h" 
#include "mattak/Pedestals.h" 
#include "mattak/DAQStatus.h" 
#include "mattak/RunInfo.h" 
#include <algorithm>
#include <random> 
#include <vector> 
#include <cstdlib> 
#include <string> 
#include <fstream> 
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
    std::cout << " Usage: rno-g-combine output.root waveforms.root headers.root daqstatus.root runinfo.root [fraction=1, or a filelist of events]" << std::endl; 
    return 1; 
  }


  bool use_file_list = false; 


  double frac = 1; 
  if (nargs > 6)
  {
    const char * fraction_or_filelist = args[6]; 
    //is this a file list or a float? 
    char *endp = 0; 
    frac = std::strtod(fraction_or_filelist,&endp); 
    if (endp)
    {
      if (access(fraction_or_filelist, R_OK))
      {
        std::cerr << "ERROR: " << fraction_or_filelist <<
        " looks like a file list but we can't read it :(" << std::endl; 
        return 1; 
      }
      use_file_list = true; 
    }
  }
  if (frac > 1 ) frac = 1; 
  if (frac <= 0 && !use_file_list) 
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


  TFile *wf_f = TFile::Open(args[2]); 
  TTree * wfs = wf_f ? (TTree*) wf_f->Get("wf") : 0; 
  if (!wfs) 
  {
    std::cerr << "Could not open waveforms from " << args[2] << std::endl; 
    return 1; 
  }
  wfs->SetBranchAddress("wf",&wf); 

  int nevents = wfs->GetEntries(); 
  std::vector<int> entries; 

  if (use_file_list) 
  {
    std::ifstream list(args[6]); 
    std::string line; 
    while (std::getline(list,line))
    {
      char * endp = 0; 
      int entry = std::strtol(line.c_str(), &endp, 10); 
      if (endp == line.c_str()) //not a number
        continue; 

      entries.push_back(entry); 
    }
    if (entries.size() == 0) 
    {
      std::cerr << "File list is empty, no output" << std::endl; 
      return 1; 
    }
    //sort just in case 
    std::sort(entries.begin(), entries.end()); 
  }

  int N = use_file_list ? entries.size() : 
               frac < 1 ? nevents * frac :
               nevents; 

  //TODO: std::sample is probalby better here
  if (!use_file_list && frac < 1) 
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

  TFile *hd_f = TFile::Open(args[3]); 
  TTree * hds = hd_f ? (TTree*) hd_f->Get("hdr") : 0; 
  if (!hds) 
  {
    std::cerr << " Could not open headers from " << args[3] << std::endl; 
    return 1; 
  }
  hds->SetBranchAddress("hdr",&hd); 
  hds->GetEntry(0); 

  TFile *ds_f = TFile::Open(args[4]); 
  TTree * dss = ds_f ? (TTree*) ds_f->Get("ds") : 0; 

  if (!dss) 
  {
    std::cerr << "Could not open status from " << args[4] << std::endl; 
    std::cerr << "Continuing without status"; 
  }
  else
  {
    dss->SetBranchAddress("ds",&ds); 
    dss->BuildIndex("readout_time_radiant"); 
  }


  TFile * ri_f = TFile::Open(args[5]); 
  mattak::RunInfo * ri = ri_f ? (mattak::RunInfo*) ri_f->Get("info") : 0; 
  if (!ri) 
  {
    std::cerr << "Could not open runinfo from " << args[5] << std::endl; 
    std::cerr << "Continuing without runinfo." <<std::endl; 
  }
  else
  {
    of.cd(); 
    ri->Write("info"); 
  }

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
    int entry = frac < 1  || use_file_list ? entries[i] : i; 

    wfs->GetEntry(entry); 
    hds->GetEntry(entry); 
    if (dss) 
    {
      int status_entry = dss->GetEntryNumberWithBestIndex(hd->readout_time); 
      if (status_entry >= 0) 
      {
        dss->GetEntry(status_entry); 
      }
      else
      {
        std::cerr << "Warning, entry " << status_entry << ": events predate daqstatus" <<std::endl ;  
      }
    }
    of.cd(); 
    out->Fill(); 
  }

  of.Write(); 
  return 0; 
}
