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
#include "TEventList.h"
#include "TRandom3.h" 



int main (int nargs, char ** args) 
{
  if (nargs < 5) 
  {
    std::cout << " Usage: rno-g-combine output.root waveforms.root headers.root daqstatus.root [runinfo.root=None] [FILTER=1]" << std::endl; 
    std::cout << "If runinfo.root is None, an empty runinfo will be used" << std::endl; 
    std::cout << "FILTER can be a number in (0,1] in which case it's interpreted as the fraction to keep, or a filename, in which case it's interpreted asd an event list, or a string like CUT:XXXXX where the part after CUT: is interpreted as a cut on the header file (e.g. CUT:trigger_info.force_trigger)." << std::endl; 
    return 1; 
  }


  bool use_file_list = false; 
  const char * cut = NULL;


  double frac = 1; 
  if (nargs > 6)
  {
    const char * filter_string = args[6]; 
    //is this a file list or a float? 
    char *endp = 0; 
    frac = std::strtod(filter_string,&endp); 
    if (*endp)
    {
      if (!strncasecmp("CUT:", filter_string, 4))
      {
        cut = filter_string+4; 
      }
      
      else
      {
        if (access(filter_string, R_OK))
        {
          std::cerr << "ERROR: " << filter_string <<
          " looks like a file list but we can't read it :(" << std::endl; 
          return 1; 
        }
        use_file_list = true; 
      }
    }
  }
  if (frac > 1 ) frac = 1; 
  if (frac <= 0 && !use_file_list && !cut) 
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
  TTree * wfs = nullptr; 
  if (wf_f) 
  {
    wfs = (TTree*) wf_f->Get("wf");
    if (!wfs) wfs = (TTree*) wf_f->Get("waveforms");
  }
  if (!wfs) 
  {
    std::cerr << "Could not open waveforms from " << args[2] << std::endl; 
    return 1; 
  }
  wfs->SetBranchAddress(wfs->GetName(),&wf); 

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

  

  TFile *hd_f = TFile::Open(args[3]); 
  TTree * hds = nullptr; 
  if (hd_f) 
  {
    hds =(TTree*) hd_f->Get("hdr"); 
    if (!hds) hds = (TTree*) hd_f->Get("header"); 
  }
  if (!hds) 
  {
    std::cerr << " Could not open headers from " << args[3] << std::endl; 
    return 1; 
  }
  hds->SetBranchAddress(hds->GetName(),&hd); 
  hds->GetEntry(0); 

  TFile *ds_f = TFile::Open(args[4]); 
  TTree * dss = nullptr; 
  if (ds_f) 
  {
    dss = (TTree*) ds_f->Get("ds"); 
    if (!dss) dss = (TTree*) ds_f->Get("daqstatus"); 

  }

  if (!dss) 
  {
    std::cerr << "Could not open status from " << args[4] << std::endl; 
    std::cerr << "Continuing without status"; 
  }
  else
  {
    dss->SetBranchAddress(dss->GetName(),&ds); 
    dss->BuildIndex("readout_time_radiant"); 
  }


  if (cut) 
  {
    TEventList evlist("evlist"); 
    hds->Draw(">>evlist",cut); 
    int N = evlist.GetN(); 

    if (N == 0) 
    {
      std::cerr << "File list is empty, no output" << std::endl; 
      return 1; 
    }
   
    entries.resize(N); 
    for (int i = 0; i < N; i++) 
    {
      entries[i] = evlist.GetList()[i]; 
    }
    //sort just in case 
    std::sort(entries.begin(), entries.end()); 
  }

  if (!use_file_list && !cut && frac < 1) 
  {
    TRandom3 r(hd->station_number * 1e8+hd->run_number); 
    entries.reserve(frac*nevents + 3*sqrt(frac*nevents)); 
    int i = 0;
    while ( i < nevents) 
    {
      int I = 1 + floor(-log(r.Rndm()) / frac); 
      if (I < 1) I = 1; 
      i = (i+I); 
      printf("%d %d\n",i,I); 
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


  int N = (use_file_list || cut || frac < 1) ? entries.size() : nevents; 

  for (int i = 0; i < N; i++) 
  {
    int entry = (frac < 1  || use_file_list || cut) ? entries[i] : i; 

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

  if (nargs < 6 || !strcasecmp(args[5],"None"))
  {

    std::cerr << "Not using a runinfo"  << std::endl; 
  }
  else
  {
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
  }


  of.Write(); 
  return 0; 
}
