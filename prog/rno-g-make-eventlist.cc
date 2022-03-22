/** 
 * This just prints out all the events in a combined.root, useful for remaking files on the server
 */

#include <iostream> 
#include "TFile.h" 
#include "TTree.h" 
#include "mattak/Header.h" 

int main (int nargs, char ** args) 
{
  if (nargs < 2) 
  {
    std::cerr << "Usage: rno-g-make-eventlist combined.root [treename=combined]" << std::endl;
    return 1; 
  }


  TFile * f = TFile::Open(args[1]); 
  if (!f) 
  {
    std::cerr << "Could not open " << args[1] << std::endl; 
    return 1; 
  }

  const char * treename = nargs >2 ? args[2] : "combined"; 
  TTree * t = (TTree*) f->Get(treename); 
  if (!t) 
  {
    std::cerr << "Could not open tree " << treename << " in " << args[1] << std::endl; 
    return 1; 
  }
     
  mattak::Header * h = 0; 
//  t->SetBranchStatus("*",0); 
  t->SetBranchStatus("waveforms",0); 
  t->SetBranchAddress("header",&h); 
  
  int iprocessed= 0; 
  for (int i = 0; i < t->GetEntries(); i++) 
  {
    t->GetEntry(i); 
    std::cout << h->event_number << std::endl; 
    iprocessed++; 
  }

  return iprocessed == 0; 
}
