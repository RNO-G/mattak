#include "rno-g.h" 
#include <iostream> 
#include "TH1.h" 
#include "TFile.h" 

void usage() 
{
  std::cout << "Usage: rno-g-readout-elapsed OUTFILE INFILE1 [INFILE2 ...]" << std::endl; 
  exit(1); 
} 


int main (int nargs, char ** args) 
{

  if (nargs < 3) usage(); 
 
  TFile f(args[1],"CREATE"); 

  if (!f.IsOpen()) return 1; 

  TH1  * hist = new TH1I ("elapsed_hist","Elapsed Time ; ms ; number", 1000,0,200); 


  for (int i = 2; i < nargs; i++) 
  {

    rno_g_file_handle_t h; 
    rno_g_init_handle(&h, args[i], "r"); 
    rno_g_header_t hd; 

    while (rno_g_header_read(h, &hd) > 0)
    {
      hist->Fill(hd.readout_elapsed_nsecs/1e6); 
    }
  }

  hist->Write(); 
}
