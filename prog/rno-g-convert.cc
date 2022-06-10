#include <iostream> 
#include "mattak/Waveforms.h" 
#include "mattak/Header.h" 
#include "mattak/Pedestals.h" 
#include "mattak/DAQStatus.h" 
#include "mattak/Converter.h" 
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void usage() 
{
  std::cout << "Usage: rno-g-converter TYPE OUTFILE INFILE1 [INFILE2 ...]" << std::endl; 
  std::cout << "   or rno-g-converter TYPE OUTFILE INDIR << std::endl" << std::endl; 
  std::cout << "   or rno-g-converter runinfo OUTFILE  auxdir [stationoverride] [runoverride] << std::endl" << std::endl; 
  std::cout << "   TYPE can be waveforms wfs header hd daqstatus ds pedestal ped runinfo" << std::endl; 
  std::cout << "   OUTFILE is the output ROOT file " << std::endl;
  std::cout << "   INFILE is one or more input files in order " << std::endl; 
  std::cout << "   INDIR is an input directory, an attempt will be made to get any number from a string and sorted. " << std::endl; 

  exit(1); 

} 

int isdir(const char * f) 
{
  struct stat st; 
  if (stat(f,   &st) != 0) 
  {
    return 0; 
  }
  return S_ISDIR(st.st_mode); 
}



int main (int nargs, char ** args) 
{

  if (nargs < 4) usage(); 
  const char * type       = args[1]; 
  const char * outfile    = args[2]; 
  //check if outfile ends with .root 
  if (!TString(outfile).EndsWith(".root"))
  {
    std::cerr << "Outfile should end with .root" << std::endl; 
    return 1; 
  }

  
  const char * firstinput = args[3]; 
  int first_input_dir = isdir(firstinput); 
  int N = 0;

  if (!strcmp(type,"wf") || !strcmp(type,"waveforms"))
  {
    if (first_input_dir) 
    {
      N = mattak::convert::convertWaveformDir(firstinput, outfile,0,0); 
    }
    else 
    {
      N = mattak::convert::convertWaveformFiles(nargs-3, (const char**) (args+3), outfile); 
    }
  }
  else if (!strcmp(type,"hd") || !strcmp(type,"header"))
  {
    if (first_input_dir) 
    {
      N = mattak::convert::convertHeaderDir(firstinput, outfile,0,0); 
    }
    else 
    {
      N = mattak::convert::convertHeaderFiles(nargs-3, (const char**) (args+3), outfile); 
    }
  }
  else if (!strcmp(type,"ds") || !strcmp(type,"daqstatus"))
  {
    if (first_input_dir) 
    {
      N = mattak::convert::convertDAQStatusDir(firstinput, outfile,0,0); 
    }
    else 
    {
      N = mattak::convert::convertDAQStatusFiles(nargs-3, (const char**) (args+3), outfile); 
    }
  }
  else if (!strcmp(type,"ped") || !strcmp(type,"pedestal"))
  {
    if (first_input_dir) 
    {
      N = mattak::convert::convertPedestalDir(firstinput, outfile,0,0); 
    }
    else 
    {
      N = mattak::convert::convertPedestalFiles(nargs-3, (const char**) (args+3), outfile); 
    }
  }
  else if (!strcmp(type,"runinfo"))
  {
    return mattak::convert::makeRunInfo(firstinput,outfile, nargs > 4 ? atoi(args[4]) : -1, nargs > 5 ? atoi(args[5]) : -1 ); 
  }
 
  else
  {
    std::cerr << "Unkown type: " << type << std::endl; 
    return 1; 
  }

  printf("Processed %d entries\n", N); 
  return 0; 
}





