#include <iostream>
#include <cstring>
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
  std::cout << "Usage: rno-g-converter [--update] TYPE OUTFILE INFILE1 [INFILE2 ...]" << std::endl;
  std::cout << "   or rno-g-converter [--update] TYPE OUTFILE INDIR << std::endl" << std::endl;
  std::cout << "   or rno-g-converter runinfo OUTFILE  auxdir [stationoverride] [runoverride] << std::endl" << std::endl;
  std::cout << "   TYPE can be waveforms wfs header hd daqstatus ds pedestal ped runinfo" << std::endl;
  std::cout << "   OUTFILE is the output ROOT file " << std::endl;
  std::cout << "   INFILE is one or more input files in order " << std::endl;
  std::cout << "   INDIR is an input directory, an attempt will be made to get any number from a string and sorted. " << std::endl;
  std::cout << "   --update extends an existing OUTFILE with the inputs it does not hold yet," << std::endl;
  std::cout << "            instead of converting everything again. Inputs already converted" << std::endl;
  std::cout << "            are skipped; anything unexpected falls back to a full conversion." << std::endl;

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
  bool update = false;
  int argi = 1;

  while (argi < nargs && args[argi][0] == '-')
  {
    if (!strcmp(args[argi], "--update") || !strcmp(args[argi], "-u"))
    {
      update = true;
      argi++;
    }
    else
    {
      if (strcmp(args[argi], "--help") && strcmp(args[argi], "-h"))
      {
        std::cerr << "Unknown option: " << args[argi] << std::endl;
      }
      usage();
    }
  }

  if (nargs - argi < 3) usage();

  const char * type       = args[argi++];
  const char * outfile    = args[argi++];
  //check if outfile ends with .root
  if (!TString(outfile).EndsWith(".root"))
  {
    std::cerr << "Outfile should end with .root" << std::endl;
    return 1;
  }


  const char * firstinput = args[argi];
  int first_input_dir = isdir(firstinput);
  int nfiles = nargs - argi;
  int N = 0;

  if (!strcmp(type,"wf") || !strcmp(type,"waveforms"))
  {
    if (first_input_dir)
    {
      N = mattak::convert::convertWaveformDir(firstinput, outfile,0,0,update);
    }
    else
    {
      N = mattak::convert::convertWaveformFiles(nfiles, (const char**) (args+argi), outfile, 0, -1, update);
    }
  }
  else if (!strcmp(type,"hd") || !strcmp(type,"header"))
  {
    if (first_input_dir)
    {
      N = mattak::convert::convertHeaderDir(firstinput, outfile,0,0,update);
    }
    else
    {
      N = mattak::convert::convertHeaderFiles(nfiles, (const char**) (args+argi), outfile, 0, -1, update);
    }
  }
  else if (!strcmp(type,"ds") || !strcmp(type,"daqstatus"))
  {
    if (first_input_dir)
    {
      N = mattak::convert::convertDAQStatusDir(firstinput, outfile,0,0,update);
    }
    else
    {
      N = mattak::convert::convertDAQStatusFiles(nfiles, (const char**) (args+argi), outfile, 0, -1, update);
    }
  }
  else if (!strcmp(type,"ped") || !strcmp(type,"pedestal"))
  {
    if (first_input_dir)
    {
      N = mattak::convert::convertPedestalDir(firstinput, outfile,0,0,update);
    }
    else
    {
      N = mattak::convert::convertPedestalFiles(nfiles, (const char**) (args+argi), outfile, 0, -1, update);
    }
  }
  else if (!strcmp(type,"runinfo"))
  {
    return mattak::convert::makeRunInfo(firstinput,outfile, nfiles > 1 ? atoi(args[argi+1]) : -1, nfiles > 2 ? atoi(args[argi+2]) : -1 );
  }

  else
  {
    std::cerr << "Unkown type: " << type << std::endl;
    return 1;
  }

  printf("Processed %d entries\n", N);
  return 0;
}
