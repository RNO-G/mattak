#include <iostream>
#include "mattak/Dataset.h"
#include <stdlib.h>
#include <strings.h>


void usage()
{

    std::cerr << "Usage: rno-g-count-events stationnum runnum[:endrunnum] [data_dir=$RNO_G_DATA]" << std::endl;
}

int main(int nargs, char ** args)
{

  if (nargs < 3)
  {
    usage();
    return 1;
  }

  int station = atoi(args[1]);
  int start_run = -1, end_run = -1;
  if (strchr(args[2],':'))
  {
    sscanf(args[2],"%d:%d", &start_run, &end_run);
  }
  else
  {
    start_run = atoi(args[2]);
    end_run = start_run;
  }

  if (start_run < 0)
  {
    usage();
    return -1;
  }

  mattak::DatasetOptions opt;
  if (nargs >3) opt.base_data_dir = args[3];

  opt.partial_skip_incomplete = false;

  for (int run = start_run; run <= end_run ; run++)
  {
    mattak::Dataset d(station,run, opt);

    if (d.N() <=0)
    {
      std::cout << "Station " << station << " Run " << run << ": Could not find dataset or dataset empty" << std::endl;
      continue;
    }

    int nforce = 0;
    int nrad  = 0;
    int nrad0 = 0;
    int nrad1 = 0;
    int nlt = 0;
    int npps = 0;

    for (int i = 0; i < d.N(); i++)
    {
      d.setEntry(i);
      if (d.header()->trigger_info.force_trigger) nforce++;
      if (d.header()->trigger_info.lt_trigger) nlt++;
      if (d.header()->trigger_info.pps_trigger) npps++;
      if (d.header()->trigger_info.radiant_trigger) nrad++;
      if (d.header()->trigger_info.which_radiant_trigger == 0) nrad0++;
      if (d.header()->trigger_info.which_radiant_trigger == 1) nrad1++;
    }

    std::cout << "Station " << station << " Run " << run << ":" << std::endl;
    std::cout << "    FORCE_TRIGGERS:   " << nforce << std::endl;
    std::cout << "    LT_TRIGGERS:      " << nlt << std::endl;
    std::cout << "    RADIANT_TRIGGERS: " << nrad<< std::endl;
    std::cout << "                 RF0: " << nrad0 << std::endl;
    std::cout << "                 RF1: " << nrad1 << std::endl;
    std::cout << "    PPS_TRIGGERS:     " << npps << std::endl;
    std::cout << "    TOTAL_TRIGGERS:   " << d.N() << std::endl;
  }

  return 0;

}
