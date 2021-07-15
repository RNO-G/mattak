#ifndef _MATTAK_CONVERTER_H
#define _MATTAK_CONVERTER_H

namespace mattak
{
  namespace convert
  {

    int convertWaveformFile(const char * infile, const char *outfile, const char * treename =0, int station_override=-1); 
    int convertWaveformFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename =0, int station_override=-1); 
    int convertWaveformDir(const char * dir, const char * outfile, const char * treename =0, int station_override=-1); 


  }
}

#endif
