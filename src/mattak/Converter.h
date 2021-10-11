#ifndef _MATTAK_CONVERTER_H
#define _MATTAK_CONVERTER_H

namespace mattak
{
  namespace convert
  {
#ifdef LIBRNO_G_SUPPORT

    int convertWaveformFile(const char * infile, const char *outfile, const char * treename =0, int station_override=-1); 
    int convertWaveformFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename =0, int station_override=-1); 
    int convertWaveformDir(const char * dir, const char * outfile, const char * treename =0, int station_override=-1); 

    int convertHeaderFile(const char * infile, const char *outfile, const char * treename =0, int station_override=-1); 
    int convertHeaderFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename =0, int station_override=-1); 
    int convertHeaderDir(const char * dir, const char * outfile, const char * treename =0, int station_override=-1); 

    int convertDAQStatusFile(const char * infile, const char *outfile, const char * treename =0, int station_override=-1); 
    int convertDAQStatusFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename =0, int station_override=-1); 
    int convertDAQStatusDir(const char * dir, const char * outfile, const char * treename =0, int station_override=-1); 

    int convertPedestalFile(const char * infile, const char *outfile, const char * treename =0, int station_override=-1); 
    int convertPedestalFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename =0, int station_override=-1); 
    int convertPedestalDir(const char * dir, const char * outfile, const char * treename =0, int station_override=-1); 

#endif
    int makeRunInfo(const char * auxdir, const char * outfile, int station_override = -1, int run_override =-1); 
  }
}

#endif
