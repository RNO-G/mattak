#ifndef _MATTAK_CONVERTER_H
#define _MATTAK_CONVERTER_H

namespace mattak
{
  namespace convert
  {
#ifdef LIBRNO_G_SUPPORT

    /* All of these take an `update` flag. With it set, an output file which
     * already exists is extended with whatever input files are not in it yet,
     * instead of being written from scratch. That is what makes converting a
     * run which is still being taken cheap, since each pass then only pays for
     * the raw files which arrived since the last one.
     *
     * Which input files an output file holds is recorded inside it, so this is
     * safe to repeat: files already converted are skipped, and anything the
     * record does not explain (an input file which changed, one which appeared
     * before a file already converted, a mismatched entry count, an output file
     * which cannot be opened) falls back to converting from scratch.
     */

    int convertWaveformFile(const char * infile, const char *outfile, const char * treename =0, int station_override=-1, bool update=false);
    int convertWaveformFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename =0, int station_override=-1, bool update=false);
    int convertWaveformDir(const char * dir, const char * outfile, const char * treename =0, int station_override=-1, bool update=false);

    int convertHeaderFile(const char * infile, const char *outfile, const char * treename =0, int station_override=-1, bool update=false);
    int convertHeaderFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename =0, int station_override=-1, bool update=false);
    int convertHeaderDir(const char * dir, const char * outfile, const char * treename =0, int station_override=-1, bool update=false);

    int convertDAQStatusFile(const char * infile, const char *outfile, const char * treename =0, int station_override=-1, bool update=false);
    int convertDAQStatusFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename =0, int station_override=-1, bool update=false);
    int convertDAQStatusDir(const char * dir, const char * outfile, const char * treename =0, int station_override=-1, bool update=false);

    int convertPedestalFile(const char * infile, const char *outfile, const char * treename =0, int station_override=-1, bool update=false);
    int convertPedestalFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename =0, int station_override=-1, bool update=false);
    int convertPedestalDir(const char * dir, const char * outfile, const char * treename =0, int station_override=-1, bool update=false);

#endif
    int makeRunInfo(const char * auxdir, const char * outfile, int station_override = -1, int run_override =-1);
  }
}

#endif
