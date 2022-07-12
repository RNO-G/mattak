#include <sys/types.h> 
#include <dirent.h> 
#include <sys/stat.h> 
#include <unistd.h> 
#include <algorithm> 
#include <vector>
#include "TTree.h" 
#include "TFile.h" 
#include "mattak/Waveforms.h"
#include "mattak/Converter.h"
#include "mattak/Header.h"
#include "mattak/DAQStatus.h"
#include "mattak/Pedestals.h"
#include "mattak/RunInfo.h"


#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h" 

template <typename T> const char * getName() { return "unnamed"; } 
template<> const char * getName<mattak::Waveforms>() { return "waveforms"; } 
template<> const char * getName<mattak::Header>() { return "header"; } 
template<> const char * getName<mattak::DAQStatus>() { return "daqstatus"; } 
template<> const char * getName<mattak::Pedestals>() { return "pedestals"; } 



template <typename Traw, int(*ReaderFn)(rno_g_file_handle_t, Traw*), typename Troot>
static int convert_impl(int N, const char ** infiles, const char * outfile, const char * treename, int station)

{
  TFile * f = 0; 
  TTree *t = 0; 
  Troot * b = 0; 

  int nprocessed = 0; 
  for (int i = 0; i < N; i++) 
  {

    rno_g_file_handle h; 
    if (0 == rno_g_init_handle(&h, infiles[i], "r")) 
    {
      Traw raw; 
      while (ReaderFn(h, &raw) > 0) 
      {
        if (!f) 
        {
          f = new TFile(outfile,"RECREATE"); 
          if (!treename) treename = getName<Troot>(); 
          t = new TTree(treename,treename); 
          b = new Troot; 
          t->Branch(treename, &b); 
        }

        nprocessed++; 
        b = new (b) Troot(&raw); 
        if (station > 0) b->station_number = station; 
        t->Fill(); 
      }
      rno_g_close_handle(&h); 
    }

  }

  if (f)
  {
    f->Write(); 
    delete f;
    f = 0; 
  }

  return nprocessed; 
}




template <typename Traw, int (*ReaderFn)(rno_g_file_handle_t, Traw*), typename Troot>
static int convert_dir(const char * dir, const char * outfile, const char * treename, int station)
{
  std::vector<std::string> files; 
  std::vector<const char *> file_ptrs; 
  DIR * dirp = opendir(dir); 
  if (!dirp) 
  {
    fprintf(stderr,"Could not open dir %s\n", dir); 
    return 0; 
  }

  struct dirent * dent; 


  while ((dent = readdir(dirp)))
  {
    std::string fname = dir; 
    fname+="/"; 
    fname+=dent->d_name; 
    files.push_back(fname); 
  }

  closedir(dirp); 
  std::sort(files.begin(), files.end()); 

  file_ptrs.reserve(files.size()); 
  for (auto f : files) 
  {
    file_ptrs.push_back(f.c_str()); 
  }

  return convert_impl<Traw,ReaderFn,Troot>(file_ptrs.size(), &file_ptrs[0], outfile, treename,station); 
}


int mattak::convert::convertWaveformFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename,int station) 
{
  return convert_impl<rno_g_waveform_t, rno_g_waveform_read, mattak::Waveforms>(nfiles, infiles, outfile, treename,station); 
}

int mattak::convert::convertWaveformFile(const char * infile, const char * outfile, const char * treename,int station)
{
  return convert_impl<rno_g_waveform_t, rno_g_waveform_read, mattak::Waveforms>(1, &infile, outfile, treename,station); 
}

int mattak::convert::convertWaveformDir(const char * dir, const char * outfile, const char * treename,int station) 
{
  return convert_dir<rno_g_waveform_t, rno_g_waveform_read, mattak::Waveforms>(dir, outfile, treename,station); 
}

int mattak::convert::convertHeaderFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename,int station) 
{
  return convert_impl<rno_g_header_t, rno_g_header_read, mattak::Header>(nfiles, infiles, outfile, treename,station); 
}

int mattak::convert::convertHeaderFile(const char * infile, const char * outfile, const char * treename,int station)
{
  return convert_impl<rno_g_header_t, rno_g_header_read, mattak::Header>(1, &infile, outfile, treename,station); 
}

int mattak::convert::convertHeaderDir(const char * dir, const char * outfile, const char * treename,int station) 
{
  return convert_dir<rno_g_header_t, rno_g_header_read, mattak::Header>(dir, outfile, treename,station); 
}

int mattak::convert::convertDAQStatusFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename,int station) 
{
  return convert_impl<rno_g_daqstatus_t, rno_g_daqstatus_read, mattak::DAQStatus>(nfiles, infiles, outfile, treename,station); 
}

int mattak::convert::convertDAQStatusFile(const char * infile, const char * outfile, const char * treename,int station)
{
  return convert_impl<rno_g_daqstatus_t, rno_g_daqstatus_read, mattak::DAQStatus>(1, &infile, outfile, treename,station); 
}

int mattak::convert::convertDAQStatusDir(const char * dir, const char * outfile, const char * treename,int station) 
{
  return convert_dir<rno_g_daqstatus_t, rno_g_daqstatus_read, mattak::DAQStatus>(dir, outfile, treename,station); 
}

int mattak::convert::convertPedestalFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename,int station) 
{
  return convert_impl<rno_g_pedestal_t, rno_g_pedestal_read, mattak::Pedestals>(nfiles, infiles, outfile, treename,station); 
}

int mattak::convert::convertPedestalFile(const char * infile, const char * outfile, const char * treename,int station)
{
  return convert_impl<rno_g_pedestal_t, rno_g_pedestal_read, mattak::Pedestals>(1, &infile, outfile, treename,station); 
}

int mattak::convert::convertPedestalDir(const char * dir, const char * outfile, const char * treename,int station) 
{
  return convert_dir<rno_g_pedestal_t, rno_g_pedestal_read, mattak::Pedestals>(dir, outfile, treename,station); 
}


#endif

int mattak::convert::makeRunInfo(const char *auxdir, const char * outfile, int station_override, int run_override) 
{
  TFile of(outfile,"RECREATE"); 
  RunInfo * ri = new RunInfo(auxdir); 
  if (station_override > 0) 
  {
    ri->station = station_override; 
  }

  if (run_override > 0) 
  {
    ri->run = run_override; 
  }

//  ri->Dump(); 
  ri->Write("info"); 
  of.Close();
  return 0; 
}

