#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <string>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TError.h"
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


/* Incremental (update) conversion.
 *
 * A run is converted over and over while it is being taken, each time gaining
 * only the handful of raw files which arrived since. Waveforms and headers carry
 * the run's sequential event number, so an existing output file already says
 * where to continue: everything up to the event number of its last entry is
 * converted. Raw events at or below that are skipped and the rest appended.
 *
 * That needs no bookkeeping beside the data, and it recovers by itself from an
 * interrupted conversion: ROOT commits a tree by rewriting its key, so a
 * converter which is killed part way leaves the tree at its last committed
 * length, and the next pass simply resumes from the event number it finds there.
 * It also guarantees the result stays ordered, since an event is only ever
 * appended if it comes after everything already in the tree.
 */

/** The event number of a raw record, or -1 for types which do not have one
 * (pedestals, daqstatus) and are therefore always converted in full. */
template <typename Traw> static long long rawEventNumber(const Traw &) { return -1; }
template<> long long rawEventNumber(const rno_g_waveform_t & raw) { return raw.event_number; }
template<> long long rawEventNumber(const rno_g_header_t & raw) { return raw.event_number; }

template <typename Troot> static long long eventNumberOf(const Troot *) { return -1; }
template<> long long eventNumberOf(const mattak::Waveforms * wf) { return wf->event_number; }
template<> long long eventNumberOf(const mattak::Header * hd) { return hd->event_number; }


/** The event number of the last entry already in outfile, or -1 if it cannot be
 * appended to (no such file, no such tree, empty, or a type without event
 * numbers) and has to be written from scratch.
 */
template <typename Troot>
static long long lastConvertedEvent(const char * outfile, const char * treename)
{
  if (access(outfile, R_OK)) return -1;

  // Read-only on purpose: opening a TFile for writing rewrites its header, and
  // most passes turn out to have nothing to add at all.
  TFile * f = TFile::Open(outfile, "READ");
  if (!f || f->IsZombie())
  {
    ::Error("mattak::convert", "Could not open %s, converting from scratch", outfile);
    delete f;
    return -1;
  }

  long long last = -1;
  TTree * t = (TTree *) f->Get(treename);

  if (t && t->GetEntries() > 0)
  {
    Troot * b = 0;
    long long first = -1;
    t->SetBranchAddress(treename, &b);

    if (t->GetEntry(0) > 0) first = eventNumberOf<Troot>(b);
    if (t->GetEntry(t->GetEntries() - 1) > 0) last = eventNumberOf<Troot>(b);

    // Event numbers are sequential, so a complete tree holds every one of them
    // between its first and its last. If any are missing -- a raw file which
    // could not be read on an earlier pass, one which turned up only after later
    // events had been converted, or one which was truncated -- then they all sit
    // below `last`, where appending would never reach them again, and only
    // converting from scratch brings them back.
    if ((first >= 0 && last >= 0 && t->GetEntries() != last - first + 1) || first > 10)
    {
      ::Error("mattak::convert",
              "%s holds %lld entries but spans events %lld..%lld, so some are missing; converting from scratch",
              outfile, (long long) t->GetEntries(), first, last);
      last = -1;
    }
  }

  delete f;
  return last;
}


template <typename Traw, int(*ReaderFn)(rno_g_file_handle_t, Traw*), typename Troot>
static int convert_impl(
    int N, const char ** infiles, const char * outfile,
    const char * treename, int station, bool update)
{
  if (!treename) treename = getName<Troot>();

  TFile * f = 0;
  TTree * t = 0;
  Troot * b = 0;

  const long long last = update ? lastConvertedEvent<Troot>(outfile, treename) : -1;
  const bool appending = last >= 0;

  int nprocessed = 0;
  for (int i = 0; i < N; i++)
  {
    rno_g_file_handle h;
    if (0 == rno_g_init_handle(&h, infiles[i], "r"))
    {
      Traw raw;
      while (ReaderFn(h, &raw) > 0)
      {
        if (appending && rawEventNumber<Traw>(raw) <= last) continue;

        // Opened only once there is something to write, so that a pass which
        // finds nothing new leaves the file untouched.
        if (!f)
        {
          if (appending)
          {
            f = TFile::Open(outfile, "UPDATE");
            t = f && !f->IsZombie() ? (TTree *) f->Get(treename) : 0;
            if (!t)
            {
              ::Error("mattak::convert", "Could not reopen %s to append to", outfile);
              delete f;
              rno_g_close_handle(&h);
              return 0;
            }
            b = new Troot;
            t->SetBranchAddress(treename, &b);
          }
          else
          {
            f = new TFile(outfile, "RECREATE");
            t = new TTree(treename, treename);
            b = new Troot;
            t->Branch(treename, &b);
          }
        }

        nprocessed++;
        b = new (b) Troot(&raw);
        if (station > 0) b->station_number = station;
        t->Fill();
      }
      rno_g_close_handle(&h);
    }
    else
    {
      ::Error("mattak::convert", "Could not open %s", infiles[i]);
    }
  }

  if (f)
  {
    if (appending) t->Write("", TObject::kOverwrite);
    else f->Write();
    delete f;
    f = 0;
  }

  return nprocessed;
}


template <typename Traw, int (*ReaderFn)(rno_g_file_handle_t, Traw*), typename Troot>
static int convert_dir(const char * dir, const char * outfile, const char * treename, int station, bool update)
{
  std::vector<std::string> files;
  std::vector<const char *> file_ptrs;
  DIR * dirp = opendir(dir);
  if (!dirp)
  {
    ::Error("mattak::convert::convert_dir", "Could not open dir %s", dir);
    return 0;
  }

  struct dirent * dent;

  while ((dent = readdir(dirp)))
  {
    std::string fname = dir;
    fname += "/";
    fname += dent->d_name;
    files.push_back(fname);
  }

  closedir(dirp);
  std::sort(files.begin(), files.end());

  file_ptrs.reserve(files.size());
  // by reference: a copy would leave file_ptrs full of dangling pointers
  for (const auto & f : files)
  {
    file_ptrs.push_back(f.c_str());
  }

  return convert_impl<Traw, ReaderFn, Troot>(file_ptrs.size(), &file_ptrs[0], outfile, treename, station, update);
}


int mattak::convert::convertWaveformFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename, int station, bool update)
{
  return convert_impl<rno_g_waveform_t, rno_g_waveform_read, mattak::Waveforms>(nfiles, infiles, outfile, treename, station, update);
}

int mattak::convert::convertWaveformFile(const char * infile, const char * outfile, const char * treename, int station, bool update)
{
  return convert_impl<rno_g_waveform_t, rno_g_waveform_read, mattak::Waveforms>(1, &infile, outfile, treename, station, update);
}

int mattak::convert::convertWaveformDir(const char * dir, const char * outfile, const char * treename, int station, bool update)
{
  return convert_dir<rno_g_waveform_t, rno_g_waveform_read, mattak::Waveforms>(dir, outfile, treename, station, update);
}

int mattak::convert::convertHeaderFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename, int station, bool update)
{
  return convert_impl<rno_g_header_t, rno_g_header_read, mattak::Header>(nfiles, infiles, outfile, treename, station, update);
}

int mattak::convert::convertHeaderFile(const char * infile, const char * outfile, const char * treename, int station, bool update)
{
  return convert_impl<rno_g_header_t, rno_g_header_read, mattak::Header>(1, &infile, outfile, treename, station, update);
}

int mattak::convert::convertHeaderDir(const char * dir, const char * outfile, const char * treename, int station, bool update)
{
  return convert_dir<rno_g_header_t, rno_g_header_read, mattak::Header>(dir, outfile, treename, station, update);
}

int mattak::convert::convertDAQStatusFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename, int station, bool update)
{
  return convert_impl<rno_g_daqstatus_t, rno_g_daqstatus_read, mattak::DAQStatus>(nfiles, infiles, outfile, treename, station, update);
}

int mattak::convert::convertDAQStatusFile(const char * infile, const char * outfile, const char * treename, int station, bool update)
{
  return convert_impl<rno_g_daqstatus_t, rno_g_daqstatus_read, mattak::DAQStatus>(1, &infile, outfile, treename, station, update);
}

int mattak::convert::convertDAQStatusDir(const char * dir, const char * outfile, const char * treename, int station, bool update)
{
  return convert_dir<rno_g_daqstatus_t, rno_g_daqstatus_read, mattak::DAQStatus>(dir, outfile, treename, station, update);
}

int mattak::convert::convertPedestalFiles(int nfiles, const char ** infiles, const char * outfile, const char * treename, int station, bool update)
{
  return convert_impl<rno_g_pedestal_t, rno_g_pedestal_read, mattak::Pedestals>(nfiles, infiles, outfile, treename, station, update);
}

int mattak::convert::convertPedestalFile(const char * infile, const char * outfile, const char * treename, int station, bool update)
{
  return convert_impl<rno_g_pedestal_t, rno_g_pedestal_read, mattak::Pedestals>(1, &infile, outfile, treename, station, update);
}

int mattak::convert::convertPedestalDir(const char * dir, const char * outfile, const char * treename, int station, bool update)
{
  return convert_dir<rno_g_pedestal_t, rno_g_pedestal_read, mattak::Pedestals>(dir, outfile, treename, station, update);
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

  //ri->Dump();
  ri->Write("info");
  of.Close();
  return 0;
}
