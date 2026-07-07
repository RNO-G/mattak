#ifndef _MATTAK_DATASET_H
#define _MATTAK_DATASET_H


#include "TFile.h"
#include "TTree.h"

#include "mattak/Header.h"
#include "mattak/Waveforms.h"
#include "mattak/DAQStatus.h"
#include "mattak/Pedestals.h"
#include "mattak/RunInfo.h"

namespace mattak
{
  class VoltageCalibration;


  /** Options controlling reading of the dataset
   * This can be passed to a constructor or each time you load a run or directory.
   * */
  struct DatasetOptions
  {
    /** The voltage calibration to use. Will eventually do something sensible if nullptr (though now does nothing)*/
    const VoltageCalibration * calib = nullptr;

    /** The base data directory, used when loading runs by station/run number. If empty, will use getenv */
    std::string base_data_dir = "";

    /** Controls behavior for skipping incomplete on partial datasets where the entire
     * metadata is available.
     *
     * If true, will only return complete events when indexing. If false, it will loop over all events, but the waveforms may not be available for all.
     *
     * */
    bool partial_skip_incomplete = true;


    /** Controls which files get read
     * If empty, will prefer full datasets followed by combined.root
     *
     * You can e.g. set to "combined" to always go for combined.root or to some
     * other string if you want someother file that will be treated as a combined.root
     *
     * */
    std::string file_preference = "";


    /** If true, emit **debug-level** messages while loading.
     *
     * All Dataset messages go through ROOT's message system (TError.h):
     * debug messages are emitted via Info, recoverable problems via Warning
     * and load failures via Error. They can therefore be suppressed or
     * redirected globally, e.g. gErrorIgnoreLevel = kError; (also from
     * PyROOT: ROOT.gErrorIgnoreLevel = ROOT.kError) or SetErrorHandler().
     * */
    bool verbose = false;
  };


  /** Lazy-loading interface to one run of RNO-G data.
   *
   * A dataset can be loaded from two kinds of sources:
   *
   *  - A run directory (loadDir / loadRun), which is either
   *      * a full dataset: waveforms.root + headers.root + daqstatus.root
   *        (and optionally pedestal.root / runinfo.root), or
   *      * a partial dataset: combined.root containing only a subset of events,
   *        optionally accompanied by the full headers.root / daqstatus.root.
   *        In this case opt.partial_skip_incomplete decides
   *        whether to index only the complete events in combined.root (true)
   *        or all events in headers.root (false; waveforms may then be
   *        unavailable for some entries).
   *
   *  - A single combined file (loadCombinedFile), in which case only the
   *    information inside that file is read and partial_skip_incomplete is
   *    forced to true.
   *
   * All load methods return 0 on success and leave the dataset unloaded on
   * failure. Getters like raw()/header()/status() lazily read the entry
   * selected with setEntry() and return nullptr if unavailable.
   *
   * Canonical usage, by station/run number (resolved as
   * base_data_dir/station<S>/run<R>, where base_data_dir defaults to
   * $RNO_G_ROOT_DATA or $RNO_G_DATA):
   *
   *   mattak::Dataset d(23, 1144);
   *   for (int i = 0; i < d.N(); i++)
   *   {
   *     d.setEntry(i);
   *     mattak::Header * hdr = d.header();
   *     mattak::Waveforms * wfs = d.raw(); // nullptr if this event has no waveforms
   *   }
   *
   * or from an explicit path to a run directory:
   *
   *   mattak::Dataset d;
   *   if (d.loadDir("/data/handcarry22/station23/run1144"))
   *   {
   *     // nonzero return: no loadable run files in that directory
   *   }
   *
   * To load a specific (combined-style) file inside the run directory instead
   * of the default lookup, set opt.file_preference to its name without the
   * ".root" extension, e.g. to force the partial combined.root even when the
   * full dataset is present:
   *
   *   mattak::DatasetOptions opt;
   *   opt.file_preference = "combined";  // reads <run dir>/combined.root
   *   mattak::Dataset d(23, 1144, opt);
   */
  class Dataset
  {

    public:
      /** Create an empty dataset; call one of the load methods afterwards. */
      Dataset(const DatasetOptions & opt = DatasetOptions());

      /** Create a dataset and immediately load a run (see loadRun). */
      Dataset(int station, int run, const DatasetOptions & opt = DatasetOptions());

      // Dataset owns its TFiles and cached objects, so it must not be copied
      Dataset(const Dataset &) = delete;
      Dataset & operator=(const Dataset &) = delete;

      virtual ~Dataset() { unload() ; }

      /** Load opt.base_data_dir/station<station>/run<run> (see loadDir). */
      int loadRun(int station, int run);
      int loadRun(int station, int run, const DatasetOptions & opt);

      /** Load a run directory (full or partial dataset, see class docs). */
      int loadDir(const char * dir);
      int loadDir(const char * dir, const DatasetOptions & opt);

      /** Load a single combined file; only information within this file is
       * read, so opt.partial_skip_incomplete is forced to true. */
      int loadCombinedFile(const char * file);
      int loadCombinedFile(const char * file, const DatasetOptions & opt);

      /** Replace the options; takes effect on the next load. If
       * opt.base_data_dir is empty it is resolved from $RNO_G_ROOT_DATA /
       * $RNO_G_DATA (or "."). */
      void setOpt(const DatasetOptions & opt = DatasetOptions());
      const DatasetOptions & getOpt() const { return opt; }
      void setVerbose(bool v) { opt.verbose = v; }
      void setCalibration(const VoltageCalibration * calib);
      const VoltageCalibration * getCalibration() const { return opt.calib; }

      /* ---------- entry selection ---------- */

      bool setEntry(int entry); //returns true if in range
      int currentEntry() const { return current_entry; }
      int N() const;

      /* ---------- per-entry data (lazily loaded, may return nullptr) ---------- */

      mattak::Header * header(bool force_reload = false);
      mattak::Waveforms * raw(bool force_reload = false);

      // is the raw data available for currentEntry? (mostly used by PyROOT backend)
      bool rawAvailable(bool force_reload = false)
      {
        return  wf.tree &&
         ( full_dataset || opt.partial_skip_incomplete ||
           wf.tree->GetEntryNumberWithIndex(header(force_reload)->event_number) >=0);
      }

      // these methods are useful if you want to read waveform metadata without reading the waveforms
      // if you are reading the waveforms, they are less efficient than getting what you want from raw
      float radiantSampleRate(bool force_reload = false);
      const float * radiantReadoutDelays(bool force_reload = false);  //size is mattak::k::num_radiant_channels, returning a float* since cppyyy doesn't seem to be able to deal with std::array properly

      mattak::CalibratedWaveforms * calibrated(bool force_reload = false); //will be nullptr if no calibration is passed
      mattak::DAQStatus * status(bool force_reload = false);
      mattak::RunInfo * info() const { return runinfo.ptr; }
      mattak::Pedestals * peds(bool force_reload = false, int entry = 0);

      TTree * daqStatusTree() { return ds.tree; }
      TTree * headTree() { return hd.tree; }
      TTree * wfTree() { return wf.tree; }

      bool isFullDataset() const { return full_dataset; }

      /* ---------- deprecated interface (kept for backwards compatibility) ---------- */

      [[deprecated("use Dataset(station, run, DatasetOptions) instead")]]
      Dataset (int station, int run, const VoltageCalibration * calib, const char * base_data_dir = nullptr, bool partial_skip_incomplete = true, bool verbose = false);

      [[deprecated("use Dataset(DatasetOptions) with base_data_dir set instead")]]
      Dataset (const char * data_dir);

      [[deprecated("use loadRun(station, run, DatasetOptions) instead")]]
      int loadRun(int station, int run, bool partial_skip_incomplete);

      [[deprecated("use loadDir(dir, DatasetOptions) instead")]]
      int loadDir(const char * dir, bool partial_skip_incomplete);

      [[deprecated("set DatasetOptions::base_data_dir instead")]]
      void setDataDir(const char * dir);

      /* ---------- implementation details (public for the loader helpers and
       * ROOT dictionary; do not use directly) ---------- */

      template <typename D>
      struct field
      {
        D * ptr = nullptr;
        int loaded_entry = -1;
        bool missing_entry = false;
      };

      template <typename D>
      struct tree_field : public field<D>
      {
        TFile * file = nullptr;
        TTree * tree = nullptr;
        TBranch * branch = nullptr;
      };

      template <typename D>
      struct file_field : public field<D>
      {
        TFile * file = nullptr;
      };

      //0 terminated arrays
      static const char ** getWaveformTreeNames();
      static const char ** getHeaderTreeNames();
      static const char ** getDAQStatusTreeNames();
      static const char ** getPedestalTreeNames();
    private:
      tree_field<Waveforms> wf;
      tree_field<Waveforms> wf_meta;
      tree_field<Header> hd;
      tree_field<DAQStatus> ds;
      tree_field<Pedestals> pd;
      file_field<RunInfo> runinfo;

      void setupRadiantMeta();

      field<CalibratedWaveforms> calib_wf;

      void unload();
      int current_entry = 0;


      bool full_dataset = false;
      DatasetOptions opt;

  };
}

#endif
