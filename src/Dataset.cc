#include "mattak/Dataset.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TError.h"



/** Delete the cached object and reset a field to its empty state
 * (as if nothing had ever been loaded into it). */
template <typename D>
static void clear(mattak::Dataset::field<D> * field)
{
  if (field->ptr)
  {
    delete field->ptr;
  }
  field->ptr = nullptr;
  field->loaded_entry = -1;
  field->missing_entry = false;
}


/** Close the owning TFile (which invalidates tree/branch) and clear the field. */
template <typename D>
static void clear(mattak::Dataset::tree_field<D> * field)
{
  if (field->file)
  {
    delete field->file;
  }
  field->file = nullptr;
  field->tree = nullptr;
  clear((mattak::Dataset::field<D>*) field);
}

/** Close the owning TFile and clear the field. */
template <typename D>
static void clear(mattak::Dataset::file_field<D> * field)
{
  if (field->file)
  {
    delete field->file;
  }
  field->file = nullptr;

  clear((mattak::Dataset::field<D>*) field);


}


#ifdef __WIN32
#define BITBUCKET "NUL"
#else
#define BITBUCKET "/dev/null"
#endif

/** Open a second, independent handle onto the already-opened waveform file/tree
 * (wf), reading only the radiant_sampling_rate and digitizer_readout_delay_ns
 * branches into wf_meta. This lets radiantSampleRate()/radiantReadoutDelays()
 * be read cheaply, without pulling in the (much larger) waveform samples that
 * a shared handle on wf itself would require. Called after every successful
 * waveform load; failures here are non-fatal (just leaves wf_meta cleared,
 * the getters then fall back to the runinfo/default sample rate). */
void mattak::Dataset::setupRadiantMeta()
{

  // HACK: separately make the sample rate from the waveform file
  clear(&wf_meta);

  //copy the file/tree so we don't double count branches
  wf_meta.file = TFile::Open(wf.file->GetName());
  if (!wf_meta.file)
  {
    ::Warning("mattak::Dataset::setupRadiantMeta", "Could not reopen %s", wf.file->GetName());
    return;
  }
  wf_meta.tree = (TTree*) wf_meta.file->Get(wf.tree->GetName());
  if (!wf_meta.tree)
  {
    ::Warning("mattak::Dataset::setupRadiantMeta", "Could not find tree %s in %s", wf.tree->GetName(), wf.file->GetName());
    clear(&wf_meta);
    return;
  }
  wf_meta.branch = nullptr; // we don't use this directly for this, since I'm not sure how it interacts with SetBranchStatus
  wf_meta.ptr = new mattak::Waveforms;
  wf_meta.tree->SetBranchAddress(wf.branch->GetName(), &wf_meta.ptr);
  wf_meta.tree->SetBranchStatus("*",0);
  UInt_t found = 0;
  wf_meta.tree->SetBranchStatus("*radiant_sampling_rate",1, &found);
  wf_meta.tree->SetBranchStatus("*digitizer_readout_delay_ns*",1, &found);
}

/** Silently check if file exists, supporting all protocols ROOT does */
static TFile * silentlyTryToOpen(const char * uri, const char * opt = "" )
{

  RedirectHandle_t rh;
  gSystem->RedirectOutput(BITBUCKET,"a",&rh);
  TFile * f= TFile::Open(uri,opt);
  gSystem->RedirectOutput(0,"",&rh);

  return f;
}

/** Open filename and attach the first matching tree/branch pair to the field.
 * Returns 0 on success. Returns -1 on failure (file missing, or no listed tree
 * has the expected branch), in which case the field is left fully cleared. */
template <typename D>
static int setup(mattak::Dataset::tree_field<D> * field, const char * filename, const char ** tree_names, const char ** branch_names = 0, bool verbose = false)
{
  clear(field);
  if (verbose) ::Info("mattak::Dataset::setup", "Trying to open %s", filename);
  field->file = !verbose ? silentlyTryToOpen(filename, "READ") : TFile::Open(filename, "READ");
  if (!field->file) return -1;

  // a "combined" tree may hold the branch under any of the candidate names
  TTree * combined_tree = (TTree*) field->file->Get("combined");

  int itry = 0;
  while(tree_names[itry])
  {
    if (verbose) ::Info("mattak::Dataset::setup", "Trying tree %s", tree_names[itry]);
    field->tree = (TTree*) field->file->Get(tree_names[itry]);
    if (!field->tree)
    {
      if (verbose) ::Info("mattak::Dataset::setup", "Trying tree combined");
      field->tree = combined_tree;
    }

    if (!field->tree)
    {
      itry++;
      continue;
    }

    const char * branch_name = branch_names ? branch_names[itry] : tree_names[itry];
    if (verbose) ::Info("mattak::Dataset::setup", "Trying branch %s", branch_name);
    if (!field->tree->GetBranch(branch_name))
    {
      itry++;
      continue;
    }

    //avoid annoying message for empty string branch
    RedirectHandle_t rh;
    if (!branch_name[0]) gSystem->RedirectOutput(BITBUCKET,"a",&rh);

    field->branch = field->tree->GetBranch(branch_name);
    field->branch->SetAddress(&field->ptr);
    if (verbose) ::Info("mattak::Dataset::setup", "Found!");

    if (!branch_name[0]) gSystem->RedirectOutput(0,"a",&rh);

    gROOT->cd();
    return 0;
  }
  if (verbose) ::Info("mattak::Dataset::setup", "Could not find a valid tree/branch pair in %s", filename);
  clear(field); // don't leave a half-initialized field (open file, tree without branch/address)
  return -1;
}

/** Open filename and read the object obj_name from it into the field.
 * Returns 0 on success. Returns -1 on failure (file missing or object not
 * found), in which case the field is left fully cleared.
 * Same convention as the tree_field overload above. */
template <typename D>
static int setup(mattak::Dataset::file_field<D> * field, const char * filename, const char * obj_name, bool verbose = false)
{
  clear(field);
  field->file = !verbose ? silentlyTryToOpen(filename,"READ") : TFile::Open(filename,"READ");
  if (!field->file) return -1;
  field->ptr = (D*) field->file->Get(obj_name);
  gROOT->cd();
  if (!field->ptr)
  {
    clear(field);
    return -1;
  }
  return 0;
}


/** Close every open file and reset every field to empty. Called at the start
 * of every load and from the destructor, so a Dataset is always either fully
 * loaded or fully empty, never a mix of the two. */
void mattak::Dataset::unload()
{
  clear(&wf);
  clear(&wf_meta);
  clear(&hd);
  clear(&ds);
  clear(&pd);
  clear(&runinfo);
  clear(&calib_wf);
}


mattak::Dataset::Dataset(int station, int run, const DatasetOptions & opt)
{
  loadRun(station, run, opt);
}

mattak::Dataset::Dataset(const DatasetOptions & opt)
{
  setOpt(opt);
}


/** Resolve the default base data directory from the environment */
static std::string defaultDataDir()
{
  const char * env = getenv("RNO_G_ROOT_DATA");
  if (!env) env = getenv("RNO_G_DATA");
  return env ? env : ".";
}

void mattak::Dataset::setOpt(const DatasetOptions & opt)
{
  this->opt = opt;
  if (this->opt.base_data_dir == "")
  {
    this->opt.base_data_dir = defaultDataDir();
  }
}

/* deprecated, forwards to the DatasetOptions interface */
mattak::Dataset::Dataset(int station, int run, const VoltageCalibration * calib, const char * data_dir, bool partial_skip, bool v)
{
  DatasetOptions o;
  o.verbose = v;
  o.calib = calib;
  o.partial_skip_incomplete = partial_skip;
  if (data_dir)
    o.base_data_dir = data_dir;

  loadRun(station, run, o);
}

/* deprecated, forwards to the DatasetOptions interface */
mattak::Dataset::Dataset(const char* data_dir)
{
  opt.base_data_dir = data_dir ? data_dir : defaultDataDir();
}

/* deprecated, forwards to the DatasetOptions interface */
void mattak::Dataset::setDataDir(const char * dir)
{
  opt.base_data_dir = dir ? dir : defaultDataDir();
}

void mattak::Dataset::setCalibration(const VoltageCalibration * c)
{
  opt.calib = c;
}

const char * waveform_tree_names[] = {"waveforms","wfs","wf","waveform",0};
const char * header_tree_names[] = {"hdr","header","hd","hds","headers",0};
const char * daqstatus_tree_names[] = {"daqstatus","ds","status",0};
const char * pedestal_tree_names[] = {"pedestals","pedestal","ped","peds","",0};

const char ** mattak::Dataset::getWaveformTreeNames()
{
  return waveform_tree_names;
}

const char ** mattak::Dataset::getHeaderTreeNames()
{
  return header_tree_names;
}

const char ** mattak::Dataset::getDAQStatusTreeNames()
{
  return daqstatus_tree_names;
}

const char ** mattak::Dataset::getPedestalTreeNames()
{
  return pedestal_tree_names;
}




int mattak::Dataset::loadRun(int station, int run, bool partial_skip)
{
  opt.partial_skip_incomplete = partial_skip;
  return loadRun(station, run);
}

int mattak::Dataset::loadRun(int station, int run, const DatasetOptions & opt)
{
  setOpt(opt);
  return loadRun(station, run);
}

/** Build the canonical run directory path (base_data_dir/station<S>/run<R>)
 * and delegate to loadDir. This is the only place that path convention is
 * encoded, so callers never need to construct it themselves. */
int mattak::Dataset::loadRun(int station, int run)
{
  TString dir;
  dir.Form("%s/station%d/run%d", opt.base_data_dir.c_str(), station, run);
  return loadDir(dir.Data());
}



int mattak::Dataset::loadCombinedFile(const char * f, const DatasetOptions & opt)
{
  setOpt(opt);
  return loadCombinedFile(f);

}


int mattak::Dataset::loadDir(const char * dir, const DatasetOptions & opt)
{
  setOpt(opt);
  return loadDir(dir);

}

int mattak::Dataset::loadDir(const char * dir, bool partial_skip)
{
  opt.partial_skip_incomplete = partial_skip;
  return loadDir(dir);
}

/** Load every field from a single combined-style file: waveforms and headers
 * are mandatory (failure unloads and returns -1), daqstatus/pedestals/runinfo
 * are optional (their getters just return nullptr if absent). Because only
 * this one file is read, there is no way to distinguish complete from
 * incomplete events, so partial_skip_incomplete is always forced to true. */
int mattak::Dataset::loadCombinedFile(const char * f)
{
  if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "loadCombinedFile(%s) called", f);

  unload();
  current_entry = 0;
  full_dataset = false;

  // only information within this file is read, so we can never iterate over incomplete events
  if (!opt.partial_skip_incomplete)
  {
    ::Warning("mattak::Dataset::loadCombinedFile", "partial_skip_incomplete=false is incompatible with loadCombinedFile, forcing to true");
    opt.partial_skip_incomplete  = true;
  }

  if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Opening %s", f);
  if (setup(&wf, f, waveform_tree_names, nullptr, opt.verbose) != 0
   || setup(&hd, f, header_tree_names, nullptr, opt.verbose) != 0)
  {
    ::Error("mattak::Dataset::loadCombinedFile", "Could not load waveforms and headers from %s", f);
    unload();
    return -1;
  }

  setupRadiantMeta();

  if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Found waveforms and headers in %s", f);

  // daqstatus, pedestals and run info are optional: their getters return nullptr if missing
  if (setup(&ds, f, daqstatus_tree_names, nullptr, opt.verbose) == 0)
  {
    if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Found daqstatus in %s", f);
  }
  else
  {
    ::Warning("mattak::Dataset::loadCombinedFile", "Could not load daqstatus from %s (this is ok if you don't use them)", f);
  }

  if (setup(&pd, f, pedestal_tree_names, nullptr, opt.verbose) == 0)
  {
    if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Found pedestals in %s", f);
  }

  // the run info may be stored under either name, try both
  bool found_runinfo = setup(&runinfo, f, "info", opt.verbose) == 0
                    || setup(&runinfo, f, "runinfo", opt.verbose) == 0;
  if (found_runinfo)
  {
    if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Found runinfo in %s", f);
  }
  else
  {
    ::Warning("mattak::Dataset::loadCombinedFile", "Could not load run info from %s", f);
  }

  return 0;
}

/** Load every field from a run directory. This is the workhorse behind
 * loadRun/loadCombinedFile/the public loadDir overloads, and does three
 * things in order:
 *   1. Figure out where the waveforms live: opt.file_preference if set,
 *      else waveforms.root (a "full dataset"), else combined.root (a
 *      "partial dataset" containing only the events with waveforms).
 *   2. Load headers and daqstatus. For a full dataset (or when the caller
 *      wants incomplete events too, partial_skip_incomplete=false) these
 *      come from their own standalone files; otherwise they are read out of
 *      the same combined-style file as the waveforms, which limits
 *      iteration to the complete events it contains.
 *   3. Load the optional pedestal/runinfo fields, whose getters simply
 *      return nullptr if unavailable.
 * Waveforms/headers/daqstatus are mandatory: any failure to load them
 * unloads the dataset and returns -1. */
int mattak::Dataset::loadDir(const char * dir)
{

  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "loadDir(%s, skip_incomplete=%d) called", dir, opt.partial_skip_incomplete);

  unload();
  current_entry = 0;
  full_dataset = false;

  /* Waveforms. Preference order: opt.file_preference, then a full dataset
   * (waveforms.root), then a partial one (combined.root). partial_file holds
   * the name (without .root) of the combined-style file if we use one. */
  const char * partial_file = nullptr;

  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "Load waveforms ...");

  if (opt.file_preference != "")
  {
    if (setup(&wf, Form("%s/%s.root", dir, opt.file_preference.c_str()), waveform_tree_names) == 0)
    {
      partial_file = opt.file_preference.c_str();
    }
    else
    {
      ::Warning("mattak::Dataset::loadDir", "Could not find preferred %s.root in %s. Reverting to default behavior", opt.file_preference.c_str(), dir);
    }
  }

  if (!partial_file)
  {
    if (setup(&wf, Form("%s/waveforms.root", dir), waveform_tree_names) == 0)
    {
      if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... full dataset found");
      full_dataset = true;
    }
    else if (setup(&wf, Form("%s/combined.root", dir), waveform_tree_names, nullptr, opt.verbose) == 0)
    {
      if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... full dataset not found, using combined.root");
      partial_file = "combined";
    }
    else
    {
      ::Error("mattak::Dataset::loadDir", "Failed to find waveforms.root or combined.root in %s", dir);
      return -1;
    }
  }

  setupRadiantMeta();

  /* Headers. The standalone headers.root is needed for a full dataset and
   * when the user wants to iterate over incomplete events too; otherwise the
   * headers come from the combined-style file. */
  bool want_full_headers = full_dataset || !opt.partial_skip_incomplete;
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "About to load headers");
  if (setup(&hd, Form("%s/%s.root", dir, want_full_headers ? "headers" : partial_file), header_tree_names, nullptr, opt.verbose) != 0)
  {
    // next to a combined-style file, headers.root may legitimately be absent:
    // fall back to the headers inside it, which limits us to complete events
    bool fell_back = false;
    if (!full_dataset && !opt.partial_skip_incomplete)
    {
      if (setup(&hd, Form("%s/%s.root", dir, partial_file), header_tree_names, nullptr, opt.verbose) == 0)
      {
        ::Warning("mattak::Dataset::loadDir", "Could not find headers.root in %s; using %s.root instead (only complete events)", dir, partial_file);
        opt.partial_skip_incomplete = true;
        fell_back = true;
      }
    }

    if (!fell_back)
    {
      ::Error("mattak::Dataset::loadDir", "Failed to load headers from %s", dir);
      unload();
      return -1;
    }
  }
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... success");

  if (!full_dataset && !opt.partial_skip_incomplete)
  {
    // index the partial waveforms by event number so raw() can find them from the header entry
    wf.tree->BuildIndex("event_number");
  }

  /* DAQ status: same source selection as the headers */
  bool want_full_daqstatus = full_dataset || !opt.partial_skip_incomplete;
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "About to load daqstatus");
  if (setup(&ds, Form("%s/%s.root", dir, want_full_daqstatus ? "daqstatus" : partial_file), daqstatus_tree_names, nullptr, opt.verbose) != 0)
  {
    ::Error("mattak::Dataset::loadDir", "Failed to load %s.root in %s", want_full_daqstatus ? "daqstatus" : partial_file, dir);
    unload();
    return -1;
  }
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... success");

  if (full_dataset)
  {
    // index the daqstatus by readout time so status() can match it to the current event
    ds.tree->BuildIndex("int(readout_time_radiant)", "1e9*(readout_time_radiant-int(readout_time_radiant))");
  }

  /* Pedestals and run info are optional: their getters return nullptr if missing */
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "About to load pedestal");
  if (setup(&pd, Form("%s/pedestal.root", dir), pedestal_tree_names, nullptr, opt.verbose) == 0)
  {
    if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... success");
  }
  else
  {
    ::Warning("mattak::Dataset::loadDir", "Failed to find pedestal.root in %s (this is usually ok if you don't need them)", dir);
  }

  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "About to load runinfo");
  if (full_dataset)
  {
    if (setup(&runinfo, Form("%s/runinfo.root", dir), "info", opt.verbose) == 0)
    {
      if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... success");
    }
    else
    {
      ::Warning("mattak::Dataset::loadDir", "Failed to read runinfo from %s/runinfo.root", dir);
    }
  }
  else
  {
    // for a partial/combined run the runinfo is stored as the "info"/"runinfo"
    // object inside the combined-style file itself (same as loadCombinedFile)
    TString combined = Form("%s/%s.root", dir, partial_file);
    if (setup(&runinfo, combined.Data(), "info") == 0 || setup(&runinfo, combined.Data(), "runinfo") == 0)
    {
      if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... success");
    }
    else
    {
      ::Warning("mattak::Dataset::loadDir", "Failed to read runinfo from %s/%s.root", dir, partial_file);
    }
  }

  return 0;
}


bool mattak::Dataset::setEntry(int entry)
{
  if (entry >= 0 && entry < N())
  {
    current_entry = entry;  //otherwise we're completely lazy!
    return true;
  }
  return false;
}

int mattak::Dataset::N() const
{
  if (hd.tree==nullptr) return -1;
  return hd.tree->GetEntries();
}

/** Return the header for currentEntry, reading it from the tree only if it
 * isn't already cached (or force is set). Headers are always present for
 * every indexed entry, so this never needs the missing-entry handling that
 * raw()/status() require. */
mattak::Header* mattak::Dataset::header(bool force)
{
  if (force || hd.loaded_entry != current_entry)
  {
    if (hd.tree == nullptr) return nullptr;
    hd.branch->GetEntry(current_entry);
    hd.loaded_entry = current_entry;
  }

  return hd.ptr;
}

/** Load the entry of `field` matching the dataset's currentEntry, caching by
 * loaded_entry the same way header() does. The difference from a plain
 * GetEntry(current) is that field may come from a partial dataset that only
 * has some events indexed by event_number: if the dataset is full (or the
 * caller asked to skip incomplete events, so field is guaranteed aligned
 * with the headers), read straight by position; otherwise look up the entry
 * via the tree's (or index_field's) event_number index and mark
 * field->missing_entry if this event isn't present at all. */
template <typename T>
static void findIncompleteEntry(mattak::Dataset::tree_field<T> * field, mattak::Dataset * d, bool force = false, mattak::Dataset::tree_field<mattak::Waveforms> * index_field = nullptr)
{
  if (force || field->loaded_entry != d->currentEntry())
  {
    if (d->isFullDataset() || d->getOpt().partial_skip_incomplete)
    {
      auto current = d->currentEntry();
      if (field->branch) field->branch->GetEntry(current);
      else field->tree->GetEntry(current);
      field->loaded_entry = current;
    }
    else
    {
      int entry= -1;
      if (index_field)
      {
        entry = index_field->tree ? index_field->tree->GetEntryNumberWithIndex(d->header(force)->event_number) : -1;
      }
      else
      {
        entry = field->tree ? field->tree->GetEntryNumberWithIndex(d->header(force)->event_number) : -1;
      }

      field->loaded_entry = d->currentEntry();

      if (entry < 0)
      {
        field->missing_entry = true;
      }
      else
      {
        field->missing_entry = false;
        if (field->branch) field->branch->GetEntry(entry);
        else field->tree->GetEntry(entry);
      }
    }
  }
}


/** RADIANT sampling rate (MHz) for currentEntry, read cheaply from wf_meta
 * (see setupRadiantMeta) without touching the full waveform. Falls back to
 * the runinfo rate, then to a hardcoded default, if wf_meta isn't set up or
 * this entry has no waveform metadata. */
float mattak::Dataset::radiantSampleRate(bool force)
{
  if (wf_meta.ptr == nullptr) {
    return (info() && info()->radiant_sample_rate) ? info()->radiant_sample_rate : k::default_radiant_sample_rate;
  }
  findIncompleteEntry(&wf_meta, this, force, &wf);
  if (wf_meta.missing_entry)
  {
    return (info() && info()->radiant_sample_rate) ? info()->radiant_sample_rate : k::default_radiant_sample_rate;
  }
  return wf_meta.ptr->radiant_sampling_rate;
}

static float zeros[mattak::k::num_radiant_channels];

/** Per-channel RADIANT digitizer readout delay (ns) for currentEntry, read
 * from wf_meta the same way as radiantSampleRate. Returns an all-zero array
 * if unavailable (there is no runinfo fallback for this one). */
const float *  mattak::Dataset::radiantReadoutDelays(bool force)
{
  if (wf_meta.ptr == nullptr)
  {
    return (const float*) zeros;
  }
  findIncompleteEntry(&wf_meta, this, force, &wf);
  return (const float*)  (
    wf_meta.missing_entry ? zeros : wf_meta.ptr->digitizer_readout_delay_ns
  );
}


/** Raw (uncalibrated) waveforms for currentEntry, or nullptr if this event
 * has none (partial dataset with an event not present in the waveform
 * file) or the loaded waveform entry doesn't line up with the current
 * header's event number (a sanity check against a misaligned index). */
mattak::Waveforms* mattak::Dataset::raw(bool force)
{
  if (wf.tree == nullptr)
    return nullptr;

  findIncompleteEntry(&wf, this, force);

  if (wf.missing_entry)
    return nullptr;

  if (wf.ptr->event_number != header(force)->event_number)
    return nullptr;

  return wf.ptr;
}


/** DAQ status for currentEntry, or nullptr if no daqstatus tree was loaded.
 * For a full dataset, daqstatus is sampled independently in time (not one
 * entry per event), so the closest entry by readout time is looked up via
 * the tree's time index; for a combined-style file, daqstatus is stored one
 * entry per event and can be read straight by position. */
mattak::DAQStatus * mattak::Dataset::status(bool force)
{
  if (!ds.tree) return nullptr;
  if (force || ds.loaded_entry != current_entry)
  {
    if (full_dataset)
    {
      // the standalone daqstatus is sampled in time, look up the entry closest to the event
      double readout_time = header(force)->readout_time;
      int ds_entry = ds.tree->GetEntryNumberWithBestIndex(readout_time, 1e9 * (readout_time - int(readout_time)));
      if (ds_entry < 0) ds_entry = 0;  // this should only happen if it's the first one?
      ds.branch->GetEntry(ds_entry);
      ds.missing_entry = false;

    }
    else
    {
      // combined-style files store one daqstatus entry per event, aligned with the headers.
      // NB: this assumes the daqstatus was NOT loaded from a standalone daqstatus.root
      // (currently possible for partial datasets with partial_skip_incomplete=false)
      ds.branch->GetEntry(current_entry);
      ds.missing_entry = false;
    }

    ds.loaded_entry = current_entry;
  }

  return ds.missing_entry ? nullptr: ds.ptr;
}

/** Voltage-calibrated waveforms for currentEntry, computed on demand from
 * raw() and header() using opt.calib. Returns nullptr if no calibration was
 * set, or if raw()/header() are unavailable for this entry. The result is
 * cached (and reused via placement-new into the same buffer) until
 * currentEntry changes or force is set, so repeated calls for the same
 * entry don't recompute the calibration. */
mattak::CalibratedWaveforms * mattak::Dataset::calibrated(bool force)
{
  //no calibration? we can't calibrate.
  if (!opt.calib) return nullptr;

  if (force || calib_wf.loaded_entry != current_entry)
  {
    mattak::Waveforms * raw_wf = raw(force);
    mattak::Header * head = header(force);

    // if there is no raw waveform, we can't do this
    if (!raw_wf || !head)
    {
      calib_wf.missing_entry = true;
    }
    else
    {
      calib_wf.missing_entry = false;

      if (!calib_wf.ptr)
      {
        calib_wf.ptr = new CalibratedWaveforms(*raw_wf, *head, *opt.calib);
      }
      else
      {
        new (calib_wf.ptr) CalibratedWaveforms(*raw_wf, *head, *opt.calib);
      }
    }

    calib_wf.loaded_entry = current_entry;
  }

  return calib_wf.missing_entry ? nullptr : calib_wf.ptr;
}


/** Pedestals for a given pedestal-tree entry (not a dataset event: pedestal
 * runs are recorded separately and indexed independently, hence the explicit
 * `entry` argument rather than using currentEntry). Returns nullptr if no
 * pedestal tree was loaded or entry is out of range. */
mattak::Pedestals * mattak::Dataset::peds(bool force, int entry)
{
  if (! pd.tree) return nullptr;
  if (entry < 0 || entry >= pd.tree->GetEntries()) return nullptr;

  if (force || entry != pd.loaded_entry)
  {
    pd.branch->GetEntry(entry);
    pd.loaded_entry = entry;
  }
  return pd.ptr;
}


//HACK HACK HACK
//davix  seems to choke here for some reason, at least for me.
//check for an environmental variable called "MATTAK_SUPPRESS_DAVIX"
__attribute__((constructor))
static void maybe_kill_davix()
{
  char * suppress = getenv("MATTAK_SUPPRESS_DAVIX");

  if (!suppress || !strcmp(suppress,"0")) return;

 // tell ROOT to load all of its plugin handlers, otherwise the first time you open a file this will happen again and override what you are about to do after this
 gPluginMgr->LoadHandlersFromPluginDirs();

  // Override the plugin handler for web files to use the legacy TWebFile instead of the newer davix which seems to be buggy
 gPluginMgr->AddHandler("TFile", "^http[s]?:", "TWebFile","Net", "TWebFile(const char*,Option_t*)");

}
