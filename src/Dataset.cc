#include "mattak/Dataset.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TError.h"



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

void mattak::Dataset::setupRadiantMeta()
{

  // HACK: separately make the sample rate from the waveform file
  clear(&wf_meta);

  if (!wf.file || !wf.tree) return;  // headers-only dataset: no waveform metadata to read

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
  // Returns 0 on success (found a valid tree + branch pair).
  // Returns -1 on failure (file couldn't open, or no matching tree/branch).
  clear(field);
  if (verbose) ::Info("mattak::Dataset::setup", "Trying to open %s", filename);
  field->file = !verbose ? silentlyTryToOpen(filename,"READ") : TFile::Open(filename,"READ");
  if (!field->file) return -1;

  int itry = 0;
  while(tree_names[itry])
  {
    field->tree = (TTree*) field->file->Get(tree_names[itry]);
    if (verbose) ::Info("mattak::Dataset::setup", "Trying tree %s", tree_names[itry]);
    if (!field->tree)
    {
      field->tree = (TTree*) field->file->Get("combined");
      if (verbose) ::Info("mattak::Dataset::setup", "Trying tree combined");
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
  // Returns 0 on success (object found).
  // Returns -1 on failure (file couldn't open, or object not there).
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


void mattak::Dataset::setOpt(const DatasetOptions & opt)
{

  this->opt = opt;
  if (opt.base_data_dir == "")
  {
    setDataDir(nullptr);
  }

}

mattak::Dataset::Dataset(int station, int run, const VoltageCalibration * calib, const char * data_dir, bool partial_skip, bool v)
{
  setVerbose(v);  // should be first
  setDataDir(data_dir);
  setCalibration(calib);
  loadRun(station, run, partial_skip);
}

mattak::Dataset::Dataset(const char* data_dir)
{
  setDataDir(data_dir);
}


void mattak::Dataset::setDataDir(const char * dir)
{
  if (dir)
  {
    opt.base_data_dir = dir;
  }
  else
  {
    const char * env = getenv("RNO_G_ROOT_DATA");
    if (!env) env = getenv("RNO_G_DATA");
    opt.base_data_dir = env ? env : ".";
  }
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
  return loadRun(station,run);
}

int mattak::Dataset::loadRun(int station, int run, const DatasetOptions & opt)
{
  setOpt(opt);
  return loadRun(station,run);
}

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

int mattak::Dataset::loadCombinedFile(const char * f)
{
  if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "loadCombinedFile(%s) called", f);

  unload();
  current_entry = 0;
  full_dataset = false;

  if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Opening %s", f);
  if (setup(&hd, f, header_tree_names, 0, opt.verbose))
  {
    ::Error("mattak::Dataset::loadCombinedFile", "Could not load headers from %s", f);
    return -1;
  }

  if (setup(&wf, f, waveform_tree_names, 0, opt.verbose))
  {
    // headers-only file: same deal as in loadDir, index all events and let the waveform
    // accessors return nullptr.
    ::Warning("mattak::Dataset::loadCombinedFile",
              "Could not load waveforms from %s, loading headers only "
              "(forcing partial_skip_incomplete=false)", f);
    opt.partial_skip_incomplete = false;
  }
  else if (!opt.partial_skip_incomplete)
  {
    ::Warning("mattak::Dataset::loadCombinedFile", "partial_skip_incomplete=false is incompatible with loadCombinedFile, forcing to true");
    opt.partial_skip_incomplete  = true;
  }

  setupRadiantMeta();

  if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Found headers in %s", f);

  // Try some optionalish things
  if (setup(&ds, f, daqstatus_tree_names, 0, opt.verbose))
  {
    ::Warning("mattak::Dataset::loadCombinedFile", "Could not load daqstatus from %s (this is ok if you don't use them)", f);
  }
  else
  {
    if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Found daqstatus in %s", f);
  }

  // without waveforms the daqstatus is looked up by readout time (see status()), so it needs an index
  if (!wf.tree && ds.tree)
  {
    ds.tree->BuildIndex("int(readout_time_radiant)", "1e9*(readout_time_radiant-int(readout_time_radiant))");
  }

  // we probably don't have pedetals, but we could try I guess?
  if (!setup(&pd, f, pedestal_tree_names, 0, opt.verbose))
  {
    if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Found pedestals in %s", f);
  }

  if ( !setup(&runinfo, f, "info", opt.verbose) || !setup(&runinfo, f, "runinfo", opt.verbose) )
  {
    if (opt.verbose) ::Info("mattak::Dataset::loadCombinedFile", "Found runinfo in %s", f);
  }
  else {
    ::Warning("mattak::Dataset::loadCombinedFile", "Could not load run info from %s", f);
  }

  return 0;
}

int mattak::Dataset::loadDir(const char * dir)
{

  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "loadDir(%s, skip_incomplete=%d) called", dir, opt.partial_skip_incomplete);

  //first clear all
  unload();
  current_entry = 0;

  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "Load waveforms ...");

  const char * partial_file = NULL;
  if (opt.file_preference != "" && !setup(&wf, Form("%s/%s.root",dir,opt.file_preference.c_str()), waveform_tree_names))
  {
    full_dataset = false;
    partial_file = opt.file_preference.c_str();
  }
  else
  {
    if (opt.file_preference != "")
    {
      ::Warning("mattak::Dataset::loadDir", "Could not find preferred %s.root in %s. Reverting to default behavior", opt.file_preference.c_str(), dir);
    }

    //we need to figure out if this is a full run or partial run, so check for existence of waveforms.root
    if (setup(&wf, Form("%s/waveforms.root", dir), waveform_tree_names, 0))
    {
      //no waveforms file!
      full_dataset = false;
      if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... full dataset not found");

      //let's load from combined file instead
      if (setup(&wf, Form("%s/combined.root", dir), waveform_tree_names, nullptr, opt.verbose))
      {
        // No waveforms anywhere. Fall back to a headers-only dataset: headers.root (loaded
        // below) is then the only hard requirement and all waveform accessors return nullptr.
        // That is precisely partial_skip_incomplete=false (index all events, waveforms may be
        // missing), which also makes the file name selections below pick the standalone files.
        ::Warning("mattak::Dataset::loadDir",
                  "Failed to find waveforms.root or combined.root in %s, loading headers only "
                  "(forcing partial_skip_incomplete=false)", dir);
        opt.partial_skip_incomplete = false;
      }
      else
      {
        partial_file = "combined";
      }
    }
    else
    {
      if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... full dataset found");
      full_dataset = true;
    }
  }

  setupRadiantMeta();

  //now load the header files
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "About to load headers");
  const char * hd_file = (full_dataset || !opt.partial_skip_incomplete) ? "headers" : partial_file;
  if (setup(&hd, Form("%s/%s.root", dir, hd_file), header_tree_names, nullptr, opt.verbose))
  {
    ::Error("mattak::Dataset::loadDir", "Failed to load %s.root in %s", hd_file, dir);
    return -1;
  }
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... success");

  if (wf.tree && !full_dataset && !opt.partial_skip_incomplete)
  {
    //set up an index on event number the events
    wf.tree->BuildIndex("event_number");
  }

  //and the status files
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "About to load daqstatus");
  const char * ds_file = (full_dataset || !opt.partial_skip_incomplete) ? "daqstatus" : partial_file;
  if (setup(&ds, Form("%s/%s.root", dir, ds_file), daqstatus_tree_names, nullptr, opt.verbose))
  {
    ::Warning("mattak::Dataset::loadDir", "Failed to load %s.root in %s (this is ok if you don't use them)", ds_file, dir);
  }
  else if (opt.verbose)
  {
    ::Info("mattak::Dataset::loadDir", " ... success");
  }

  // In a headers-only run the daqstatus is a standalone, asynchronously written file just
  // like in a full run, so it must be looked up by readout time (see status()).
  if ((full_dataset || !wf.tree) && ds.tree)
  {
    ds.tree->BuildIndex("int(readout_time_radiant)", "1e9*(readout_time_radiant-int(readout_time_radiant))");
  }

  //and the pedestal files
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "About to load pedestal");
  if (setup(&pd, Form("%s/pedestal.root", dir), pedestal_tree_names, nullptr, opt.verbose))
  {
    ::Warning("mattak::Dataset::loadDir", "Failed to find pedestal.root in %s (this is usually ok if you don't need them)", dir);
  }
  else
  {
    if (opt.verbose) ::Info("mattak::Dataset::loadDir", " ... success");
  }

  //and try the runinfo file
  if (opt.verbose) ::Info("mattak::Dataset::loadDir", "About to load runinfo");
  if (full_dataset || !wf.tree)
  {
    if (setup(&runinfo, Form("%s/runinfo.root", dir), "info",  opt.verbose) == 0)
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
    // For a partial/combined run the runinfo is stored as the "info"/"runinfo"
    // object inside the combined file itself, so read it from there (same as
    // loadCombinedFile).
    const char * combined = Form("%s/%s.root", dir, partial_file);
    if (setup(&runinfo, combined, "info") == 0 || setup(&runinfo, combined, "runinfo") == 0)
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


float mattak::Dataset::radiantSampleRate(bool force)
{
  if (wf_meta.ptr == nullptr) {
    return (info() && info()->radiant_sample_rate) ? info()->radiant_sample_rate : 3200;
  }
  findIncompleteEntry(&wf_meta, this, force, &wf);
  if (wf_meta.missing_entry)
  {
    return (info() && info()->radiant_sample_rate) ? info()->radiant_sample_rate : 3200;
  }
  return wf_meta.ptr->radiant_sampling_rate;
}

static float zeros[mattak::k::num_radiant_channels];

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


mattak::DAQStatus * mattak::Dataset::status(bool force)
{
  if (!ds.tree) return nullptr;
  if (force || ds.loaded_entry != current_entry)
  {
    if (full_dataset || !wf.tree)
    {
      double readout_time = header(force)->readout_time;
      int ds_entry = ds.tree->GetEntryNumberWithBestIndex(readout_time, 1e9 * (readout_time - int(readout_time)));
      if (ds_entry < 0) ds_entry = 0;  // this should only happen if it's the first one?
      ds.branch->GetEntry(ds_entry);
      ds.missing_entry = false;

    }
    else
    {
      ds.branch->GetEntry(current_entry);
      ds.missing_entry = false;
    }

    ds.loaded_entry = current_entry;
  }

  return ds.missing_entry ? nullptr: ds.ptr;
}

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
