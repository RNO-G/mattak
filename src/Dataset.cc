#include "mattak/Dataset.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include <iostream>



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

  //copy the file/tree so we don't double count branches
  wf_meta.file = TFile::Open(wf.file->GetName());
  if (!wf_meta.file)
  {
    std::cerr << "setupRadiantMeta: could not reopen " << wf.file->GetName() << std::endl;
    return;
  }
  wf_meta.tree = (TTree*) wf_meta.file->Get(wf.tree->GetName());
  if (!wf_meta.tree)
  {
    std::cerr << "setupRadiantMeta: could not find tree " << wf.tree->GetName() << " in " << wf.file->GetName() << std::endl;
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
  if (verbose) std::cout << "Trying to open " << filename << std::endl;
  field->file = !verbose ? silentlyTryToOpen(filename, "READ") : TFile::Open(filename, "READ");
  if (!field->file) return -1;

  // a "combined" tree may hold the branch under any of the candidate names
  TTree * combined_tree = (TTree*) field->file->Get("combined");

  int itry = 0;
  while(tree_names[itry])
  {
    if (verbose) std::cout << "trying tree " << tree_names[itry] << std::endl;
    field->tree = (TTree*) field->file->Get(tree_names[itry]);
    if (!field->tree)
    {
      if (verbose) std::cout << "trying tree combined" << std::endl;
      field->tree = combined_tree;
    }

    if (!field->tree)
    {
      itry++;
      continue;
    }

    const char * branch_name = branch_names ? branch_names[itry] : tree_names[itry];
    if (verbose) std::cout << "trying branch " << branch_name << std::endl;
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
    if (verbose) std::cout << "Found!" << std::endl;

    if (!branch_name[0]) gSystem->RedirectOutput(0,"a",&rh);

    gROOT->cd();
    return 0;
  }
  if (verbose) std::cerr << "Could not find a valid tree/branch pair in " << filename << std::endl;
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
  if (opt.verbose) std::cout << "mattak::Dataset::loadCombinedFile (" << f  << ") called" << std::endl;

  unload();
  current_entry = 0;
  full_dataset = false;

  // only information within this file is read, so we can never iterate over incomplete events
  if (!opt.partial_skip_incomplete)
  {
    std::cerr << "partial_skip_incomplete is incompatible with loadCombinedFile, forcing it to true" << std::endl;
    opt.partial_skip_incomplete  = true;
  }

  if (opt.verbose) std::cout << "Opening " << f << std::endl;
  if (setup(&wf, f, waveform_tree_names, nullptr, opt.verbose) != 0
   || setup(&hd, f, header_tree_names, nullptr, opt.verbose) != 0)
  {
    std::cerr << "Could not load waveforms and headers from " << f << std::endl;
    unload();
    return -1;
  }

  setupRadiantMeta();

  if (opt.verbose) std::cout << "Found waveforms and headers in " << f << std::endl;

  // daqstatus, pedestals and run info are optional: their getters return nullptr if missing
  if (setup(&ds, f, daqstatus_tree_names, nullptr, opt.verbose) == 0)
  {
    if (opt.verbose) std::cout << "Found daqstatus in " << f << std::endl;
  }
  else
  {
    std::cerr << "Could not load daqstatus from " << f << " (this is ok if you don't use them)" << std::endl;
  }

  if (setup(&pd, f, pedestal_tree_names, nullptr, opt.verbose) == 0)
  {
    if (opt.verbose) std::cout << "Found pedestals in " << f << std::endl;
  }

  // the run info may be stored under either name, try both
  bool found_runinfo = setup(&runinfo, f, "info", opt.verbose) == 0
                    || setup(&runinfo, f, "runinfo", opt.verbose) == 0;
  if (found_runinfo)
  {
    if (opt.verbose) std::cout << "Found runinfo in " << f << std::endl;
  }
  else
  {
    std::cerr << "Could not load run info from " << f << std::endl;
  }

  return 0;
}

int mattak::Dataset::loadDir(const char * dir)
{

  if (opt.verbose) std::cout << "mattak::Dataset::loadDir (" << dir  << ", skip_incomplete=" << opt.partial_skip_incomplete << ") called" << std::endl;

  unload();
  current_entry = 0;
  full_dataset = false;

  /* Waveforms. Preference order: opt.file_preference, then a full dataset
   * (waveforms.root), then a partial one (combined.root). partial_file holds
   * the name (without .root) of the combined-style file if we use one. */
  const char * partial_file = nullptr;

  if (opt.verbose) std::cout << "Load waveforms ..." << std::endl;

  if (opt.file_preference != "")
  {
    if (setup(&wf, Form("%s/%s.root", dir, opt.file_preference.c_str()), waveform_tree_names) == 0)
    {
      partial_file = opt.file_preference.c_str();
    }
    else
    {
      std::cerr << "Warning, could not find preferred " << opt.file_preference << ".root in " << dir << ". Reverting to default behavior" << std::endl;
    }
  }

  if (!partial_file)
  {
    if (setup(&wf, Form("%s/waveforms.root", dir), waveform_tree_names) == 0)
    {
      if (opt.verbose) std::cout << " ... full dataset found" << std::endl;
      full_dataset = true;
    }
    else if (setup(&wf, Form("%s/combined.root", dir), waveform_tree_names, nullptr, opt.verbose) == 0)
    {
      if (opt.verbose) std::cout << " ... full dataset not found, using combined.root" << std::endl;
      partial_file = "combined";
    }
    else
    {
      std::cerr << "Failed to find waveforms.root or combined.root in " << dir << std::endl;
      return -1;
    }
  }

  setupRadiantMeta();

  /* Headers. The standalone headers.root is needed for a full dataset and
   * when the user wants to iterate over incomplete events too; otherwise the
   * headers come from the combined-style file. */
  bool want_full_headers = full_dataset || !opt.partial_skip_incomplete;
  if (opt.verbose) std::cout << "About to load headers ...";
  if (setup(&hd, Form("%s/%s.root", dir, want_full_headers ? "headers" : partial_file), header_tree_names, nullptr, opt.verbose) != 0)
  {
    // next to a combined-style file, headers.root may legitimately be absent:
    // fall back to the headers inside it, which limits us to complete events
    bool fell_back = false;
    if (!full_dataset && !opt.partial_skip_incomplete)
    {
      if (setup(&hd, Form("%s/%s.root", dir, partial_file), header_tree_names, nullptr, opt.verbose) == 0)
      {
        std::cerr << "Could not find headers.root in " << dir << "; using " << partial_file << ".root instead (only complete events)" << std::endl;
        opt.partial_skip_incomplete = true;
        fell_back = true;
      }
    }

    if (!fell_back)
    {
      std::cerr << "Failed to load headers from " << dir << std::endl;
      unload();
      return -1;
    }
  }
  if (opt.verbose) std::cout << " success" << std::endl;

  if (!full_dataset && !opt.partial_skip_incomplete)
  {
    // index the partial waveforms by event number so raw() can find them from the header entry
    wf.tree->BuildIndex("event_number");
  }

  /* DAQ status: same source selection as the headers */
  bool want_full_daqstatus = full_dataset || !opt.partial_skip_incomplete;
  if (opt.verbose) std::cout << "about to load daqstatus" << std::endl;
  if (setup(&ds, Form("%s/%s.root", dir, want_full_daqstatus ? "daqstatus" : partial_file), daqstatus_tree_names, nullptr, opt.verbose) != 0)
  {
    std::cerr << "Failed to find " << (want_full_daqstatus ? "daqstatus" : partial_file) << ".root in " << dir << std::endl;
    unload();
    return -1;
  }
  if (opt.verbose) std::cout << " success" << std::endl;

  if (full_dataset)
  {
    // index the daqstatus by readout time so status() can match it to the current event
    ds.tree->BuildIndex("int(readout_time_radiant)", "1e9*(readout_time_radiant-int(readout_time_radiant))");
  }

  /* Pedestals and run info are optional: their getters return nullptr if missing */
  if (opt.verbose) std::cout << "about to load pedestal" << std::endl;
  if (setup(&pd, Form("%s/pedestal.root", dir), pedestal_tree_names, nullptr, opt.verbose) == 0)
  {
    if (opt.verbose) std::cout << " success" << std::endl;
  }
  else
  {
    std::cerr << "Failed to find pedestal.root in " << dir << " (This is usually ok if you don't need them)" << std::endl;
  }

  if (opt.verbose) std::cout << "about to load runinfo" << std::endl;
  if (setup(&runinfo, Form("%s/runinfo.root", dir), "info", opt.verbose) == 0)
  {
    if (opt.verbose) std::cout << " success" << std::endl;
  }
  else
  {
    std::cerr << "Failed to read runinfo ..." << std::endl;
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
