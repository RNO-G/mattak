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
static void clear(mattak::Dataset::tree_field<D> * field, bool borrowed = false)
{
  if (!borrowed && field->file)
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
  wf_meta.tree = (TTree*) wf_meta.file->Get(wf.tree->GetName());
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

template <typename D>
static int setup(mattak::Dataset::tree_field<D> * field, const char * filename, const char ** tree_names, const char ** branch_names = 0, bool verbose = false)
{
  clear(field);
  if (verbose) std::cout << "Trying to open " << filename << std::endl;
  field->file = !verbose ? silentlyTryToOpen(filename,"READ") : TFile::Open(filename,"READ");
  if (!field->file) return -1;

  int itry = 0;
  while(tree_names[itry])
  {
    field->tree = (TTree*) field->file->Get(tree_names[itry]);
    if (verbose) std::cout << "trying tree " << tree_names[itry] << std::endl;
    if (!field->tree)
    {
      field->tree = (TTree*) field->file->Get("combined");
      if (verbose) std::cout << "trying tree combined " << std::endl;
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
  return -1;
}

template <typename D>
static int setup(mattak::Dataset::file_field<D> * field, const char * filename, const char * obj_name, bool verbose = true)
{
  clear(field);
  field->file = !verbose ? silentlyTryToOpen(filename,"READ") : TFile::Open(filename,"READ");
  if (!field->file) return -1;
  field->ptr = (D*) field->file->Get(obj_name);
  gROOT->cd();
  return field->ptr!=nullptr;
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
    opt.base_data_dir = dir;
  else
    opt.base_data_dir = getenv("RNO_G_ROOT_DATA") ?: getenv("RNO_G_DATA") ?: ".";
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
  if (opt.verbose) std::cout << "mattak::Dataset::loadCombinedFile ( " << f  << ") called" << std::endl;
  full_dataset = false;

  if (!opt.partial_skip_incomplete)
  {
    std::cerr << "partial_skip_incomplete is incompatible with loadCombinedFile " << std::endl;
    opt.partial_skip_incomplete  = true;
  }

  if (opt.verbose) std::cout << "Opening " << f << std::endl;
  if (setup(&wf, f, waveform_tree_names, 0, opt.verbose) || setup(&hd, f, header_tree_names, 0, opt.verbose))
  {
    std::cerr << "Could not load waveforms and headers from " << f << std::endl;
    return -1;
  }

  setupRadiantMeta();

  if (opt.verbose) std::cout << "Found waveforms and headers in" << f << std::endl;

  // Try some optionalish things
  if (setup(&ds, f, daqstatus_tree_names, 0, opt.verbose))
  {
    std:: cerr << "Could not load daqstatus from " << f << " (this is ok if you don't use them) " << std::endl;
  }
  else
  {
    if (opt.verbose) std::cout << "Found daqstatus in" << f << std::endl;
  }

  // we probably don't have pedetals, but we could try I guess?
  if (!setup(&pd, f, pedestal_tree_names, 0, opt.verbose))
  {
    if (opt.verbose) std::cout << "Found pedestals in" << f << std::endl;
  }

  if ( !(setup(&runinfo, f, "info") || setup(&runinfo, f, "runinfo")) )
  {
    if (opt.verbose) std::cout << "Found runinfo in" << f << std::endl;
  }

  return 0;
}

int mattak::Dataset::loadDir(const char * dir)
{

  if (opt.verbose) std::cout << "mattak::Dataset::loadDir (" << dir  << ", skip_incomplete=" << opt.partial_skip_incomplete << ") called" << std::endl;

  //first clear all
  unload();
  current_entry = 0;

  if (opt.verbose) std::cout << "Load waveforms ..." << std::endl;

  const char * partial_file = NULL;
  if (opt.file_preference != "" && !setup(&wf, Form("%s,%s.root",dir,opt.file_preference.c_str()), waveform_tree_names))
  {
    full_dataset = false;
    partial_file = opt.file_preference.c_str();
  }
  else
  {
    if (opt.file_preference != "")
    {
      std::cerr << "Warning, could not find preferred %s.root in %s. Reverting to default behavior" << std::endl;
    }

    //we need to figure out if this is a full run or partial run, so check for existence of waveforms.root
    if (setup(&wf, Form("%s/waveforms.root", dir), waveform_tree_names, 0))
    {
      //no waveforms file!
      full_dataset = false;
      if (opt.verbose) std::cout << " ... full dataset not found " << std::endl;

      //let's load from combined file instead
      if (setup(&wf, Form("%s/combined.root", dir), waveform_tree_names, nullptr, opt.verbose))
      {
        //uh oh, we didn't find it there either :(
        std::cerr << "Failed to find waveforms.root or combined.root in " << dir << std::endl;
        return -1;
      }

      partial_file = "combined";
    }
    else
    {
      if (opt.verbose) std::cout << " ... full dataset found " << std::endl;
      full_dataset = true;
    }
  }

  setupRadiantMeta();
 
  //now load the header files
  if (opt.verbose) std::cout << "about to load headers " << std::endl;
  if (setup(&hd, Form("%s/%s.root", dir, (full_dataset || !opt.partial_skip_incomplete) ? "headers" : partial_file), header_tree_names, nullptr, opt.verbose))
  {
    std::cerr << "Failed to find headers.root or " << partial_file << " .root in " << dir << std::endl;
    return -1;
  }
  if (opt.verbose) std::cout << " success" << std::endl;

  if (!full_dataset && !opt.partial_skip_incomplete)
  {
    //set up an index on event number the events
    wf.tree->BuildIndex("event_number");
  }

  //and the status files
  if (opt.verbose) std::cout << "about to load daqstatus " << std::endl;
  if (setup(&ds, Form("%s/%s.root", dir, full_dataset || !opt.partial_skip_incomplete ? "daqstatus" : partial_file), daqstatus_tree_names, nullptr, opt.verbose))
  {
    std::cerr << "Failed to find daqstatus.root or " << partial_file << " in " << dir << std::endl;
    return -1;
  }
  if (opt.verbose) std::cout << " success" << std::endl;

  if (full_dataset)
  {
    ds.tree->BuildIndex("int(readout_time_radiant)", "1e9*(readout_time_radiant-int(readout_time_radiant))");
  }

  //and the pedestal files
  if (opt.verbose) std::cout << "about to load pedestal " << std::endl;
  if (setup(&pd, Form("%s/pedestal.root", dir), pedestal_tree_names, nullptr, opt.verbose))
  {
    std::cerr << "Failed to find pedestal.root in " << dir << " (This is usually ok if you don't need them) " << std::endl;
  }
  else
  {
  	if (opt.verbose) std::cout << " success" << std::endl;
  }

  //and try the runinfo file
  if (opt.verbose) std::cout << "about to load runinfo " << std::endl;
  if (setup(&runinfo, Form("%s/runinfo.root", dir), "info",  opt.verbose))
  {
     std::cerr << "Failed to read runinfo ... " << std::endl;
  }
  else
  {
    if (opt.verbose) std::cout << " success" << std::endl;
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
  if (wf_meta.ptr == nullptr) return (info() && info()->radiant_sample_rate) ? info()->radiant_sample_rate : 3200;
  findIncompleteEntry(&wf_meta,this, force, &wf);
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
  findIncompleteEntry(&wf_meta,this, force, &wf);
  return (const float*)  (  wf_meta.missing_entry ? zeros : wf_meta.ptr->digitizer_readout_delay_ns );
}


mattak::Waveforms* mattak::Dataset::raw(bool force)
{
  if (wf.tree == nullptr) return nullptr; 
  findIncompleteEntry(&wf, this, force);
  return wf.missing_entry ? nullptr: wf.ptr;
}


mattak::DAQStatus * mattak::Dataset::status(bool force)
{
  if (!ds.tree) return nullptr;
  if (force || ds.loaded_entry != current_entry)
  {
    if (full_dataset)
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
  if (entry < 0 || entry > pd.tree->GetEntries()) return nullptr;

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
