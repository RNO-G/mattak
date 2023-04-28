#include "mattak/Dataset.h" 
#include "TSystem.h" 
#include "TROOT.h" 
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
  clear ((mattak::Dataset::field<D>*) field); 
}

template <typename D> 
static void clear(mattak::Dataset::file_field<D> * field)
{
  if (field->file) 
  {
    delete field->file; 
  }
  field->file = nullptr; 

  clear ((mattak::Dataset::field<D>*) field); 

 
}


#ifdef __WIN32
#define BITBUCKET "NUL" 
#else
#define BITBUCKET "/dev/null" 
#endif

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
static int setup(mattak::Dataset::tree_field<D> * field, const char * filename, const char ** tree_names, const char ** branch_names = 0) 
{
  clear(field); 
  field->file = silentlyTryToOpen(filename,"READ"); 
  if (!field->file) return -1; 

  int itry = 0; 
  while(tree_names[itry]) 
  {
    field->tree = (TTree*) field->file->Get(tree_names[itry]); 
    if (!field->tree) 
    {
      field->tree = (TTree*) field->file->Get("combined"); 
    }

    if (!field->tree) 
    {
      itry++; 
      continue; 
    }

    const char * branch_name = branch_names ? branch_names[itry] : tree_names[itry]; 
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

    if (!branch_name[0]) gSystem->RedirectOutput(0,"a",&rh); 

    gROOT->cd(); 
    return 0; 
 }
 std::cerr << "Could not find a valid tree/branch pair in " << filename << std::endl; 
 return -1; 

}

template <typename D> 
static int setup(mattak::Dataset::file_field<D> * field, const char * filename, const char * obj_name) 
{
  clear(field); 
  field->file = silentlyTryToOpen(filename,"READ"); 
  if (!field->file) return -1; 
  field->ptr = (D*) field->file->Get(obj_name); 
  gROOT->cd();
  return field->ptr!=nullptr; 
}


void mattak::Dataset::unload() 
{
  clear(&wf); 
  clear(&hd); 
  clear(&ds); 
  clear(&pd); 
  clear(&runinfo); 
  clear(&calib_wf); 
}


mattak::Dataset::Dataset(int station, int run, const VoltageCalibration * calib, const char * data_dir, bool partial_skip, bool v) 
{
  setVerbose(v);  // should be first
  setDataDir(data_dir); 
  setCalibration(calib); 
  loadRun(station,run,partial_skip);
}

mattak::Dataset::Dataset(const char* data_dir) 
{
  setDataDir(data_dir); 
}


void mattak::Dataset::setDataDir(const char * dir) 
{
  if (dir) data_dir = dir; 
  else
  {
    data_dir = getenv("RNO_G_ROOT_DATA")
              ?: getenv("RNO_G_DATA") 
                ?: "."; 
  }
}

void mattak::Dataset::setCalibration(const VoltageCalibration * c) 
{
  calib = c; 
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
  TString dir;
  dir.Form("%s/station%d/run%d", data_dir.c_str(), station, run); 
  return loadDir(dir.Data(), partial_skip); 
}

int mattak::Dataset::loadDir(const char * dir, bool partial_skip) 
{

  if (verbose) std::cout << "mattak::Dataset::loadDir ( " << dir  << "," << partial_skip << ") called" << std::endl;
  skip_incomplete = partial_skip; 

  //first clear all 
  unload(); 
  current_entry = 0; 

  if (verbose) std::cout << "about to load waveforms " << std::endl; 
  //we need to figure out if this is a full run or partial run, so check for existence of waveforms.root
  if (setup(&wf, Form("%s/waveforms.root", dir), waveform_tree_names, 0))
  {
    //no waveforms file! 
    full_dataset = false;
    if (verbose) std::cout << " full dataset not found " << std::endl;

    //let's load from combined file instead
    if (setup(&wf, Form("%s/combined.root", dir), waveform_tree_names))
    {
      //uh oh, we didn't find it there either :( 
      std::cerr << "Failed to find waveforms.root or combined.root in " << dir << std::endl; 
      return -1; 
    }
  }
  else
  {
    if (verbose) std::cout << " full dataset found " << std::endl;
    full_dataset = true; 
  }


  if (verbose) std::cout << "about to load headers ..."; 
  //now load the header files 
  if ( setup(&hd, 
      Form("%s/%s.root", dir, (full_dataset || !partial_skip) ? "headers" : "combined"), 
      header_tree_names) )
  {
    std::cerr << "Failed to find headers.root or combined.root in " << dir << std::endl; 
    return -1; 
  } 
  if (verbose) std::cout << " success" << std::endl; 

  if (!full_dataset && !partial_skip)
  {
    //set up an index on event number the events
    wf.tree->BuildIndex("event_number"); 
  }

  if (verbose) std::cout << "about to load daqstatus ..."; 
  //and the status files
  if ( setup(&ds, 
      Form("%s/%s.root", dir, full_dataset || !partial_skip ? "daqstatus" : "combined"), 
      daqstatus_tree_names) )
  {
    std::cerr << "Failed to find daqstatus.root or combined.root in " << dir << std::endl; 
    return -1; 
  }
  if (verbose) std::cout << " success" << std::endl; 

  if (full_dataset) 
  {
    ds.tree->BuildIndex("readout_time_radiant"); 
  }

  if (verbose) std::cout << "about to load pedestal ..."; 
  //and the pedestal files
  if ( setup(&pd, 
      Form("%s/pedestal.root", dir), 
      pedestal_tree_names) )
  {
    std::cerr << "Failed to find pedestal.root in " <<dir << std::endl; 
    return -1; 
  } 
  if (verbose) std::cout << " success" << std::endl; 

  if (verbose) std::cout << "about to load runinfo ..."; 
  //and try the runinfo file 
  setup(&runinfo, Form("%s/runinfo.root", dir),"info"); 

  if (verbose) std::cout << " success" << std::endl; 
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

mattak::Waveforms* mattak::Dataset::raw(bool force) 
{
  if (force || wf.loaded_entry != current_entry)
  {
    if (wf.tree == nullptr) return nullptr; 

    if (full_dataset || skip_incomplete) 
    {
      wf.branch->GetEntry(current_entry); 
    }
    else
    {
      int wf_entry = wf.tree->GetEntryNumberWithIndex(header(force)->event_number); 
      if (wf_entry < 0) 
      {
        wf.missing_entry = true; 
      }
      else
      {
        wf.missing_entry = false; 
        wf.branch->GetEntry(wf_entry); 
      }

    }

    wf.loaded_entry = current_entry; 
  }

  return wf.missing_entry ? nullptr: wf.ptr; 
}


mattak::DAQStatus * mattak::Dataset::status(bool force) 
{
  if (force || ds.loaded_entry != current_entry) 
  {
    if (full_dataset) 
    {
      int ds_entry = ds.tree->GetEntryNumberWithBestIndex(header(force)->readout_time); 
      if (ds_entry < 0) 
      {
        ds.missing_entry = true; 
      }
      else
      {
        ds.branch->GetEntry(ds_entry); 
        ds.missing_entry = false; 

      }
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
  if (!calib) return nullptr; 

  if (force || calib_wf.loaded_entry != current_entry) 
  {
    mattak::Waveforms * r = raw(force); 
    mattak::Header * h = header(force); 

    // if tehre is no raw waveform, we can't do this
    if (!r|| !h) 
    {
      calib_wf.missing_entry = true; 
    } 
    else
    {
      calib_wf.missing_entry = false; 

      if (!calib_wf.ptr) 
      {
        calib_wf.ptr = new CalibratedWaveforms(*r, *h, *calib); 
      }
      else 
      {
        new (calib_wf.ptr) CalibratedWaveforms(*r, *h, *calib); 
      }
    }

    calib_wf.loaded_entry = current_entry; 
  }

  return calib_wf.missing_entry ? nullptr : calib_wf.ptr; 
}


mattak::Pedestals * mattak::Dataset::peds(bool force, int entry)
{
  if (entry < 0 || entry > pd.tree->GetEntries()) return nullptr; 

  if (force || entry != pd.loaded_entry) 
  {
    pd.branch->GetEntry(entry); 
    pd.loaded_entry = entry; 
  }
  return pd.ptr; 
}
