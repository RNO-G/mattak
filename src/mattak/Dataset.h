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
  class Dataset
  {

    public:
      //data_dir defaults to RNO_G_ROOT_DATA
      Dataset (int station, int run, const VoltageCalibration * calib = nullptr, const char * data_dir = nullptr, bool partial_skip_incomplete = true); 
      Dataset (const char * data_dir = nullptr); 
      void setVerbose(bool v) { verbose = v; } 
      virtual ~Dataset() { unload() ; }
      /** Loads run corresponding to station/run" in the data dir, i.e. equivalent to calling loadDir(data_dir/stationS/runR) */ 
      int loadRun(int station, int run, bool partial_skip_incomplete = true); 
      int loadDir(const char * dir, bool partial_skip_incomplete = true); 
      void setDataDir(const char * dir); 

      void setEntry(int entry); 
      int N() const; 

      mattak::Header * header(bool force_reload = false); 
      mattak::Waveforms * raw(bool force_reload = false); 
      mattak::CalibratedWaveforms * calibrated(bool force_reload = false); //will be nullptr if no calibration is passed 
      mattak::DAQStatus * status(bool force_reload = false); 
      mattak::RunInfo * info() const { return runinfo.ptr; }
      mattak::Pedestals * peds(bool force_reload = false, int entry = 0); 

      TTree * daqStatusTree() { return ds.tree; }
      TTree * headTree() { return hd.tree; }
      TTree * wfTree() { return wf.tree; }

      bool isFullDataset() const { return full_dataset; }
      void setCalibration(const VoltageCalibration * calib); 
      const VoltageCalibration * getCalibration() const { return calib; } 


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
      tree_field<Header> hd; 
      tree_field<DAQStatus> ds; 
      tree_field<Pedestals> pd; 
      file_field<RunInfo> runinfo; 
      field<CalibratedWaveforms> calib_wf; 

      void unload(); 
      int current_entry = 0; 

      const VoltageCalibration * calib = nullptr; 
      std::string data_dir; 
      bool full_dataset ; 
      bool skip_incomplete;
      bool verbose = false; 

  }; 
}

#endif 
