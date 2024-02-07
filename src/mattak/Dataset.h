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


        bool verbose = false; 
      }; 


  class Dataset
  {

    public:
      Dataset(int station, int run, const DatasetOptions & opt = DatasetOptions()); 
      Dataset(const DatasetOptions & opt = DatasetOptions()); 
      //data_dir defaults to RNO_G_ROOT_DATA
      void setVerbose(bool v) { opt.verbose = v; } 
      virtual ~Dataset() { unload() ; }

      void setOpt(const DatasetOptions & opt = DatasetOptions()); 
      const DatasetOptions & getOpt() const { return opt; } 

      int loadRun(int station, int run, const DatasetOptions & opt) ; 
      int loadDir(const char * dir, const DatasetOptions & opt); 
      int loadCombinedFile(const char * file, const DatasetOptions & opt) ; 

      int loadRun(int station, int run) ; 
      int loadDir(const char * dir); 
      int loadCombinedFile(const char * file); 


      /** 
       * Deprecated, kept for ABI compatibility
       */
      Dataset (int station, int run, const VoltageCalibration * calib, const char * base_data_dir = nullptr, bool partial_skip_incomplete = true, bool verbose = false); 
      Dataset (const char * data_dir); 
      int loadRun(int station, int run, bool partial_skip_incomplete); 
      int loadDir(const char * dir, bool partial_skip_incomplete ); 
      void setDataDir(const char * dir); 


      bool setEntry(int entry); //returns true if in range 
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
      tree_field<Header> hd; 
      tree_field<DAQStatus> ds; 
      tree_field<Pedestals> pd; 
      file_field<RunInfo> runinfo; 
      field<CalibratedWaveforms> calib_wf; 

      void unload(); 
      int current_entry = 0; 

      const VoltageCalibration * calib = nullptr; 
      bool full_dataset ; 
      DatasetOptions opt; 

  }; 
}

#endif 
