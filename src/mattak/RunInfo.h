#ifndef _MATTAK_RUNINFO_H
#define _MATTAK_RUNINFO_H 

#include "TTimeStamp.h" 
#include "TDatime.h" 
#include "TObject.h" 
#include <string> 
#include <unordered_map> 
#include <mattak/Constants.h> 




namespace mattak 
{

  struct FirmwareVersion : public TObject
  {
    uint8_t major = 0; 
    uint8_t minor = 0; 
    uint8_t rev = 0; 
    uint8_t day = 0; 
    uint8_t month = 0; 
    uint16_t year = 0; 
    ClassDef(FirmwareVersion,1); 
  }; 


  struct FlowerGainCode : public TObject 
  {
    uint8_t codes[k::num_lt_channels]  = {0}; 
    time_t when = 0; 
    int station =0; 
    int run = 0; 
    ClassDef(FlowerGainCode,1); 
  }; 

  class RunInfo : public TObject
  {

    public: 
      RunInfo()  {;}
      RunInfo(const char * auxdir); 
      int station = 0; 
      int run = 0; 

      TTimeStamp run_start_time = 0; 
      TTimeStamp run_end_time = 0; 

      float radiant_sample_rate = 0; 
      FirmwareVersion radiant_fpga;
      FirmwareVersion radiant_bm;
      FirmwareVersion flower; 

      std::string librnog_version = "";
      std::string daq_version = "";
      std::string mattak_version = ""; 
      std::string comment =""; 

      float MB_free_data_partition = 0;
      float MB_free_main_partition = 0;

      //Returns "" if none! 
      const std::string & lookup(const std::string & key) const;

      std::vector<FlowerGainCode> flower_codes; 

      //looks up the int, filling val. returns 0 if successful (return value is not the value). 
      int lookupInt(const std::string & key, int *val, int base =0) const; 
      int lookupFloat(const std::string & key, float *val) const; 
      int lookupTimeStamp(const std::string & key, TTimeStamp * v) const; 
      int lookupFirmwareVersion(const std::string & verkey, const std::string & datekey, FirmwareVersion *val) const; 

      
    private: 
      std::unordered_map<std::string,std::string> kvp; 
      ClassDef(RunInfo,1); 
  }; 



}




#endif
