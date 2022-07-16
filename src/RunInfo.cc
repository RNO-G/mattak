#include "mattak/RunInfo.h" 
#include <fstream> 
#include "mattak/Version.h" 
#include <stdlib.h> 
#include "TString.h"
#include <iostream> 

ClassImp(mattak::RunInfo); 




void trim(std::string & s)
{
  const char * ws = " \t"; 
  s.erase(s.find_last_not_of(ws)+1); 
  s.erase(0, s.find_first_not_of(ws)); 


}

mattak::RunInfo::RunInfo(const char * auxdir)
{
  mattak_version = mattak::version(); 

  //TODO: use fmt? 
  std::ifstream fruninfo(Form("%s/runinfo.txt",auxdir)); 

  if (fruninfo.good())
  {
    std::string line; 


    //First fill the kvp 
    while (std::getline(fruninfo,line))
    {
    
      std::string key; 
      std::string value; 
      //find = 
      
      size_t where = line.find("="); 
    
      //Skip lines without 
      if (where == std::string::npos) 
      {
        continue; 
      }

      key = line.substr(0,where); 
      value = line.substr(where+1); 

      //trim leading and trailing whitespace
      trim(key); 
      trim(value); 
      std::cout <<key <<":" << value << std::endl; 
      kvp[key] = value; 
    }


    //now let's fill in things we can find 
    librnog_version = lookup("LIB-RNO-G-GIT-HASH"); 
    daq_version = lookup("RNO-G-ICE-SOFTWARE-GIT-HASH"); 
    lookupInt("STATION",&station);
    lookupInt("RUN",&run);
    lookupFloat("RADIANT-SAMPLERATE",&radiant_sample_rate); 
    lookupFloat("FREE-SPACE-MB-OUTPUT-PARTITION", &MB_free_data_partition);
    lookupFloat("FREE-SPACE-MB-RUNFILE-PARTITION", &MB_free_main_partition);
    lookupTimeStamp("RUN-START-TIME", &run_start_time);
    lookupTimeStamp("RUN-END-TIME", &run_end_time);
    lookupFirmwareVersion("RADIANT-FWVER", "RADIANT-FWDATE", &radiant_fpga);
    lookupFirmwareVersion("RADIANT-BM-FWVER", "RADIANT-BM-FWDATE", &radiant_bm);
    lookupFirmwareVersion("FLOWER-FWVER", "FLOWER-FWDATE", &flower);
  }
  else
  {
    std::cerr << "Could not open runinfo.txt" << std::endl; 
  }

  //now look for flower gain codes 

  int iflower = 0; 
  while(true) 
  {
    std::ifstream ifs(Form("%s/flower_gain_codes.%d.txt", auxdir, iflower++)); 
    if (!ifs.good()) break; 
    FlowerGainCode code; 
    std::string line; 
    std::getline(ifs,line); 
    sscanf(line.c_str(),"# Flower gain codes, station=%d, run=%d,  time=%d", &code.station,&code.run, &code.when); 
    std::getline(ifs,line); 
    sscanf(line.c_str(),"%hhu %hhu %hhu %hhu", &code.codes[0],&code.codes[1], &code.codes[2], &code.codes[3]); 
    flower_codes.push_back(code); 
  }

  //now look for comment 
  std::ifstream ifscomment(Form("%s/comment.txt", auxdir), std::ios::ate); 
  if (ifscomment.good())
  {
    auto size = ifscomment.tellg(); 
    comment.resize(size); 
    ifscomment.seekg(0); 
    ifscomment.read(&comment[0], size);
  }
}

int mattak::RunInfo::lookupInt(const std::string & key, int * val, int base) const
{
  std::string str = lookup(key); 
  if (str=="") 
  {
    return 1; 
  }

  int maybe_val = 0; 
  if (1!=sscanf(str.c_str(),base == 16 ? "%x" : base == 10 ? "%d" : base == 8 ? "%o" : "%i", &maybe_val)) 
  {
    return 1; 
  }

  else if (val) *val = maybe_val; 
  return 0 ; 
}

int mattak::RunInfo::lookupFloat(const std::string & key, float * val) const
{
  std::string str = lookup(key); 
  if (str=="") 
  {
    return 1; 
  }
  float tmp; 

  if  (1!=sscanf(str.c_str(),"%f", &tmp))
  {
    return 1; 
  }

  if (val) *val = tmp; 

  return 0; 

}

int mattak::RunInfo::lookupTimeStamp(const std::string & key, TTimeStamp * val) const
{
  std::string str = lookup(key); 
  if (str=="") 
  {
    return 1; 
  }

  int sec, nsec; 
  if (2!=sscanf(str.c_str(), "%d.%09d", &sec, &nsec))
  {
    return 1; 
  }
  if (val) 
  {
    val->SetSec(sec);
    val->SetNanoSec(nsec);
  }
  return 0; 
}

int mattak::RunInfo::lookupFirmwareVersion(const std::string & verkey, const std::string & datekey, FirmwareVersion *val) const
{

  std::string ver = lookup(verkey); 
  std::string date = lookup(datekey); 
  if (ver=="" || date =="") return 1;
  if (!val) return 1; 


  return 3 != sscanf(ver.c_str(), "%02hhu.%02hhu.%02hhu", &val->major, &val->minor, &val->rev) 
  || 3 != sscanf(date.c_str(), "%hu-%02hhu.%02hhu", &val->year, &val->month, &val->day); 

}

static std::string empty(""); 
const std::string & mattak::RunInfo::lookup(const std::string & key) const
{
  auto it = kvp.find(key); 
  if (it == kvp.end()) return empty; 
  else return it->second;
}
