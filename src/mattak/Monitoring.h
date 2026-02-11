#ifndef __MONITORING_H__
#define __MONITORING_H__

// #pragma link C++ class mattak::Monitoring+;

// Include mattak classes
#include "mattak/Waveforms.h"
#include "mattak/Header.h"
#include "mattak/DAQStatus.h"
#include "mattak/Constants.h"

#include "TObject.h" 
#include <stdint.h>
#include <string>
#include <map>
#include <vector>

namespace mattak 
{
  
  class Monitoring : public TObject 
  {
    public:
      // default constructor
      Monitoring() = default;
      // constructor with parameters
      Monitoring(uint32_t run, uint16_t station);
      // destructor
      ~Monitoring() = default;
      struct runParameter
      {
        std::string name;
        std::vector<float> mean;
        std::vector<float> std;
        std::vector<float> min;
        std::vector<float> max;
      };
      
      // Basic identifying information
      uint32_t run_number = 0;
      uint16_t station_number = 0;

      // Monitoring quantities (initial minimal set)
      uint32_t num_events = 0; // total number of events (in the run or file?)
      double unixTime = 0;   // seconds since epoch associated with firtst event (of the run or file?)
      // system house-keeping or environmental parameters of interest for the run
      float system_voltage = 0.0f;
      float battery_level = 0.0f;
      float wind_speed = 0.0f;
      float temperature = 0.0f;
      float radiant_voltage = 0.0f;

      // Future-proofing: generic key-value storage
      std::vector<std::string> trigger_type;
      std::vector< Monitoring::runParameter> runParameters; // for run-level parameters
      std::map<std::string, std::vector<float>> eventParameters; // for event-level parameters

    // private:
    // everything is public for simplicity
    ClassDef(Monitoring,4);

  };

}

#endif