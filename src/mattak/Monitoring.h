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
      // explicit default constructor
      Monitoring() : run_number(0), station_number(0) {}
      // constructor with parameters
      Monitoring(uint32_t run, uint16_t station);
      // destructor
      ~Monitoring() = default;
      // small struct for extra run parameter
      struct runParameter
      {
        std::string name;
        std::vector<float> mean;
        std::vector<float> stddev;
        std::vector<float> min;
        std::vector<float> max;
      };
      
      // Basic identifying information
      uint32_t run_number = 0;
      uint16_t station_number = 0;

      // Monitoring quantities (initial minimal set)
      uint32_t num_events = 0; // total number of events
      float rms_per_channel[24]; // root mean square per channel


      // Future-proofing: generic key-value storage
      std::vector< Monitoring::runParameter> extraParameters; // for run-level parameters


    // private:
    // everything is public for simplicity
    ClassDef(Monitoring,4);

  };

}

#endif