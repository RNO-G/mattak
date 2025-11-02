#ifndef __MONITORING_H__
#define __MONITORING_H__

// #pragma link C++ class mattak::Monitoring+;

// #include "Header.h"
// #include "DAQStatus.h"

#include "TObject.h" 
#include <stdint.h>
#include <string>

namespace mattak 
{
  class Monitoring : public TObject 
  {
    // struct MonitoringData 
    // {
    //   int run_number;
    //   double radiant_voltage;
    //   double battery_voltage;
    //   double wind_speed;
    // };
    public:
      Monitoring(){;};
      Monitoring(uint32_t run, uint16_t station, std::string RNOG_DATA_DIR);
      ~Monitoring(){;};
      // run number 
      uint32_t run_number = 0;
      uint16_t station_number = 0;
    ClassDef(Monitoring,1);

  };

}

#endif