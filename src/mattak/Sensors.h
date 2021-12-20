#ifndef __MATTAK_SENSORS_H__
#define __MATTAK_SENSORS_H__

#include "TObject.h"
#include <stdint.h> 


namespace mattak 
{



  class Sensors
  {
    public: 

      Sensors() { ; } 


      //Note this is approximate, not all of these are necessarily collected at the same time
      //May have multiple times in the future... 
      double readout_time = 0; 

      float PV_voltage = 0;   //V 
      float PV_current = 0;   //mA
      float battery_voltage = 0; //V
      float battery_current = 0; 

      float T_power_board = -273;  // C
      float T_case =-273;  // C
      float T_amps =-273;  // C
      float T_controller_board = -273;  // C

      float radiant_current = 0;  // mA
      float lt_current = 0; // mA
      float sbc_current = 0;  // mA

      //TODO define mappings for these (mA) 
      float downhole_currents[3] = {0}; 
      float amp_currents[6] = {0}; 


      static int makeTreeFromDatabase(const char * outfile, const char * pgconninfo, int station, int start_time, int end_time, const char * tree_name = "sensors"); 
      static int makeInterpolatedEventTreeFromDatabase(const char * outfile, const char * pgconninfo, const char * header_file, const char * tree_name = "sensors"); 


      ClassDef(Sensors,1); 
  }; 

  class LTEStats
  {

      //LTE stats 
      double readout_time =0; 
      float lte_rsrq = 0; 
      float lte_rsrp = 0;
      float lte_rssi = 0;
      float lte_snr = 0;
      ClassDef(LTEStats,1); 
  }; 

}; 

#endif
