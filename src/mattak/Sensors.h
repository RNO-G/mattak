#ifndef __MATTAK_SENSORS_H__
#define __MATTAK_SENSORS_H__

#include "TObject.h"
#include <stdint.h> 

//not defined yet
typedef int rno_g_sensors_t; 
#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h" 
#else
#endif


namespace mattak 
{

  class Sensors
  {
    public: 

      Sensors() { ; } 
      Sensors(const rno_g_sensors_t * sensors); 

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

      //CPU stats 
      float disk_space_SD_card = -1; //GB
      float disk_space_internal = -1; //GB
      float free_mem = -1; // MB 
      float loadavg[3] = { 0,0,0}; 
      float cpu_uptime = 0; //seconds

      //micro stats 
      uint32_t mcu_uptime = 0; //seconds

      //LTE stats 
      float lte_rsrq = 0; 
      float lte_rsrp = 0;
      float lte_rssi = 0;
      float lte_snr = 0;

      ClassDef(Sensors,1); 
  }; 
}; 

#endif
