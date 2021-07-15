#ifndef __MATTAK_DAQ_STATUS_H__
#define __MATTAK_DAQ_STATUS_H__


#include <stdint.h> 
#include "TObject.h" 

#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h" 
#else
typedef int rno_g_daqstatus_t; 
#endif

#include "mattak/Constants.h"



namespace mattak 
{

  struct LTScalerGroup
  {
    uint16_t trig_coinc = 0;
    uint16_t trig_per_chan[mattak::k::num_lt_channels] = {0};
    uint16_t servo_coinc = 0;
    uint16_t servo_per_chan[mattak::k::num_lt_channels] = {0};
  }; 

  class DAQStatus : public TObject
  {
    public: 
      DAQStatus() { ; } 
      DAQStatus(const rno_g_daqstatus_t * stat); 

      double readout_time_radiant = 0; 
      double readout_time_lt = 0; 
      
      //We will also have thresholds and such here, once we have a trigger scheme defined... 
      
      uint32_t radiant_thresholds[mattak::k::num_radiant_channels] = {0}; 
      uint32_t radiant_scalers[mattak::k::num_radiant_channels] = {0};  //mHz?
      uint8_t radiant_prescalers_m1[mattak::k::num_radiant_channels] = {0};  //mHz?
      float radiant_scaler_period;

      uint32_t lt_trigger_thresholds[mattak::k::num_lt_channels] = {0}; //
      uint32_t lt_servo_thresholds[mattak::k::num_lt_channels] = {0}; //
      LTScalerGroup lt_1Hz_scalers;
      LTScalerGroup lt_1Hz_gated_scalers;
      LTScalerGroup lt_100Hz_scalers;
      uint64_t lt_ncycles = 0; 
      uint16_t lt_scaler_counter = 0; 
      uint8_t station_number = 0; 
    ClassDef(DAQStatus,2); 
  }; 

}


#endif
