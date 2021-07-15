#include "mattak/DAQStatus.h" 
#include <iostream> 


ClassImp(mattak::DAQStatus); 


mattak::DAQStatus::DAQStatus(const rno_g_daqstatus_t * status) 
  : DAQStatus() 
{

#ifndef LIBRNO_G_SUPPORT
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) status; 
#else

  this->readout_time_radiant = status->when_radiant; 
  this->readout_time_lt = status->when_lt; 

  for (int i = 0; i < mattak::k::num_radiant_channels; i++) 
  {
    this->radiant_thresholds[i] = status->radiant_thresholds[i]; 
    this->radiant_scalers[i] = status->radiant_scalers[i]; 
    this->radiant_prescalers_m1[i] = status->radiant_prescalers[i]; 
  }

  this->radiant_scaler_period = status->radiant_scaler_period; 

  for (int i = 0; i < mattak::k::num_lt_channels; i++) 
  {
    this->lt_trigger_thresholds[i] = status->lt_trigger_thresholds[i];
    this->lt_servo_thresholds[i] = status->lt_servo_thresholds[i];
    this->lt_1Hz_scalers.trig_per_chan[i] = status->lt_scalers.s_1Hz.trig_per_chan[i];
    this->lt_1Hz_scalers.servo_per_chan[i] = status->lt_scalers.s_1Hz.servo_per_chan[i];
    this->lt_1Hz_gated_scalers.trig_per_chan[i] = status->lt_scalers.s_1Hz_gated.trig_per_chan[i];
    this->lt_1Hz_gated_scalers.servo_per_chan[i] = status->lt_scalers.s_1Hz_gated.servo_per_chan[i];
    this->lt_100Hz_scalers.trig_per_chan[i] = status->lt_scalers.s_100Hz.trig_per_chan[i];
    this->lt_100Hz_scalers.servo_per_chan[i] = status->lt_scalers.s_100Hz.servo_per_chan[i];
  }

  this->lt_ncycles = status->lt_scalers.ncycles;
  this->lt_scaler_counter = status->lt_scalers.scaler_counter_1Hz;
  this->lt_1Hz_scalers.trig_coinc = status->lt_scalers.s_1Hz.trig_coinc;
  this->lt_1Hz_scalers.servo_coinc = status->lt_scalers.s_1Hz.servo_coinc;
  this->lt_1Hz_gated_scalers.trig_coinc = status->lt_scalers.s_1Hz_gated.trig_coinc;
  this->lt_1Hz_gated_scalers.servo_coinc = status->lt_scalers.s_1Hz_gated.servo_coinc;
  this->lt_100Hz_scalers.trig_coinc = status->lt_scalers.s_100Hz.trig_coinc;
  this->lt_100Hz_scalers.servo_coinc = status->lt_scalers.s_100Hz.servo_coinc;




#endif

}
