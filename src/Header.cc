#include "mattak/Header.h" 
#include <iostream> 
#include <cmath>


ClassImp(mattak::Header); 


mattak::Header::Header(const rno_g_header_t * head) 
  : mattak::Header() 
{

#ifdef LIBRNO_G_SUPPORT

  this->run_number = head->run_number; 
  this->event_number = head->event_number; 
  this->trigger_number = head->trigger_number; 
  this->station_number = head->station_number; 
  this->buffer_length = head->radiant_nsamples; 
  this->pretrigger_samples= head->pretrigger_windows*128; 
  this->readout_time = head->readout_time_secs + 1e-9 * head->readout_time_nsecs; 
  this->pps_num = head->pps_count; 
  this->sysclk = head->sys_clk; 
  this->sysclk_last_pps = head->sysclk_last_pps; 
  this->sysclk_last_last_pps = head->sysclk_last_last_pps; 

  //subsecond part
  double sysclk_diff = this->sysclk - this->sysclk_last_pps; 
  const double two_to_the_32 = 4294967296; 
  if (sysclk_diff < 0) sysclk_diff += two_to_the_32;  
  double last_sysclk_diff = this->sysclk_last_pps - this->sysclk_last_last_pps; 
  if (last_sysclk_diff < 0)  last_sysclk_diff += two_to_the_32; 
  this->trigger_time = sysclk_diff / last_sysclk_diff; 

  //readout time is always afer trigger time, so figure out the second based on what's closet
  if (1e-9 * head->readout_time_nsecs < this->trigger_time) this->trigger_time += head->readout_time_secs-1; 
  else this->trigger_time += head->readout_time_secs; 

  this->trigger_info.force_trigger = !!(head->trigger_type & RNO_G_TRIGGER_SOFT); 
  this->trigger_info.pps_trigger = !!(head->trigger_type & RNO_G_TRIGGER_PPS); 
  this->trigger_info.rf_trigger = ! (this->trigger_info.force_trigger || this->trigger_info.pps_trigger); 
  this->trigger_info.radiant_trigger = !!(head->trigger_type & (RNO_G_TRIGGER_RF_RADIANTX)); 
  this->trigger_info.lt_trigger = !!(head->trigger_type & (RNO_G_TRIGGER_RF_LT_SIMPLE | RNO_G_TRIGGER_RF_LT_PHASED)); 
  if (this->trigger_info.radiant_trigger) 
  {
    this->trigger_info.which_radiant_trigger =  
    (head->trigger_type & RNO_G_TRIGGER_RF_RADIANT0) ? 0 : 
    (head->trigger_type & RNO_G_TRIGGER_RF_RADIANT1) ? 1 : 
    -1; 
  }
  else this->trigger_info.which_radiant_trigger=-128; 


  for (int i = 0 ; i < RNO_G_NUM_RADIANT_CHANNELS; i++) 
  {
    this->trigger_info.radiant_info.start_windows[i][0] = head->radiant_start_windows[i][0]; 
    this->trigger_info.radiant_info.start_windows[i][1] = head->radiant_start_windows[i][1]; 
  }

  for (int i = 0; i  < 2; i++) 
  {
    this->trigger_info.radiant_info.RF_masks[i] = head->radiant_trigger_cfg[i].mask; 
    this->trigger_info.radiant_info.RF_ncoinc[i] = head->radiant_trigger_cfg[i].ncoinc; 
    this->trigger_info.radiant_info.RF_window[i] = head->radiant_trigger_cfg[i].window; 
  }
     
  this->trigger_info.lt_info.window = head->lt_simple_trigger_cfg.window;
  this->trigger_info.lt_info.num_coinc = head->lt_simple_trigger_cfg.num_coinc;
  this->trigger_info.lt_info.vppmode = head->lt_simple_trigger_cfg.vpp_mode;



#else
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) head; 
#endif

}
