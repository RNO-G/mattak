#ifndef __MATTAK_DAQ_STATUS_H__
#define __MATTAK_DAQ_STATUS_H__


#include <stdint.h> 
#include "TObject.h" 

#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h" 
#else
typedef int rno_g_status_t; 
#endif

#include "mattak/Constants.h"

namespace mattak 
{
  class DAQStatus : public TObject
  {
    public: 
      DAQStatus() { ; } 
      DAQStatus(const rno_g_status_t * stat); 

      double readout_time = 0; 
      float deadtime =0; 
      
      //We will also have thresholds and such here, once we have a trigger scheme defined... 
      
      uint32_t radiant_thresholds[mattak::k::num_radiant_thresholds] = {0}; 
      uint32_t radiant_l1_rates[mattak::k::num_radiant_thresholds] = {0};  //mHz?

      uint32_t lt_thresholds[mattak::k::num_lt_thresholds] = {0}; //
      uint32_t lt_rates[mattak::k::num_lt_thresholds] = {0}; //

    ClassDef(DAQStatus,1); 
  }; 

}


#endif
