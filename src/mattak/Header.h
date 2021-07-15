#ifndef __MATTAK_HEADER_H__
#define __MATTAK_HEADER_H__

#include <stdint.h> 
#include "mattak/TriggerInfo.h" 
#include "TObject.h" 

#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h" 
#else
typedef int rno_g_header_t; 
#endif

namespace mattak 
{


  class Header : public TObject 
  {

    public: 
      Header() { ; }; 
      Header(const rno_g_header_t * hdr); 

      uint32_t run_number = 0; 
      uint32_t event_number = 0; 
      uint32_t trigger_number = 0;   // May differ from event_number due to deadtime 
      uint16_t station_number = 0; 
      uint16_t buffer_length = 0; 
      uint16_t pretrigger_samples = 0; 
      double readout_time = 0; 
      uint32_t pps_num= 0;
      uint32_t sysclk= 0;
      uint32_t sysclk_last_pps= 0;
      uint32_t sysclk_last_last_pps= 0;
      double trigger_time; 

      TriggerInfo trigger_info; 


    ClassDef(Header,1); 
  }; 



}


#endif 
