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
      /** Construct an empty Header object */
      Header() { ; };

      /** Construct a Header object from an rno_g_header_t (the equivalent thing in the DAQ data format) */
      Header(const rno_g_header_t * hdr);

      /* Number of this run */
      uint32_t run_number = 0;

      /* Number of this event */
      uint32_t event_number = 0;

      /* Currently unused, in future may be used to keep track of missed triggers due to deadtime */
      uint32_t trigger_number = 0;

      /* The number of this station */
      uint16_t station_number = 0;

      /* The number of samples in the waveforms */
      uint16_t buffer_length = 0;

      /* In theory, the number of samples pretrigger, but not yet filled properly. */
      uint16_t pretrigger_samples = 0;

      /* The readout time, as a UTC double. This is the time the event made it to the SBC */
      double readout_time = 0;

      /* The number of PPS (pulse per second) received since the start of the run by the RADIANT. WARNING: this can slip */
      uint32_t pps_num= 0;

      /* The number of cycles in the nominally 100 MHz clock for the current event. Note that this wraps every 53 seconds or so.
       * WARNING: this can slip.
       * */
      uint32_t sysclk= 0;

      /** The number of cycles in the nominally 100 MHz clock at the time of the last PPS
       * WARNING: this can slip.
       * */
      uint32_t sysclk_last_pps= 0;

      /** The number of cycles in the nominally 100 MHz clock at the time of the  PPS previous to last
       * WARNING: this can slip.
       * */
      uint32_t sysclk_last_last_pps= 0;

      /* The trigger time, as a UTC double. Note that this does not have as much precision as is possible,
       * which should be fixed in the future (though you can rederive from sysclk).
       *
       * This is  currently derived from readout_time, sysclk, sysclk_last_pps and sysclk_last_last_pps,
       * though we found that this may be unreliable sometimes.
       * */
      double trigger_time;

      /* Trigger information */
      TriggerInfo trigger_info;


    ClassDef(Header,1);
  };



}


#endif
