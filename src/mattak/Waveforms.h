#ifndef __MATTAK_WAVEFORMS_H__
#define __MATTAK_WAVEFORMS_H__

#include <stdint.h> 
#include <vector> 
#include <array> 
#include "TObject.h" 

#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h" 
#else
typedef int rno_g_waveform_t; 
#endif


#include "mattak/Constants.h"

namespace mattak
{
  /** 
   * This stores "raw" (unrolled, pedestal subtracted) data. The voltage calibration
   * has not been applied yet. The waveforms are stored uncalibrated because they take up
   * less space this way, but can be calibrated on the fly eventually. 
   *
   */
  class Waveforms : public TObject
  {

    public: 
      Waveforms() { ; }
      Waveforms(const rno_g_waveform_t * wf); 
      uint32_t run_number = 0; 
      uint32_t event_number = 0; 
      uint16_t station_number = 0; 
      uint16_t buffer_length = 0; 

      //These are pedestal subtracted, so signed
      int16_t radiant_data[mattak::k::num_radiant_channels][mattak::k::num_radiant_samples] = {}; 
    ClassDef(Waveforms,2); 

  }; 


}

#endif
