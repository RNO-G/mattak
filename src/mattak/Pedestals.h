#ifndef __MATTAK_PEDESTALS_H__
#define __MATTAK_PEDESTALS_H__

#include <stdint.h> 
#include <vector> 
#include <array> 
#include "TObject.h" 
class TH1; 

#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h" 
#else
typedef int rno_g_pedestal_t; 
#endif

#include "mattak/Constants.h"

namespace mattak
{
  /** 
   * This stores the pedestals for each channel, along with the time they were made and the vbias   *
   */
  class Pedestals : public TObject
  {

    public: 
      Pedestals() { ; }
      Pedestals(const rno_g_pedestal_t * peds); 
      Pedestals(const char * pedfile); 
      uint32_t when=0; 
      uint32_t nevents=0; 
      uint32_t mask=0; 
      uint8_t flags=0; 
      uint8_t station_number=0; 
      float vbias[2];  //V, -1 means unknown
      uint16_t pedestals[mattak::k::num_radiant_channels][mattak::k::num_lab4_samples] = {}; 
      //Makes a pedestal histogram of the channel, user must delete 
      TH1 * pedestalHist(int chan, const char * name = NULL) const; 
    ClassDef(Pedestals,3); 
    private: 
      void doInit(const rno_g_pedestal_t *peds); 

  }; 


}






#endif
