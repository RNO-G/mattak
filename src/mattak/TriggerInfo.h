#ifndef __MATTAK_TRIGGER_INFO_H__
#define __MATTAK_TRIGGER_INFO_H__

#include <stdint.h> 
#include "TObject.h" 
#include "mattak/Constants.h" 

namespace mattak 
{

  struct RadiantTriggerInfo 
  {
    uint32_t channel_mask = 0; 
    uint8_t start_windows[mattak::k::num_radiant_channels][2] = {}; 
    uint32_t RF_masks[2] = {}; 
    uint8_t RF_ncoinc[2] = {}; 
    uint8_t RF_window[2] = {}; 

    ClassDef(RadiantTriggerInfo,2); 
  }; 

  struct LTTriggerInfo
  {
    uint8_t window = 0; 
    uint8_t num_coinc = 0; 
    bool vppmode = 0; 
    ClassDef(LTTriggerInfo,1); 
  }; 

  struct TriggerInfo
  {
    bool rf_trigger = false;  // True if this is a type of RF trigger
    bool force_trigger = false;  // True if this is a force trigger (software triggER)
    bool pps_trigger = false;  // True if this is a trigger made by the PPS 
    bool ext_trigger = false;  // True if this is an external trigger
    bool radiant_trigger = false;  // True if this is a trigger from the RADIANT tunnel diodes
    bool lt_trigger = false;  // True if this is a trigger from the low-threshold board
    int8_t which_radiant_trigger = -1; //which of the radiant triggers triggered. This is not reliably known so may be -1 even if radiant_trigger is true; 

    RadiantTriggerInfo radiant_info; 
    LTTriggerInfo lt_info; 
    ClassDef(TriggerInfo,2); 
  }; 


}



#endif
