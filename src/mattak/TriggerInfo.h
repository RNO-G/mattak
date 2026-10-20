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

    ClassDef(RadiantTriggerInfo, 2);
  };

  struct LTTriggerInfo
  {
    uint8_t window = 0;
    uint8_t num_coinc = 0;
    bool vppmode = 0;
    uint8_t channel_mask = 0;
    uint16_t beam_mask = 0;
    ClassDef(LTTriggerInfo, 2);
  };

  struct DidaqTriggerInfo
  {
    // DiDAQ replaces RADIANT/FLOWER; didaq_start_offsets is the union alternative to
    // RadiantTriggerInfo::start_windows (see rno_g_header_t in rno-g.h).
    uint16_t start_offsets[mattak::k::num_radiant_channels] = {};

    uint16_t type = 0;
    uint32_t channel_mask = 0;
    uint16_t beam_mask = 0;

    ClassDef(DidaqTriggerInfo, 1);
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
    int8_t which_lt_trigger = -1; //which lt trigger is used to make lt_trigger. 0 is the hi-lo and 1 is the phased. -1 for an error or not lt_trigger

    bool didaq_trigger = false;  // True if this is an RF trigger from the DiDAQ digitizer

    RadiantTriggerInfo radiant_info;
    LTTriggerInfo lt_info;
    DidaqTriggerInfo didaq_info;

    ClassDef(TriggerInfo, 4);
  };


}



#endif
