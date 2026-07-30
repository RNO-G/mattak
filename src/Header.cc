#include "mattak/Header.h"
#include "TError.h"
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
  this->trigger_info.rf_trigger = !(this->trigger_info.force_trigger || this->trigger_info.pps_trigger);

  // The DIDAQ bit is the top (last) bit of the now 16-bit trigger_type. When set, the low
  // bits are DIDAQ's own trigger sources, which reuse some RADIANT/FLOWER bit positions for an
  // unrelated meaning (e.g. bit 2 is RF_LT_SIMPLE on RADIANT/FLOWER but RF_DIDAQ_COINC0 on DIDAQ),
  // so they must not be interpreted with the RADIANT/FLOWER logic below.
  bool is_didaq = head->trigger_type & RNO_G_TRIGGER_DIDAQ;

  if (!is_didaq)
  {
    this->trigger_info.radiant_trigger = !!(head->trigger_type & (RNO_G_TRIGGER_RF_RADIANTX));
    this->trigger_info.lt_trigger = !!(head->trigger_type & (RNO_G_TRIGGER_RF_LT_SIMPLE | RNO_G_TRIGGER_RF_LT_PHASED));
    this->trigger_info.didaq_trigger = false;
    if (this->trigger_info.radiant_trigger)
    {
      this->trigger_info.which_radiant_trigger =
        (head->trigger_type & RNO_G_TRIGGER_RF_RADIANT0) ? 0 :
        (head->trigger_type & RNO_G_TRIGGER_RF_RADIANT1) ? 1 :
        -1;
    }
    else this->trigger_info.which_radiant_trigger=-128;

    if (this->trigger_info.lt_trigger)
    {
      this->trigger_info.which_lt_trigger =
        (head->trigger_type & RNO_G_TRIGGER_RF_LT_SIMPLE) ? 0 :
        (head->trigger_type & RNO_G_TRIGGER_RF_LT_PHASED) ? 1 :
        -1;
    }
    else this->trigger_info.which_lt_trigger = -1;

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
    this->trigger_info.lt_info.channel_mask = head->lt_simple_trigger_cfg.channel_mask;
    this->trigger_info.lt_info.beam_mask = head->lt_phased_trigger_cfg.beam_mask;

  }
  else
  {
    // No RADIANT/FLOWER boards on DIDAQ, so radiant_trigger/lt_trigger don't apply.
    this->trigger_info.radiant_trigger = false;
    this->trigger_info.lt_trigger = false;

    this->trigger_info.didaq_info.type = head->trigger_type;

    // RNO_G_TRIGGER_RF_DIDAQ_* are (RNO_G_TRIGGER_DIDAQ | some-bit) combinations, and we already
    // know the DIDAQ bit is set here -- so testing `trigger_type & RNO_G_TRIGGER_RF_DIDAQ_X` with
    // a plain truthiness check would always be true (satisfied by the shared DIDAQ bit alone),
    // regardless of whether bit X itself is set. Compare against the whole combined constant with
    // `==` instead, so both the DIDAQ bit and the specific source bit are required.
    bool didaq_coinc0 = (head->trigger_type & RNO_G_TRIGGER_RF_DIDAQ_COINC0) == RNO_G_TRIGGER_RF_DIDAQ_COINC0;
    bool didaq_coinc1 = (head->trigger_type & RNO_G_TRIGGER_RF_DIDAQ_COINC1) == RNO_G_TRIGGER_RF_DIDAQ_COINC1;
    bool didaq_surf_up = (head->trigger_type & RNO_G_TRIGGER_RF_DIDAQ_SURF_UP) == RNO_G_TRIGGER_RF_DIDAQ_SURF_UP;
    bool didaq_surf_down = (head->trigger_type & RNO_G_TRIGGER_RF_DIDAQ_SURF_DOWN) == RNO_G_TRIGGER_RF_DIDAQ_SURF_DOWN;
    bool didaq_deep_phased = (head->trigger_type & RNO_G_TRIGGER_RF_DIDAQ_DEEP_PHASED) == RNO_G_TRIGGER_RF_DIDAQ_DEEP_PHASED;

    this->trigger_info.didaq_trigger = didaq_coinc0 || didaq_coinc1 || didaq_surf_up || didaq_surf_down || didaq_deep_phased;

    // trigger_mask: "Which channels (or beams?) caused the trigger" -- DIDAQ has no separate
    // channel-mask/beam-mask config fields of its own, so use the per-event trigger_mask for both,
    // filed under whichever of the two DEEP_PHASED (beam-like) vs. everything else (channel-like) applies.
    if (didaq_deep_phased)
      this->trigger_info.didaq_info.beam_mask = head->trigger_mask;
    else if (this->trigger_info.didaq_trigger)
      this->trigger_info.didaq_info.channel_mask = head->trigger_mask;

    // `didaq_start_offsets` is the union alternative to `radiant_start_windows` (see comment above).
    for (int i = 0 ; i < RNO_G_NUM_RADIANT_CHANNELS; i++)
    {
      this->trigger_info.didaq_info.start_offsets[i] = head->didaq_start_offsets[i];
    }
  }
#else
  ::Error("mattak::Header::Header", "Not compiled with librno-g support");
  (void) head;
#endif

}
