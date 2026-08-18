#ifndef __MONITORING_H__
#define __MONITORING_H__

#include "TObject.h"

#include <cstdint>
#include <vector>

namespace mattak
{

  class EventSummary : public TObject
  {
    public:
      // default constructor
      EventSummary() = default;

      // destructor
      ~EventSummary() = default;

      // Basic identifying information
      uint32_t event_number = 0;

      // Per-channel information
      std::vector<float> rms; // root mean square per channel
      std::vector<uint16_t> max_abs_amplitude;
      // largest peak-to-peak amplitude found in a sliding 10 ns window (see mon_write.py)
      std::vector<uint16_t> max_peak_to_peak_amplitude;
      std::vector<float> glitching_test_statitic;
      // we decided to only store the max. abs. block offset per waveform/channel to save space
      std::vector<uint16_t> block_offset;

      ClassDef(EventSummary, 2);
  };


  class RunSummary : public TObject
  {
    public:
      // default constructor
      RunSummary() = default;

      // destructor
      ~RunSummary() = default;

      // Basic identifying information
      uint32_t run_number = 0;
      uint8_t station_number = 0;

      // Trigger information. A run is taken with a single digitizer, so only the RADIANT/LT
      // counters or only the DIDAQ counters are non-zero (see mattak::Dataset::digitizer).
      uint32_t n_events = 0;
      uint32_t n_forced_triggers = 0;
      uint32_t n_lt_triggers = 0;
      uint32_t n_rf0_triggers = 0;
      uint32_t n_rf1_triggers = 0;
      uint32_t n_didaq_coinc0_triggers = 0;
      uint32_t n_didaq_coinc1_triggers = 0;
      uint32_t n_didaq_surf_up_triggers = 0;
      uint32_t n_didaq_surf_down_triggers = 0;
      uint32_t n_didaq_deep_phased_triggers = 0;

      // Per-channel average spectrum information. A per-trigger-type spectrum is left empty if
      // the run contains no event of that trigger type (in particular, all the DIDAQ ones are
      // empty for a RADIANT run and vice versa).
      std::vector<float> frequencies; // frequencies of the FFT bins
      std::vector<std::vector<float>> avg_spectrum; // average spectrum per channel (e.g. average FFT amplitude per frequency bin)
      std::vector<std::vector<float>> avg_spectrum_force;
      std::vector<std::vector<float>> avg_spectrum_lt;
      std::vector<std::vector<float>> avg_spectrum_rf0;
      std::vector<std::vector<float>> avg_spectrum_rf1;
      std::vector<std::vector<float>> avg_spectrum_didaq_coinc0;
      std::vector<std::vector<float>> avg_spectrum_didaq_coinc1;
      std::vector<std::vector<float>> avg_spectrum_didaq_surf_up;
      std::vector<std::vector<float>> avg_spectrum_didaq_surf_down;
      std::vector<std::vector<float>> avg_spectrum_didaq_deep_phased;

    ClassDef(RunSummary, 2);
  };


} // namespace mattak

#endif
