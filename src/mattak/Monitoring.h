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

      uint32_t event_number = 0;
      std::vector<float> rms; // root mean square per channel

      std::vector<float> max_abs_amplitude;

      std::vector<float> glitching_test_statitic;

      std::vector<std::vector<float>> block_offset;

      ClassDef(EventSummary, 1);
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

      // Trigger information
      uint32_t n_events = 0;
      uint32_t n_forced_triggers = 0;
      uint32_t n_lt_triggers = 0;
      uint32_t n_rf0_triggers = 0;
      uint32_t n_rf1_triggers = 0;

      std::vector<float> frequencies; // frequencies of the FFT bins
      std::vector<std::vector<float>> avg_spectrum; // average spectrum per channel (e.g. average FFT amplitude per frequency bin)
      std::vector<std::vector<float>> avg_spectrum_force;
      std::vector<std::vector<float>> avg_spectrum_lt;
      std::vector<std::vector<float>> avg_spectrum_rf0;
      std::vector<std::vector<float>> avg_spectrum_rf1;

      std::vector<uint32_t> glitch_counts;
      std::vector<std::vector<float>> avg_block_offset;

    ClassDef(RunSummary, 1);
  };


} // namespace mattak

#endif