#ifndef __MONITORING_H__
#define __MONITORING_H__

#include "TObject.h"

#include <cstdint>
#include <array>
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
      std::array<float, 24> rms_per_channel; // root mean square per channel

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
      std::array<std::vector<float>, 24> avg_spectrum_per_channel; // average spectrum per channel (e.g. average FFT amplitude per frequency bin)
      std::array<std::vector<float>, 24> avg_spectrum_per_channel_force;
      std::array<std::vector<float>, 24> avg_spectrum_per_channel_lt;
      std::array<std::vector<float>, 24> avg_spectrum_per_channel_rf0;
      std::array<std::vector<float>, 24> avg_spectrum_per_channel_rf1;

    ClassDef(RunSummary, 1);
  };


} // namespace mattak

#endif