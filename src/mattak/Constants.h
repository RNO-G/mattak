#ifndef __MATTAK_CONSTANTS_H__
#define __MATTAK_CONSTANTS_H__


#include <stdint.h> 


namespace mattak 
{
  namespace k
  {
    constexpr uint8_t num_radiant_channels = 24; 
    constexpr uint16_t num_radiant_samples = 2048; 
    constexpr uint16_t num_lab4_samples = 4096; 
    constexpr uint16_t num_lt_channels = 4; 
    constexpr uint16_t num_lt_beams = 12;
    constexpr uint16_t radiant_window_size = 128;
    constexpr uint16_t radiant_windows_per_buffer = num_radiant_samples / radiant_window_size;

    /* fallback sampling rate (MHz) for older data that stores it neither with the waveforms nor in the runinfo */
    constexpr float default_radiant_sample_rate = 3200;
  }
}

#endif
