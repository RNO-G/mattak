#ifndef __MATTAK_CONSTANTS_H__
#define __MATTAK_CONSTANTS_H__


#include <stdint.h> 


namespace mattak 
{
  namespace k
  {
    constexpr uint8_t num_radiant_channels = 24;
    constexpr uint16_t num_radiant_samples = 2048;
    constexpr uint16_t num_didaq_samples = 4096; // must match RNO_G_MAX_DIDAQ_NSAMPLES in rno-g.h
    constexpr uint16_t num_lab4_samples = 4096;
    constexpr uint16_t num_lt_channels = 4;
    constexpr uint16_t num_lt_beams = 12;
    constexpr uint16_t num_didaq_coinc = 2; // must match RNO_G_NUM_DIDAQ_COINC in rno-g.h
    constexpr uint16_t num_didaq_beams = 10; // must match RNO_G_NUM_DIDAQ_BEAMS in rno-g.h
    constexpr uint16_t radiant_window_size = 128; 
    constexpr uint16_t radiant_windows_per_buffer = num_radiant_samples / radiant_window_size; 
  }
}

#endif
