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
  }
}

#endif
