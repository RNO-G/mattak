#ifndef __MATTAK_CONSTANTS_H__
#define __MATTAK_CONSTANTS_H__


#include <stdint.h> 


namespace mattak 
{
  namespace k
  {
    constexpr uint8_t num_radiant_channels = 24; 
    constexpr uint16_t num_radiant_samples = 2048; 
    constexpr uint16_t num_radiant_thresholds = 24; ///  completely made up 
    constexpr uint16_t num_lt_thresholds = 8; ///  completely made up 
  }
}

#endif
