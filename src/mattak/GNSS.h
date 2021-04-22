#ifndef __MATTAK_GNSS_H__
#define __MATTAK_GNSS_H__

#include <stdint.h> 
#include "TObject.h" 
#include <vector> 



namespace mattak 
{


  class GNSS
  {

    public:
      GNSS() { ; } 

      class Sat
      {
        public: 

          Sat() { ; } 
          enum SatType
          {
            SatType_NONE, 
            SatType_GPS, 
            SatType_GLONASS, 
            SatType_GALILEO, 
            SatType_BEIDOU 
          } type  = SatType_NONE; 

          int id = -1; 
          bool have_ephemeris = false; 
          bool used_in_fix = false; 
          bool have_almanac = false; 
          float azimuth = 0;
          float elevation = 0; 
          float c_n0 = 0; //C/N0 
          double pseudorange = 0;  
          ClassDef(GNSS::Sat,1); 
      }; 

      uint32_t time;  

      int nsats = 0; 
      std::vector<Sat> sats; 

      float latitude = 0;
      float longitude = 0; 
      float altitude = 0; 
      float accuracy = -1; 


      double pdop = 0; 
      double tdop = 0; 
      double dop = 0; 


    ClassDef(GNSS,1); 
  }; 


}

#endif

