#ifndef _MATTAK_DETECTOR_H
#define _MATTAK_DETECTOR_H

#include "mattak/Constants.h"
#include <array>
#include <vector>
#include "TTimeStamp.h"
#include "TVector3.h"

namespace mattak
{

  class Detector : public TObject
  {
    public:
      Detector() {;}
      std::array<char, mattak::k::num_radiant_channels> ant_type;
      std::array<double, mattak::k::num_radiant_channels> cable_delays;
      std::array<TVector3, mattak::k::num_radiant_channels> antenna_positions;
      std::array<TVector3, mattak::k::num_radiant_channels> global_antenna_positions;
      TVector3 station_position;

      int station =-1;

      int writeJSON(const char * filename);
      static const Detector* fromDB(int station, const TTimeStamp & t = TTimeStamp(), const char* connection_string = nullptr);
//      static const Detector* fromJSONStream(std::istream & istr);
//      static const Detector* fromJSONFile(const char * filename);


      TTimeStamp det_time;

    private:
      ClassDef(Detector,1);
  };

}


#endif
