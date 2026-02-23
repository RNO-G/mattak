#ifndef __MONITORING_H__
#define __MONITORING_H__

#include "TObject.h"
#include <array>


namespace mattak
{

  class EventSummary : public TObject
  {
    public:
      // default constructor
      EventSummary() = default;

      // destructor
      ~EventSummary() = default;

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
      int frun_number = 0;
      int fstation_number = 0;
      int fevent_count = 0;

    ClassDef(RunSummary, 1);

  };


} // namespace mattak

#endif