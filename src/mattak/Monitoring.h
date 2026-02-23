#ifndef __MONITORING_H__
#define __MONITORING_H__

#include "TObject.h"
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

      int rms_per_channel[24]; // root mean square per channel

      ClassDef(EventSummary, 3);
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

    ClassDef(RunSummary, 2);

  };


} // namespace mattak

#endif