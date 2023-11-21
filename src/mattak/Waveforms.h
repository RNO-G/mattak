#ifndef __MATTAK_WAVEFORMS_H__
#define __MATTAK_WAVEFORMS_H__

#include <stdint.h> 
#include <vector> 
#include <array> 
#include "TObject.h" 
class TVirtualPad; 
class TGraph; 

#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h" 
#else
typedef int rno_g_waveform_t; 
#endif


#include "mattak/Constants.h"
#include "mattak/Header.h" 
#include "mattak/VoltageCalibration.h" 

namespace mattak
{


  struct WaveformPlotOptions
  {

    bool autoscale = true; 
    bool global_scale = true; //only relevant if autoscaling 
    double scale_fact = 1.1; 
    bool symmetric_scale = true; 

    //only used when not autoscaling 
    float min = -200;
    float max = 200; 
    bool stats = true; 
    uint32_t mask = 0xffffff; 

    //used if not using an existing pad 
    int width=1800; 
    int height = 900; 

    //time base
    bool ns = true; 
  }; 


  class IWaveforms : public TObject
  {
    public: 
     uint32_t run_number = 0; 
     uint32_t event_number = 0; 
     uint16_t station_number = 0; 
     uint16_t buffer_length = 0; 
     virtual TGraph * makeGraph(int chan, bool ns = true) const = 0; 
     virtual TVirtualPad* drawWaveforms(const WaveformPlotOptions & opt = WaveformPlotOptions(), TVirtualPad * where = nullptr) const = 0; 
     ClassDef(IWaveforms, 1); 
  }; 





  /** 
   * This stores "raw" (unrolled, pedestal subtracted) data. The voltage calibration
   * has not been applied yet. The waveforms are stored uncalibrated because they take up
   * less space this way, but can be calibrated on the fly eventually. 
   *
   */
  class Waveforms : public IWaveforms
  {

    public: 
      Waveforms() { ; }
      Waveforms(const rno_g_waveform_t * wf); 
          //These are pedestal subtracted, so signed

       int16_t radiant_data[mattak::k::num_radiant_channels][mattak::k::num_radiant_samples] = {}; 
       virtual TGraph * makeGraph(int chan, bool ns = true) const; 
       virtual TVirtualPad* drawWaveforms(const WaveformPlotOptions & opt = WaveformPlotOptions(), TVirtualPad * where = nullptr) const; 
    ClassDef(Waveforms,3); 

  }; 


  /** This stores calibrated data with the voltage calibration applied. 
   *
   * Note
   * that we could inherit from Waveforms but I'm not sure that makes sense
   * since then we'd also have the raw data (which we don't need here?) We
   * could define an interface that we both inherit from, but that might
   * require custom streamers for backwards compatibility, I'm not sure. It
   * might be worth doing that anyway, but I won't do that yet.  So instead,
   * we'll duplicate some things, annoyingly. 
   *
   **/ 
  
  class CalibratedWaveforms : public IWaveforms {

    public: 

      CalibratedWaveforms() { ; }
      CalibratedWaveforms(const Waveforms & wf, const Header & h,  const VoltageCalibration & vc, bool isOldFirmware = false); 

      virtual TGraph * makeGraph(int chan, bool ns = true) const; 
      virtual TVirtualPad* drawWaveforms(const WaveformPlotOptions & opt = WaveformPlotOptions(), TVirtualPad * where = nullptr) const; 
      double radiant_data[mattak::k::num_radiant_channels][mattak::k::num_radiant_samples] = {}; 
      ClassDef(CalibratedWaveforms, 1); 
  }; 


}

#endif
