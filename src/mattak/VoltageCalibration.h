#ifndef _MATTAK_VOLTAGE_CALIBRAION_H
#define _MATTAK_VOLTAGE_CALIBRAION_H


#ifndef MATTAK_NOROOT
#include "TH2.h"
#include "TTree.h"
#include "TGraph.h"
#endif


#include "mattak/Constants.h"
#include <vector>
#include <array>
namespace mattak
{

  constexpr int max_voltage_calibration_fit_order = 9;

  // Free apply voltage calibration so we can call it without ROOT
  // @param N number of samples (must be multiple of window size!)
  // @param in the input (raw samples, pedestal subtracted)
  //
  double * applyVoltageCalibration(int N, const int16_t * in, double * out, int start_window, bool isOldFirmware, int fit_order,
                        int nResidPoints, const double * packed_fit_params, bool isUsingResid, const double * packed_aveResid_volt, const double * packed_aveResid_adc);


#ifndef MATTAK_NOROOT
  class VoltageCalibration : public TObject
  {
    public:
      VoltageCalibration()  { ; }

      /**
       * Create a Voltage Calibration from bias_scan_file. If bias_scan_file ends with .root, it's assumed to contain a tree named pedestals with branch pedestals.
       * Otherwise it's assumed to be afile which consists of a bunch of raw pedestals with associated vbiases.
       * This does a polynomial fit of order fit_order . for each sample in each channel for voltages betwee fit_min_V and fit_max_V.
       *
       * If fit_Vref is non-zero then the fit is fixed to cross the scan at that vref (useful if that's your nominal pedestal!)
       *
       *
       */
      VoltageCalibration(const char * bias_scan_file, double fit_Vref = 1.5, int fit_order = 9, double fit_min_V  = 0.2, double fit_max_V = 2.2, bool isUsingResid = true);
      VoltageCalibration(TTree * bias_scan_tree, const char * branch_name = "pedestals",  double fit_Vref = 1.5, int fit_order = 9, double fit_min_V  = 0.2, double fit_max_V = 2.2, bool isUsingResid = true);
      void recalculateFits(int fit_order, double fit_min_V, double fit_max_V, double fit_Vref = 1.5, bool isUsingResid = true, uint32_t mask = 0xffffff, int turnover_threshold = 20);
      void saveFitCoeffsInFile();
      void readFitCoeffsFromFile(const char * inFile);

      int getNresidPoints(int chan) const { return nResidPoints[chan>=mattak::k::num_radiant_channels/2]; }
      const double * getPackedAveResid_volt(int chan) const { return &resid_volt[chan>=mattak::k::num_radiant_channels/2][0]; }
      const double * getPackedAveResid_adc(int chan) const { return &resid_adc[chan>=mattak::k::num_radiant_channels/2][0]; }
      int getFitOrder() const { return fit_order; }
      double getFitMin() const { return fit_min; }
      double getFitMax() const { return fit_max; }
      double getFitVref() const { return fit_vref; }
      const double * getFitCoeffs(int chan, int sample) const { return getPackedFitCoeffs(chan) + sample * (getFitOrder()+1); }
      double getFitCoeff(int chan, int sample, int coeff) const { return getFitCoeffs(chan,sample)[coeff]; }
      const double * getPackedFitCoeffs(int chan) const { return &fit_coeffs[chan][0]; }
      double * apply(int chan, int N, const int16_t * in, int start_window, double * out = 0, bool isOldFirmware = false) const
      {
        return applyVoltageCalibration(N, in, out, start_window, isOldFirmware, getFitOrder(), getNresidPoints(chan),
                                       getPackedFitCoeffs(chan), isResid(), getPackedAveResid_volt(chan), getPackedAveResid_adc(chan));
      }
      TH2S * makeHist(int channel) const;
      TGraph * makeAdjustedInverseGraph(int channel, int sample, bool resid=false) const;
      TGraph * makeSampleGraph(int channel, int sample, bool resid=false) const;
      TGraph * getAveResidGraph_dac1() const { return graph_residAve[0]; }
      TGraph * getAveResidGraph_dac2() const { return graph_residAve[1]; }
      TH2S * getResidHist_dac1() const { return hist_resid[0]; }
      TH2S * getResidHist_dac2() const { return hist_resid[1]; }
      int getFitNdof(int channel, int samp) const { return fit_ndof[channel][samp]; }
      double getFitChisq(int channel, int samp) const { return fit_chisq[channel][samp]; }
      double getFitMaxErr(int channel, int samp) const { return fit_maxerr[channel][samp]; }
      int getStationNumber() const { return station_number; }
      uint32_t getStartTime() const { return start_time; }
      uint32_t getEndTime() const { return end_time; }
      int scanSize() const { return vbias[0].size() ; }
      const int16_t * scanADCVals(int channel, int samp) const { return &scan_result[channel][samp][0]; }
      const double * scanBias(int chan) const  { return &vbias[chan>=mattak::k::num_radiant_channels/2][0]; }
      int scanTurnover(int chan, int samp) { return turnover_index[chan][samp]; }
      bool isResid() const { return fit_isUsingResid; }

    private:
      std::array<std::vector<double>, 2> vbias;  //Left, Right
      std::vector<std::array<std::array<int16_t, mattak::k::num_lab4_samples>, mattak::k::num_radiant_channels>> scan_result;
      std::array<std::vector<double>, mattak::k::num_radiant_channels> fit_coeffs;  //packed format, per channel
      std::array<std::array<int, mattak::k::num_lab4_samples>, mattak::k::num_radiant_channels> fit_ndof; //number of degree of freedom
      std::array<std::array<double, mattak::k::num_lab4_samples>, mattak::k::num_radiant_channels> fit_chisq; //sum of difference squared, really...
      std::array<std::array<double, mattak::k::num_lab4_samples>, mattak::k::num_radiant_channels> fit_maxerr; //maximum error
      std::array<std::array<int, mattak::k::num_lab4_samples>, mattak::k::num_radiant_channels> turnover_index; //where we start turning over
      std::array<double, mattak::k::num_radiant_channels> adc_offset;
      std::array<std::array<TGraph*, mattak::k::num_lab4_samples>, mattak::k::num_radiant_channels> graph;
      std::array<std::vector<double>, 2> resid_volt;  // 2 DACs
      std::array<std::vector<double>, 2> resid_adc;  // 2 DACs
      std::array<TGraph*, 2> graph_residAve; // 2 DACs
      std::array<int, 2> nResidPoints; // 2 DACs
      std::array<TH2S*, 2> hist_resid; // 2 DACs
      std::array<bool, mattak::k::num_radiant_channels> isBad_channelAveChisqPerDOF;
      std::array<std::array<bool, mattak::k::num_lab4_samples>, mattak::k::num_radiant_channels> isBad_sampChisqPerDOF;
      std::array<std::array<bool, 4>, 2> isResidOutOfBoxFrame; // 4 thresholds for each DAC
      int fit_order;
      int station_number;
      double fit_vref;
      double fit_min;
      double fit_max;
      uint32_t start_time;
      int turnover_threshold;
      void setupFromTree(TTree*t, const char * branch_name, double vref, int order, double min, double max, bool isUsingResid);
      uint32_t end_time;
      bool hasBiasScanData;
      bool fit_isUsingResid = true;
      bool left_equals_right;
    ClassDef(VoltageCalibration, 2);
  };
#endif



}


#endif
