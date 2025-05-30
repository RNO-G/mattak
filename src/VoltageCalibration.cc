#include "mattak/VoltageCalibration.h"
#include <iostream>
#include <stdio.h>

#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h"
#endif


static double evalPars(double x, int order, const double * p)
{
  // calculate a polynominal of n order backwards
  // a_n = p_n; a_n-1 = a_n * x + p_n-1; ...
  double ans = p[order];
  int i = order - 1;
  while (i >= 0)
  {
    ans = x * ans + p[i--];
  }

  return ans;
}

static double* adcTablePerSample(int order, int npoints, const double * par, const double * resid_volt, const double * resid_adc)
{
  double *adcTable = new double[npoints];

  for (int i = 0; i < npoints; i++)
  {
    // When we perform a calibration with residuals, we fit f(V) -> ADC
    adcTable[i] = evalPars(resid_volt[i], order, par) + resid_adc[i];
  }

  return adcTable;
}

static double adcToVolt(double in_adc, int npoints, const double * volt_array, const double * adc_array)
{
  double m;
  double out_volt = 0;

  // If in_adc is zero, out_volt is zero
  if (in_adc == 0) return out_volt;

  // If in_adc is out of range...
  if (in_adc < adc_array[0])
  {
    m = (volt_array[1] - volt_array[0]) / (adc_array[1] - adc_array[0]);
    out_volt = volt_array[0] + (in_adc - adc_array[0]) * m;
    return out_volt;
  }
  if (in_adc > adc_array[npoints-1])
  {
    m = (volt_array[npoints-1] - volt_array[npoints-2]) / (adc_array[npoints-1] - adc_array[npoints-2]);
    out_volt = volt_array[npoints-1] + (in_adc - adc_array[npoints-1]) * m;
    return out_volt;
  }

  for (int i = 0; i < npoints; i++)
  {
    if (in_adc == adc_array[i])
    {
      out_volt = volt_array[i]; // Lucky if this happens!
      return out_volt;
    }

    if (i < npoints-1)
    {
      // Most likely we will get out_volt from interpolation
      if (in_adc > adc_array[i] && in_adc < adc_array[i+1])
      {
        m = (volt_array[i+1] - volt_array[i]) / (adc_array[i+1] - adc_array[i]);
        out_volt = volt_array[i] + (in_adc - adc_array[i]) * m;
        return out_volt;
      }
    }
  }

  return out_volt;
}


#ifndef MATTAK_NOROOT

#include "mattak/Pedestals.h"

#include "TLinearFitter.h"
#include "TFile.h"
#include "TF1.h"
#include "TList.h"
#include "TString.h"
#include "TBox.h"

ClassImp(mattak::VoltageCalibration);


mattak::VoltageCalibration::VoltageCalibration(TTree * tree, const char * branch_name, double vref, int fit_order, double min, double max, bool isUsingResid)
  : TObject()
{
  setupFromTree(tree, branch_name, vref, fit_order, min, max, isUsingResid);
}

// almost copy-pasted from raw file reader...
void mattak::VoltageCalibration::setupFromTree(TTree * tree, const char * branch_name, double vref, int fit_order, double min, double max, bool isUsingResid)
{

  mattak::Pedestals * ped = 0;
  tree->SetBranchAddress(branch_name,&ped);

  for (int i = 0; i < tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
    double bias_l = ped->vbias[0];
    double bias_r = ped->vbias[1];
    if (!scanSize())
    {
      start_time = ped->when;
      station_number = ped->station_number;
    }
    else
    {
      end_time = ped->when;
    }
    vbias[0].push_back(bias_l);
    vbias[1].push_back(bias_r);
    scan_result.emplace_back();
    memcpy(&scan_result.back()[0][0], ped->pedestals, sizeof(ped->pedestals));
  }
  delete ped;
  ped = 0;
  scan_result.shrink_to_fit();
  hasBiasScanData = true;
  recalculateFits(fit_order, min, max, vref, isUsingResid);
}

mattak::VoltageCalibration::VoltageCalibration(const char * bias_scan_file, double vref, int fit_order, double min, double max, bool isUsingResid)
: TObject()
{

  /*
     Work out the extension of the file bias_scan_file.
     If the extension after the last period ('.') is '.root',
     assume we are loading coefficients or the bias scan file has been loaded in ROOT format.
     Otherwise they must be in a file format that required librno-g to read.
  */

  const char * suffix = strrchr(bias_scan_file,'.');

  if (!strcmp(suffix,".root") && !(*(suffix + sizeof(".root")-1)))
  {
    TFile *f = TFile::Open(bias_scan_file);
    // check to see if the file opened properly
    if (!f)
    {
      std::cerr << "Could not open apparenty ROOT file in " << bias_scan_file << std::endl;
      return;
    }

   //check if we are likely loading coeffients

    if (f->Get("general_tree"))
    {
      if (!readFitCoeffsFromFile(f))
      {
        std::cerr << "File looks like it has fit coeffs, but I problem reading" << std::endl;
      }
      delete f;
      return;
    }


    TTree * t = (TTree*) f->Get("pedestals");
    if (!t)
    {
      std::cerr << "Could not open tree pedestals in " << bias_scan_file << std::endl;
      delete f;
      return;
    }
    hasBiasScanData = true;
    setupFromTree(t,"pedestals", vref, fit_order, min, max, isUsingResid);
    delete f;
    return;
  }
#ifndef LIBRNO_G_SUPPORT

  (void) bias_scan_file;
  throw std::runtime_error(
    "Bias Scan File has not been loaded in ROOT format.\n"
    "In order to read raw bias scan data, librno-g support is required in mattak.\n"
    "But librno-g support is not detected. "
    "Likely LIBRNO_G_SUPPORT = OFF in the mattak CMakeLists.txt."
  );

#else


  rno_g_pedestal_t ped;
  rno_g_file_handle_t h;
  if (rno_g_init_handle(&h, bias_scan_file,"r"))
  {
    std::cerr <<"Trouble opening "<< bias_scan_file << std::endl;
    return;
  }

  while  (rno_g_pedestal_read(h, &ped) )
  {
    double bias_l = ped.vbias[0] / 4095. *3.3;
    double bias_r = ped.vbias[1] / 4095. *3.3;
    if (!scanSize())
    {
      start_time = ped.when;
      station_number= ped.station;
    }
    else
    {
      end_time = ped.when;
    }

    vbias[0].push_back(bias_l);
    vbias[1].push_back(bias_r);
    scan_result.emplace_back();
    memcpy(&scan_result.back()[0][0], ped.pedestals, sizeof(ped.pedestals));
  }
  scan_result.shrink_to_fit();
  hasBiasScanData = true;
  recalculateFits(fit_order, min, max, vref, isUsingResid);
  rno_g_close_handle(&h);
#endif

}

mattak::VoltageCalibration::~VoltageCalibration()
{
  delete graphs;

  for (auto h : hist_resid) delete h;

}
static TString formula[1+mattak::max_voltage_calibration_fit_order] = {"pol0","pol1","pol2","pol3","pol4","pol5","pol6","pol7","pol8","pol9"};


/// THIS IS IN cASE WE WANTED TO SHIFT BACK IN TERMS OF ABSOLUTE ADU. BUT WE DON'T!

//static void taylor_shift(int n, double * a, double x0)
//{
//  if (!x0) return;
//
//  double t[n+1][n+1]; //VLA should be fine here since order is constrained
//  //algorithm taken from here: https://planetcalc.com/7726/
//
//
//  double a_n_x0n = a[n] * pow(x0,n);
//  for (int i = 0; i < n; i++)
//  {
//    t[i][0] = a[n-i-1] * pow(x0,n-i-1);
//    t[i][i+1] = a_n_x0n;
//  }
//
//  for (int j = 0; j < n; j++)
//  {
//    for (int i = j+1; i <=n; i++)
//    {
//
//      t[i][j+1] = t[i-1][j] + t[i-1][j+1];
//    }
//  }
//
//  for (int i = 0; i < n; i++)
//  {
//    a[i] = t[n][i+1] * pow(x0,-i);
//  }
//}


void mattak::VoltageCalibration::recalculateFits(int order, double min, double max, double vref, bool isUsingResid, uint32_t mask, int turnover_threshold)
{

  if (!hasBiasScanData)
  {
    std::cerr << "Cannot recalculate fits without bias scan data " << std::endl;
  }

  if (!graphs) graphs = new std::array<std::array<TGraph, mattak::k::num_lab4_samples>, mattak::k::num_radiant_channels>{};
  gErrorIgnoreLevel = kFatal;


  fit_order = order < max_voltage_calibration_fit_order ? order : max_voltage_calibration_fit_order;
  order = fit_order;

  fit_isUsingResid = isUsingResid;

  if (order < max_voltage_calibration_fit_order)
  {
    std::cout << "\nYou are using " << order << "-degree polynomials for the fitting..." << std::endl;
    std::cout << "SUGGESTION: Using order 9 (default) is highly suggested for getting better calibration results." << std::endl;
  }

  if (!fit_isUsingResid) std::cout << "\nSUGGESTION: Extra term the residual function is NOT USED, you may turn it on to improve voltage calibration." << std::endl;

  fit_min = min;
  this->turnover_threshold = turnover_threshold;
  fit_max = max;
  fit_vref = vref;

  if (min > max)
  {
    fit_min = max;
    fit_max = min;
  }

  int nResidSets[mattak::k::num_radiant_channels];
  std::array<std::vector<double>, mattak::k::num_radiant_channels> residAve_volt;
  std::array<std::vector<double>, mattak::k::num_radiant_channels> residAve_adc;
  std::array<std::vector<double>, mattak::k::num_radiant_channels> residVar_adc;

  TLinearFitter fit(1, formula[order]);
  for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++)
  {
    nResidSets[ichan] = 0;

    fit_coeffs[ichan].clear();
    if ( (mask & (1 << ichan))  == 0) continue;

    //check if channel seems broken by looking for zeroes

    fit_coeffs[ichan].resize((order+1) * mattak::k::num_lab4_samples, 0);

    int nbroken = 0;

    for (int i = 0; i < mattak::k::num_lab4_samples; i++)
    {
      fit.ClearPoints();

      int jmin = 0;
      int jmax = 0;

      if (vref)
      {
        double last_v = 0;
        for (unsigned j = 0; j < scanSize(); j++)
        {
          double v = ichan < mattak::k::num_radiant_channels / 2 ? vbias[0][j]: vbias[1][j];
          if (v >= vref)
          {
            if (v == vref || j == 0)
            {
              adc_offset[ichan] = scan_result[j][ichan][i];
            }
            else //have to interpolate
            {
              double frac_low = (v - vref) / (v-last_v);
              adc_offset[ichan] = scan_result[j][ichan][i] * (1-frac_low) + frac_low * scan_result[j-1][ichan][i];
            }
            break;
          }
          last_v = v;

        }
      }


      double last_adc = -4096;
      int nzero = 0;
      for (unsigned j = 0; j < scanSize(); j++)
      {
        double v = ichan < mattak::k::num_radiant_channels / 2 ? vbias[0][j]: vbias[1][j];
        if (v >= fit_min && v <= fit_max)
        {
          if (vref)
          {
            v-=vref;
          }
          double adc = scan_result[j][ichan][i] - adc_offset[ichan];
          if (adc == 0) nzero++;

          //don't include if it turns over
          if (!jmin ||  adc  > last_adc - turnover_threshold)
          {
            if (fit_isUsingResid) (*graphs)[ichan][i].SetPoint((*graphs)[ichan][i].GetN(),v,adc);
            else (*graphs)[ichan][i].SetPoint((*graphs)[ichan][i].GetN(),adc,v);
            jmax = j;
          }
          last_adc = adc;

          if (!jmin) jmin = j;
        }
      }

      int npoints = (*graphs)[ichan][i].GetN();
      double *data_adc;
      double *data_v;

      if (fit_isUsingResid)
      {
        data_adc = (*graphs)[ichan][i].GetY();
        data_v = (*graphs)[ichan][i].GetX();
      }
      else
      {
        data_adc = (*graphs)[ichan][i].GetX();
        data_v = (*graphs)[ichan][i].GetY();
      }

      for (int i = 0; i < npoints; i++)
      {
        if (fit_isUsingResid) fit.AddPoint(&data_v[i], data_adc[i]);
        else fit.AddPoint(&data_adc[i], data_v[i]);
      }

      if (vref) fit.FixParameter(0, 0);

      if (nzero > 1)
      {
        isSampleBroken[ichan][i] = true;
        nbroken++;
      }
      else
      {
        isSampleBroken[ichan][i] = false;
        fit.Eval();
      }

      if (vref) fit.ReleaseParameter(0);

      for (int iorder = 0; iorder <= fit_order; iorder++)
      {
        fit_coeffs[ichan][i * (order+1) + iorder] = fit.GetParameter(iorder);
      }

/*** leave it in terms of vref ***/

//      if (vref)
//      {
//        //we need to unshift our polynomial.
//        //first let's grab the coefficient
//        double * a = &fit_coeffs[ichan][i*(order+1)];
//
//        //add y offset back
////        a[0] += vref;
//
//        //use taylor shift to do the x offset (Shaw/Traub method)
////        taylor_shift(order, a, -adc_offset);
//      }
//
//
//

      fit_ndof[ichan][i] = npoints;  // number of points including start and end
      fit_ndof[ichan][i] -= (fit_order + 1);   //  subtract number of parameters of the fit function
      if (vref) fit_ndof[ichan][i] += 1;  // if parameter 0 is fixed, add back one degree of freedom

      turnover_index[ichan][i] = jmax+1;

      if (fit_isUsingResid && !isSampleBroken[ichan][i])
      { // Sum up all the ADC(V) residuals

        nResidSets[ichan] ++;

        if (residAve_volt[ichan].size() == 0 || residAve_adc[ichan].size() == 0 || residVar_adc[ichan].size() == 0)
        {
          residAve_volt[ichan].resize(npoints);
          residAve_adc[ichan].resize(npoints);
          residVar_adc[ichan].resize(npoints);
        }

        for (int j = 0 ; j < npoints; j++)
        {
          double adc = data_adc[j];
          double v = data_v[j];
          double adcResid = adc - evalPars(v, fit_order, &fit_coeffs[ichan][i * (order+1)]);
          residAve_volt[ichan][j] += v;
          residAve_adc[ichan][j] += adcResid;
          residVar_adc[ichan][j] += pow(adcResid, 2);
        }
      }

    }
    if (nbroken) printf("WARNING: Channel %d seems to have %d broken samples?\n", ichan, nbroken);
  }


  //
  // Calculate max deviation and chi squared
  //
  std::cout << "Calculating max deviation and chi squared for fit quality validation..." << std::endl;

  std::vector<double> aveChisq(mattak::k::num_radiant_channels);

  TBox *smallBox[mattak::k::num_radiant_channels];
  TBox *bigBox[mattak::k::num_radiant_channels];

  // Constants for the residual histograms
  const int nBinsY = 22;
  const int histLowY = -55;
  const int histHighY = 55;
  int nBinsX;
  double histLowX;
  double histHighX;

  for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++)
  {
    std::cout << "Channel " << ichan;
    if ( (mask & (1 << ichan))  == 0) continue;

    // Calculate the average ADC(V) residuals and make a TGraph for each channel
    TString graphNameTitle = TString::Format("aveResid_s%d_c%d_%d_%d", station_number, ichan, start_time, end_time);
    graph_residAve[ichan] = new TGraphErrors();
    graph_residAve[ichan]->SetNameTitle(graphNameTitle);
    graph_residAve[ichan]->SetTitle(graphNameTitle);
    graph_residAve[ichan]->GetXaxis()->SetTitle("VBias [Volt]");
    graph_residAve[ichan]->GetYaxis()->SetTitle("ADC Residual");

    if (fit_isUsingResid)
    {
      int npoints_residGraph = residAve_volt[ichan].size();
      for (int ipoint = 0; ipoint < npoints_residGraph; ipoint++)
      {
        residAve_volt[ichan][ipoint] /= nResidSets[ichan];
        residAve_adc[ichan][ipoint] /= nResidSets[ichan];
        residVar_adc[ichan][ipoint] = abs(residVar_adc[ichan][ipoint]/nResidSets[ichan] - pow(residAve_adc[ichan][ipoint], 2));
        graph_residAve[ichan]->SetPoint(ipoint, residAve_volt[ichan][ipoint], residAve_adc[ichan][ipoint]);
        graph_residAve[ichan]->SetPointError(ipoint, 0, sqrt(residVar_adc[ichan][ipoint]));
      }

      resid_volt[ichan].resize(npoints_residGraph*2-1);
      resid_adc[ichan].resize(npoints_residGraph*2-1);
      nResidPoints[ichan] = resid_volt[ichan].size();

      // Interpolating the average residuals
      for (int i = 0; i < npoints_residGraph; i++)
      {
        resid_volt[ichan][i*2] = graph_residAve[ichan]->GetPointX(i);
        resid_adc[ichan][i*2] = graph_residAve[ichan]->GetPointY(i);
      }
      for (int i = 0; i < npoints_residGraph-1; i++)
      {
        resid_volt[ichan][i*2+1] = (resid_volt[ichan][i*2] + resid_volt[ichan][i*2+2])/2;
        resid_adc[ichan][i*2+1] = resid_adc[ichan][i*2] + (resid_volt[ichan][i*2+1] - resid_volt[ichan][i*2])*(resid_adc[ichan][i*2+2] - resid_adc[ichan][i*2])/(resid_volt[ichan][i*2+2] - resid_volt[ichan][i*2]);
      }

      // Residual histograms for each channel
      TString histNameTitle = TString::Format("residHist_s%d_c%d_%d_%d", station_number, ichan, start_time, end_time);
      nBinsX = npoints_residGraph;
      histLowX = graph_residAve[ichan]->GetPointX(0);
      histHighX = graph_residAve[ichan]->GetPointX(npoints_residGraph-1);

      hist_resid[ichan] = new TH2S(histNameTitle, histNameTitle, nBinsX, histLowX, histHighX, nBinsY, histLowY, histHighY);
      hist_resid[ichan]->SetDirectory(nullptr);
      hist_resid[ichan]->GetXaxis()->SetTitle("VBias [Volt]");
      hist_resid[ichan]->GetYaxis()->SetTitle("ADC Residual");
    }

    int sampCount = 0;
    aveChisq[ichan] = 0;

    std::vector<int> badFit;

    bool channelHasBadFit = false;
    isBad_channelAveChisqPerDOF[ichan] = false;

    for (int i = 0; i < mattak::k::num_lab4_samples; i++)
    {
      if (!(i%128))
      {
        printf(".");
        fflush(stdout); // Print a dot every 128 samples processed
      }

      fit_chisq[ichan][i] = 0;
      fit_maxerr[ichan][i] = 0;
      isBad_sampChisqPerDOF[ichan][i] = false;

      int npoints = (*graphs)[ichan][i].GetN();
      double *data_adc;
      double *data_v;
      double *adcTable = nullptr;

      if (fit_isUsingResid)
      {
        data_adc = (*graphs)[ichan][i].GetY();
        data_v = (*graphs)[ichan][i].GetX();
        adcTable = adcTablePerSample(fit_order, nResidPoints[ichan], &fit_coeffs[ichan][i * (order+1)], &resid_volt[ichan][0], &resid_adc[ichan][0]);
      }
      else
      {
        data_adc = (*graphs)[ichan][i].GetX();
        data_v = (*graphs)[ichan][i].GetY();
      }

      // Calculate max deviation and chi squared
      for (int j = 0 ; j < npoints; j++)
      {
        double adc = data_adc[j];
        double v_meas = data_v[j];
        double v_pred;

        if (fit_isUsingResid) v_pred = adcToVolt(adc, nResidPoints[ichan], &resid_volt[ichan][0], &adcTable[0]);
        else v_pred = evalPars(adc, fit_order, &fit_coeffs[ichan][i * (order+1)]);

        double delta = fabs(v_meas-v_pred);
        if (delta > fit_maxerr[ichan][i]) fit_maxerr[ichan][i] = delta;
        fit_chisq[ichan][i] += ( delta*delta/(0.002*0.002) );

        if (fit_isUsingResid)
        {
          double histX = v_meas;
          double histY = adc - (evalPars(v_meas, fit_order, &fit_coeffs[ichan][i * (order+1)]) + graph_residAve[ichan]->GetPointY(j));
          hist_resid[ichan]->Fill(histX, histY);
        }
      }

      if (fit_chisq[ichan][i]/fit_ndof[ichan][i] > 30.0)
      {
        channelHasBadFit = true;
        isBad_sampChisqPerDOF[ichan][i] = true;
        badFit.push_back(i);
      }

      if (!isSampleBroken[ichan][i])
      {
        aveChisq[ichan] = aveChisq[ichan] + (fit_chisq[ichan][i]/fit_ndof[ichan][i]);
        sampCount++;
      }

      delete [] adcTable;
    }

    // chi2 check for fit quality validation
    aveChisq[ichan] /= sampCount;
    if (aveChisq[ichan] > 6.0)
    {
      isBad_channelAveChisqPerDOF[ichan] = true;
      printf("\nBAD FITTING WARNING: The average chi2/DOF over all samples of CH%d is %f (> 6.0)!!!", ichan, aveChisq[ichan]);
    }

    if (aveChisq[ichan] <= 6.0 && channelHasBadFit)
    {
      for (int samp = 0; samp < badFit.size(); samp++)
      {
        int bad = badFit[samp];
        printf("\nBAD FITTING WARNING: chi2/DOF of sample %d in CH%d is %f (> 30.0)!!!", bad, ichan, fit_chisq[ichan][bad]/fit_ndof[ichan][bad]);
      }
    }

    // Residual histograms and the box frame check for fit quality validation
    if (fit_isUsingResid)
    {
      for (int j = 0; j < 4; j++)
      {
        isResidOutOfBoxFrame[ichan][j] = false;
      }

      nBinsX = hist_resid[ichan]->GetNbinsX();
      histLowX = hist_resid[ichan]->GetXaxis()->GetXmin();
      histHighX = hist_resid[ichan]->GetXaxis()->GetXmax();

      int smallBoxBinX1 = (int) nBinsX*0.55;
      int smallBoxBinX2 = (int) nBinsX*0.75;
      int smallBoxBinY1 = (int) nBinsY*0.25 + 1;
      int smallBoxBinY2 = (int) nBinsY*0.75 + 1;

      double smallBoxX1 = histLowX + (histHighX - histLowX)*smallBoxBinX1/nBinsX;
      double smallBoxX2 = histLowX + (histHighX - histLowX)*smallBoxBinX2/nBinsX;
      double smallBoxY1 = histLowY + (histHighY - histLowY)*smallBoxBinY1/nBinsY;
      double smallBoxY2 = histLowY + (histHighY - histLowY)*(smallBoxBinY2-1)/nBinsY;

      double bigBoxX1 = histLowX;
      double bigBoxX2 = histHighX;
      double bigBoxY1 = histLowY + (histHighY - histLowY)/nBinsY;
      double bigBoxY2 = histHighY - (histHighY - histLowY)/nBinsY;

      smallBox[ichan] = new TBox(smallBoxX1, smallBoxY1, smallBoxX2, smallBoxY2);
      smallBox[ichan]->SetFillStyle(0);
      smallBox[ichan]->SetLineStyle(9);
      smallBox[ichan]->SetLineWidth(3);
      smallBox[ichan]->SetLineColor(2);
      bigBox[ichan] = new TBox(bigBoxX1, bigBoxY1, bigBoxX2, bigBoxY2);
      bigBox[ichan]->SetFillStyle(0);
      bigBox[ichan]->SetLineStyle(9);
      bigBox[ichan]->SetLineWidth(3);
      bigBox[ichan]->SetLineColor(2);

      bool aboveSmallBoxY2 = false;
      bool belowSmallBoxY1 = false;
      bool aboveBigBoxY2 = false;
      bool belowBigBoxY1 = false;

      for (int binNumberX = 1; binNumberX <= graph_residAve[ichan]->GetN(); binNumberX++)
      {
        // Small Box Check
        if (binNumberX > smallBoxBinX1 && binNumberX <= smallBoxBinX2)
        {
          if (hist_resid[ichan]->GetBinContent(binNumberX, smallBoxBinY2) > 1 && !aboveSmallBoxY2)
          {
            aboveSmallBoxY2 = true;
            isResidOutOfBoxFrame[ichan][0] = true;
            printf("\nBAD FITTING WARNING: Some residuals in CH%d go beyond the small box upper threshold (> 25 adu)!!!", ichan);
          }
          if (hist_resid[ichan]->GetBinContent(binNumberX, smallBoxBinY1) > 1 && !belowSmallBoxY1)
          {
            belowSmallBoxY1 = true;
            isResidOutOfBoxFrame[ichan][1] = true;
            printf("\nBAD FITTING WARNING: Some residuals in CH%d go below the small box lower threshold (< -25 adu)!!!", ichan);
          }
        }

        // Big Box Check
        if (hist_resid[ichan]->GetBinContent(binNumberX, nBinsY) > 1 && !aboveBigBoxY2)
        {
          aboveBigBoxY2 = true;
          isResidOutOfBoxFrame[ichan][2] = true;
          printf("\nBAD FITTING WARNING: Some residuals in CH%d go beyond the big box upper threshold (> 50 adu)!!!", ichan);
        }
        if (hist_resid[ichan]->GetBinContent(binNumberX, 1) > 1 && !belowBigBoxY1)
        {
          belowBigBoxY1 = true;
          isResidOutOfBoxFrame[ichan][3] = true;
          printf("\nBAD FITTING WARNING: Some residuals in CH%d go below the big box lower threshold (< -50 adu)!!!", ichan);
        }

        if (aboveBigBoxY2 && belowBigBoxY1 && aboveSmallBoxY2 && belowSmallBoxY1) break;
      }

      hist_resid[ichan]->GetListOfFunctions()->Add(bigBox[ichan]);
      hist_resid[ichan]->GetListOfFunctions()->Add(smallBox[ichan]);
    }

    std::cout << std::endl;
  }

}


TH2S * mattak::VoltageCalibration::makeHist(int chan) const
{
  if (!hasBiasScanData) { printf("\nWARNING: Need to get data from a bias scan file in order to plot histograms!\n"); return 0; }

  int nV = scanSize();
  if (nV < 2) return 0; //makes no sense!

  bool right =  chan >= mattak::k::num_radiant_channels / 2  ;
  double min_V =  vbias[right][0]- fit_vref;
  double dmin_V = vbias[right][1]- min_V - fit_vref;
  double max_V = vbias[right][nV-1] - fit_vref;
  double dmax_V =  max_V - (vbias[right][nV-2] - fit_vref);

  //set up the bin edges
  std::vector<double> Vs;
  Vs.reserve(nV+1);
  Vs.push_back(min_V-dmin_V/2);
  for (int i = 0; i < nV-1; i++)
  {
    Vs.push_back( (vbias[right][i]+ vbias[right][i+1])/2 - fit_vref);
  }
  Vs.push_back(max_V + dmax_V/2);


  TH2S * h = new TH2S(Form("hbias_s%d_c%d_%d_%d", station_number, chan, start_time,end_time),
                      Form("Bias Scan, Station %d, Channel %d, Time [%d-%d], VRef=%g ; Sample ; Vbias [Volt] ; adu", station_number, chan, start_time, end_time,fit_vref),
                      mattak::k::num_lab4_samples, 0, mattak::k::num_lab4_samples,
                      nV, min_V, max_V
                      );
  h->GetYaxis()->Set(nV,  &Vs[0]);


  for (int j = 1; j <= nV; j++)
  {
    for (int i = 1; i <=  mattak::k::num_lab4_samples; i++)
    {
      h->SetBinContent(i,j, scan_result[j-1][chan][i-1] - adc_offset[chan]);
    }
  }

  h->SetStats(0);
  h->SetDirectory(0);

  return h;
}

TGraph * mattak::VoltageCalibration::makeAdjustedInverseGraph(int chan, int samp, bool resid) const
{
  if (!graphs)
  {
    std::cerr << "Cannot use makeAdjustedInverseGraph from  saved coefficients" << std::endl;
    return nullptr;
  }

  if (!hasBiasScanData)
  {
    printf("\nWARNING: Need to get data from a bias scan file in order to make graphs!\n");
    return nullptr;
  }

  if (!fit_isUsingResid)
  {
    printf("\nWARNING: Plots can only be made with function 'makeAdjustedInverseGraph()' when 'fit_isUsingResid' is TRUE!\n");
    return nullptr;
  }

  TGraph *g = new TGraph();
  g->SetName(Form("gsample_inverse_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time));
  g->SetTitle(Form("Station %d Ch %d sample %d [%d-%d]   %s", station_number, chan, samp, start_time, end_time, resid ? "(residuals)" : ""));
  g->GetXaxis()->SetTitle("VBias [Volt]");
  g->GetYaxis()->SetTitle(resid ? "ADC Residual" : "ADC");

  int npoints = (*graphs)[chan][samp].GetN();
  double *data_adc = (*graphs)[chan][samp].GetY();
  double *data_v = (*graphs)[chan][samp].GetX();

  TF1 * fn = new TF1(Form("fsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time), formula[fit_order], fit_min, fit_max, TF1::EAddToList::kNo);
  fn->SetParameters(getFitCoeffs(chan,samp));
  fn->SetRange(data_v[0],data_v[npoints-1]);
  fn->SetLineColor(2);

  for (int j = 0; j < npoints; j++)
  {
    double v = data_v[j];
    double adc = data_adc[j] - graph_residAve[chan]->GetPointY(j);

    if (resid) adc -= fn->Eval(v);
    g->SetPoint(g->GetN(),v,adc);
  }

  if (!resid)
  {
    g->GetListOfFunctions()->Add(fn);
    fn->SetParent(g);
    fn->Save(data_v[0],data_v[npoints-1],0,0,0,0);
  }

  return g;
}

TGraph * mattak::VoltageCalibration::makeSampleGraph(int chan, int samp, bool resid) const
{

  if (!graphs)
  {
    std::cerr << "Cannot use makeSampleGraph from  saved coefficients" << std::endl;
    return nullptr;
  }

  if (!hasBiasScanData)
  {
    printf("\nWARNING: Need to get data from a bias scan file in order to make graphs!\n");
    return nullptr;
  }

  TGraph *g = new TGraph();
  g->SetName(Form("gsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time));
  g->SetTitle(Form("Station %d Ch %d sample %d [%d-%d], #chi^{2}= %g   %s", station_number, chan, samp, start_time, end_time, fit_chisq[chan][samp], resid ? "(residuals)" : ""));
  g->GetXaxis()->SetTitle("ADC");
  g->GetYaxis()->SetTitle(resid ? "(VBias - Predicted VBias) [Volt]" : "VBias [Volt]");

  int npoints = (*graphs)[chan][samp].GetN();
  double *data_adc;
  double *data_v;
  double *adcTable = nullptr;

  TF1 *fn;

  if (fit_isUsingResid)
  {
    data_adc = (*graphs)[chan][samp].GetY();
    data_v = (*graphs)[chan][samp].GetX();
    adcTable = adcTablePerSample(fit_order, nResidPoints[chan], getFitCoeffs(chan,samp), getPackedAveResid_volt(chan), getPackedAveResid_adc(chan));
  }
  else
  {
    data_adc = (*graphs)[chan][samp].GetX();
    data_v = (*graphs)[chan][samp].GetY();

    fn = new TF1(Form("fsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time), formula[fit_order], fit_min, fit_max, TF1::EAddToList::kNo);
    fn->SetParameters(getFitCoeffs(chan,samp));
    fn->SetRange(data_adc[0],data_adc[npoints-1]);
    fn->SetLineColor(2);
  }

  for (int j = 0; j < npoints; j++)
  {
    double v = data_v[j];
    double adc = data_adc[j];

    if (resid)
    {
      if (fit_isUsingResid) v -= adcToVolt(adc, nResidPoints[chan], getPackedAveResid_volt(chan), adcTable);
      else v -= fn->Eval(adc);
    }

    g->SetPoint(g->GetN(),adc,v);
  }

  if (!resid && !fit_isUsingResid)
  {
    g->GetListOfFunctions()->Add(fn);
    fn->SetParent(g);
    fn->Save(data_adc[0],data_adc[npoints-1],0,0,0,0);
  }

  delete [] adcTable;

  return g;
}

TH2S * mattak::VoltageCalibration::getResidHist(int chan) const
{
  TH2S * h;

  if (!fit_isUsingResid)
  {
    printf("\nWARNING: Plots can only be made with function 'getResidHist()' when 'fit_isUsingResid' is TRUE!\n");
    return nullptr;
  }
  else
  {
    h = hist_resid[chan];
  }

  return h;
}

TGraphErrors * mattak::VoltageCalibration::getAveResidGraph(int chan) const
{
  TGraphErrors * g;

  if (!fit_isUsingResid)
  {
    printf("\nWARNING: Plots can only be made with function 'getAveResidGraph()' when 'fit_isUsingResid' is TRUE!\n");
    return nullptr;
  }
  else
  {
    g = graph_residAve[chan];
  }

  return g;
}

void mattak::VoltageCalibration::saveFitCoeffsInFile()
{
  TString outFileName;
  if (fit_isUsingResid) outFileName = TString::Format("volCalConsts_pol%d_s%d_%d-%d.root", getFitOrder(), getStationNumber(), getStartTime(), getEndTime());
  else outFileName = TString::Format("volCalConsts_pol%d_noResid_s%d_%d-%d.root", getFitOrder(), getStationNumber(), getStartTime(), getEndTime());
  TFile f(outFileName, "RECREATE", "", ROOT::RCompressionSetting::EDefaults::kUseSmallest);

  TTree *general_tree = new TTree("general_tree", "general_tree");
  general_tree->Branch("fitOrder", &fit_order, "fitOrder/I");
  general_tree->Branch("stationNumber", &station_number, "stationNumber/I");
  general_tree->Branch("startTime", &start_time, "startTime/i");
  general_tree->Branch("endTime", &end_time, "endTime/i");
  general_tree->Branch("fit_isUsingResid", &fit_isUsingResid, "fit_isUsingResid/O");
  general_tree->SetDirectory(&f);
  general_tree->Fill();
  general_tree->Write();

  TTree fitCoeffs_tree("coeffs_tree", "coeffs_tree");
  std::vector<float> coeff(fit_order+1);
  std::vector<float> *p_coeff = &coeff;
  fitCoeffs_tree.Branch("coeff", "std::vector<float>", &p_coeff);
  fitCoeffs_tree.SetDirectory(&f);

  TTree aveResidGraph_tree("aveResidGraph_tree", "aveResidGraph_tree");
  TGraphErrors *p_aveResidGraph = nullptr;
  aveResidGraph_tree.Branch("aveResidGraph", "TGraphErrors", &p_aveResidGraph);
  aveResidGraph_tree.SetDirectory(&f);

  TTree chisqValidation_tree("chisqValidation_tree", "chisqValidation_tree");
  std::vector<bool> sampChisqPerDOF(mattak::k::num_lab4_samples);
  std::vector<bool> *p_sampChisqPerDOF = &sampChisqPerDOF;
  bool channelAveChisqPerDOF;
  chisqValidation_tree.Branch("sampChisqPerDOF", "std::vector<bool>", &p_sampChisqPerDOF);
  chisqValidation_tree.Branch("channelAveChisqPerDOF", &channelAveChisqPerDOF, "channelAveChisqPerDOF/O");
  chisqValidation_tree.SetDirectory(&f);

  const int nThresholds = 4; // 4 box frame thresholds for each channel histogram
  TTree residValidation_tree("residValidation_tree", "residValidation_tree");
  std::vector<bool> residOutOfBoxFrame(nThresholds);
  std::vector<bool> *p_residOutOfBoxFrame = &residOutOfBoxFrame;
  residValidation_tree.Branch("residOutOfBoxFrame", "std::vector<bool>", &p_residOutOfBoxFrame);
  residValidation_tree.SetDirectory(&f);

  for (int iChan = 0; iChan < mattak::k::num_radiant_channels; iChan++)
  {
    channelAveChisqPerDOF = isBad_channelAveChisqPerDOF[iChan];

    for (int iSamp = 0; iSamp < mattak::k::num_lab4_samples; iSamp++)
    {
      sampChisqPerDOF[iSamp] = isBad_sampChisqPerDOF[iChan][iSamp];

      for (int iOrder = 0; iOrder <= fit_order; iOrder++)
      {
        coeff[iOrder] = getFitCoeff(iChan, iSamp, iOrder);
      }
      fitCoeffs_tree.Fill();
    }

    p_aveResidGraph = graph_residAve[iChan];
    aveResidGraph_tree.Fill();

    for (int i = 0; i < nThresholds; i++)
    {
      residOutOfBoxFrame[i] = isResidOutOfBoxFrame[iChan][i];
    }
    residValidation_tree.Fill();

    chisqValidation_tree.Fill();
  }
  fitCoeffs_tree.Write();
  aveResidGraph_tree.Write();
  residValidation_tree.Write();
  chisqValidation_tree.Write();

  std::cout << "\nAll voltage calibration constants saved in file: " << outFileName << "\n\n" << std::endl;
  f.Close();
}

void mattak::VoltageCalibration::readFitCoeffsFromFile(const char * inFile, bool cache_tables)
{
  //
  // Get information about the bias scan and fit coefficients from the input root file
  //
  hasBiasScanData = false;

   TFile * f  = TFile::Open(inFile);
   if (!readFitCoeffsFromFile(f, cache_tables))
   {
     std::cerr << "Trouble reading from " << inFile;
   }
   delete f ;
}

bool mattak::VoltageCalibration::readFitCoeffsFromFile(TFile * inputFile, bool cache_tables)
{
  has_cache_tables_ = cache_tables;

  if (!inputFile->IsOpen()) return false;

  TTree *general_tree = (TTree*)inputFile->Get("general_tree");
  if (!general_tree) return false;


  if (general_tree->SetBranchAddress("fitOrder", &fit_order) < 0 ) return false;
  if (general_tree->SetBranchAddress("stationNumber", &station_number) < 0) return false ;
  if (general_tree->SetBranchAddress("startTime", &start_time) < 0 ) return false;
  if (general_tree->SetBranchAddress("endTime", &end_time) < 0 ) return false;
  if (general_tree->SetBranchAddress("fit_isUsingResid", &fit_isUsingResid) < 0 ) return false;
  general_tree->GetEntry(0);

  TTree *fitCoeffs_tree = (TTree*)inputFile->Get("coeffs_tree");
  if (!fitCoeffs_tree) return false;
  std::vector<float> coeff(fit_order + 1);
  std::vector<float> *p_coeff = &coeff;
  if (fitCoeffs_tree->SetBranchAddress("coeff", &p_coeff) < 0) return false;

  TTree *aveResidGraph_tree = (TTree*)inputFile->Get("aveResidGraph_tree");
  if (!aveResidGraph_tree)
  {
    printf("\nFILE READING ERROR: You are probably using an old calibration file, please use one with the newest version of voltage calibration instead!\n");
    return false;
  }
  TGraphErrors *p_aveResidGraph;
  if (aveResidGraph_tree->SetBranchAddress("aveResidGraph", &p_aveResidGraph) < 0) return false;
  aveResidGraph_tree->GetEntry(0);

  TTree *chisqValidation_tree = (TTree*)inputFile->Get("chisqValidation_tree");
  if (!chisqValidation_tree) return false;
  std::vector<bool> sampChisqPerDOF(mattak::k::num_lab4_samples);
  std::vector<bool> *p_sampChisqPerDOF = &sampChisqPerDOF;
  bool channelAveChisqPerDOF;
  if (chisqValidation_tree->SetBranchAddress("sampChisqPerDOF", &p_sampChisqPerDOF) < 0) return false;
  if (chisqValidation_tree->SetBranchAddress("channelAveChisqPerDOF", &channelAveChisqPerDOF) < 0) return false;

  const int nThresholds = 4; // 4 box frame thresholds for each channel histogram
  TTree *residValidation_tree = (TTree*)inputFile->Get("residValidation_tree");
  if (!residValidation_tree) return false;
  std::vector<bool> residOutOfBoxFrame(nThresholds);
  std::vector<bool> *p_residOutOfBoxFrame = &residOutOfBoxFrame;
  if (residValidation_tree->SetBranchAddress("residOutOfBoxFrame", &p_residOutOfBoxFrame) < 0) return false;

  if (fit_order < max_voltage_calibration_fit_order)
  {
    printf("\n%d-degree polynomials were used for the fitting...\n", fit_order);
    printf("SUGGESTION: Using order 9 (default) is highly suggested for getting better calibration results.\n");
  }

  if (!fit_isUsingResid) printf("\nNOTICE: 'fit_isUsingResid' is FALSE => The extra term residual function is not used!\n");

  bool isBadFit = false;
  int nBadChannels = 0;
  int nBadChannels_box = 0;
  int nChannels_badSamplesFound = 0;
  std::vector<bool> badSamplesFound(mattak::k::num_radiant_channels);

  for(int iChan = 0; iChan < mattak::k::num_radiant_channels; iChan++)
  {
    int nBadSamples = 0;
    badSamplesFound[iChan] = false;

    chisqValidation_tree->GetEntry(iChan);
    isBad_channelAveChisqPerDOF[iChan] = channelAveChisqPerDOF;

    if (isBad_channelAveChisqPerDOF[iChan])
    {
      printf("BAD FITTING WARNING: The average chi2/DOF over all samples of CH%d is greater than 6.0!!!\n", iChan);
      nBadChannels++;
    }

    if (nBadChannels > 2)
    {
      isBadFit = true;
    }

    fit_coeffs[iChan].clear();
    fit_coeffs[iChan].resize((fit_order + 1) * mattak::k::num_lab4_samples, 0);

    for(int iSamp = 0; iSamp < mattak::k::num_lab4_samples; iSamp++)
    {
      isBad_sampChisqPerDOF[iChan][iSamp] = sampChisqPerDOF[iSamp];
      if (!isBad_channelAveChisqPerDOF[iChan] && isBad_sampChisqPerDOF[iChan][iSamp])
      {
        printf("BAD FITTING WARNING: chi2/DOF of sample %d in CH%d is greater than 30.0!!!\n", iSamp, iChan);
        nBadSamples++;
      }

      fitCoeffs_tree->GetEntry(iChan*mattak::k::num_lab4_samples+iSamp);
      for(int iOrder = 0; iOrder < fit_order+1; iOrder++)
      {
        fit_coeffs[iChan][iSamp * (fit_order + 1) + iOrder] = coeff[iOrder];
      }
    }

    if (nBadSamples > 20)
    {
      badSamplesFound[iChan] = true;
    }

    nChannels_badSamplesFound += badSamplesFound[iChan];

    if (nChannels_badSamplesFound > 2)
    {
      isBadFit = true;
    }

    residValidation_tree->GetEntry(iChan);
    for (int i = 0; i < nThresholds; i++)
    {
      isResidOutOfBoxFrame[iChan][i] = residOutOfBoxFrame[i];
      if (isResidOutOfBoxFrame[iChan][i])
      {
        if (i == 0) printf("BAD FITTING WARNING: Some residuals in CH%d are beyond the SMALL BOX FRAME upper threshold (25 adu)!!!\n", iChan);
        else if (i == 1) printf("BAD FITTING WARNING: Some residuals in CH%d are below the SMALL BOX FRAME lower threshold (-25 adu)!!!\n", iChan);
        else if (i == 2) printf("BAD FITTING WARNING: Some residuals in CH%d are beyond the BIG BOX FRAME upper threshold (50 adu)!!!\n", iChan);
        else printf("BAD FITTING WARNING: Some residuals in CH%d are below the BIG BOX FRAME lower threshold (-50 adu)!!!\n", iChan);
        nBadChannels_box++;
      }
    }

    if (nBadChannels_box > 2)
    {
      isBadFit = true;
    }

    // Average residuals for each channel
    graph_residAve[iChan] = p_aveResidGraph;

    if (fit_isUsingResid)
    {
      // Interpolating the average residuals
      int npoints_residGraph = graph_residAve[iChan]->GetN();
      graph_residAve[iChan]->SetBit(TGraph::kIsSortedX);  // We can do that because our data are sorted. Makes later Eval calls faster
      const double dV = graph_residAve[iChan]->GetPointX(1) - graph_residAve[iChan]->GetPointX(0);
      resid_volt[iChan].resize(npoints_residGraph * 2 - 1);
      resid_adc[iChan].resize(npoints_residGraph * 2 - 1);
      nResidPoints[iChan] = resid_volt[iChan].size();

      for (int i = 0; i < npoints_residGraph; i++)
      {
        resid_volt[iChan][i*2] = graph_residAve[iChan]->GetPointX(i);
        resid_adc[iChan][i*2] = graph_residAve[iChan]->GetPointY(i);
      }
      for (int i = 0; i < npoints_residGraph-1; i++)
      {
        //usampling by a factor of 2
        resid_volt[iChan][i*2+1] = resid_volt[iChan][i*2] + dV/2;
        resid_adc[iChan][i*2+1] = graph_residAve[iChan]->Eval(resid_volt[iChan][i*2+1]);
      }
    }

  }

  hasBiasScanData = false;
  inputFile->Close();

  if (isBadFit)
  {
    printf("BAD CALIBRATION FIT detected, refuse to use this file... ABORT...");
    return false;
  }

  if (has_cache_tables_)
  {
    for(int channel = 0; channel < mattak::k::num_radiant_channels; channel++)
    {
      for(int sample = 0; sample < mattak::k::num_lab4_samples; sample++)
      {
        cached_adc_tables_[channel][sample].resize(nResidPoints[channel]);
        const double *params = &fit_coeffs[channel][sample * (fit_order + 1)];
        const double* adcTable = adcTablePerSample(fit_order, nResidPoints[channel], params, &resid_volt[channel][0], &resid_adc[channel][0]);
        for (int idx = 0; idx < nResidPoints[channel]; idx++)
        {
          cached_adc_tables_[channel][sample][idx] = adcTable[idx];
        }
        delete [] adcTable;
      }
    }
  }
  return true;
}


#else
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

double * mattak::applyVoltageCalibration (int nSamples_wf, const int16_t * in, double * out, int start_window, bool isOldFirmware, int fit_order,
                            int nResidPoints, const double * packed_fit_params, bool isUsingResid, const double * packed_aveResid_volt, const double * packed_aveResid_adc)

{
  if (!out) out = new double[nSamples_wf];

  if (nSamples_wf % mattak::k::radiant_window_size)
  {
    std::cerr << "Not multiple of window size!" << std::endl;
    return 0;
  }

  if (fit_order > mattak::max_voltage_calibration_fit_order)
  {
    std::cerr << "Fit Order Too Large!" << std::endl;
    return 0;
  }

  int nSamplesPerGroup = mattak::k::num_radiant_samples;
  int nWindowsPerGroup = mattak::k::radiant_windows_per_buffer;
  int isamp_lab4, isamp_A, isamp_B;

  isamp_A = (start_window >= nWindowsPerGroup) * nSamplesPerGroup;

  if (isOldFirmware)
  {
    nSamplesPerGroup /= 2;
    nWindowsPerGroup /= 2;
  }

  for (int i = 0; i < nSamples_wf; i++)
  {
    isamp_B = (i + start_window * mattak::k::radiant_window_size) % nSamplesPerGroup;

    if (isOldFirmware)
    {
      isamp_B += (i >= nSamplesPerGroup) * nSamplesPerGroup;
    }

    isamp_lab4 = isamp_A + isamp_B;

    // creating a pointer which points to the first parameter of this particular sample
    const double *params = packed_fit_params + isamp_lab4 * (fit_order + 1);
    // const double *params = &packed_fit_params[isamp_lab4 * (fit_order + 1)][0];

    double adc = in[i];
    double *adcTable = nullptr;
    if (isUsingResid)
    {
      adcTable = adcTablePerSample(fit_order, nResidPoints, params, packed_aveResid_volt, packed_aveResid_adc);
      out[i] = adcToVolt(adc, nResidPoints, packed_aveResid_volt, adcTable);
    }
    else
    {
      // When we perform a calibration without residuals, we directly fit f(ADC) -> V
      out[i] = evalPars(adc, fit_order, params);
    }

    if (adcTable) delete [] adcTable;
  }

  return out;
}

double * mattak::applyVoltageCalibration(
  int nSamples_wf, const int16_t * in, double * out, int start_window, bool isOldFirmware, int nResidPoints, const double * voltage_table, const std::array<std::vector<double>, 4096>* adc_table)

{
  if (!out) out = new double[nSamples_wf];

  if (nSamples_wf % mattak::k::radiant_window_size)
  {
    std::cerr << "Not multiple of window size!" << std::endl;
    return 0;
  }

  int nSamplesPerGroup = mattak::k::num_radiant_samples;
  int nWindowsPerGroup = mattak::k::radiant_windows_per_buffer;
  int sample_lab4, sample_A, sample_B;

  sample_A = (start_window >= nWindowsPerGroup) * nSamplesPerGroup;

  if (isOldFirmware)
  {
    nSamplesPerGroup /= 2;
    nWindowsPerGroup /= 2;
  }


  for (int sample_wf = 0; sample_wf < nSamples_wf; sample_wf++)
  {
    sample_B = (sample_wf + start_window * mattak::k::radiant_window_size) % nSamplesPerGroup;

    if (isOldFirmware)
    {
      sample_B += (sample_wf >= nSamplesPerGroup) * nSamplesPerGroup;
    }

    sample_lab4 = sample_A + sample_B;

    // Get table for correct sample from pointer. Get pointer to point to address of first ele of table
    const double * adc_table_sample = &adc_table->at(sample_lab4)[0];
    out[sample_wf] = adcToVolt(in[sample_wf], nResidPoints, voltage_table, adc_table_sample);
  }

  return out;
}


// PYTHON BINDING

#ifdef MATTAK_NOROOT

namespace py = pybind11;

static py::array_t<double> apply_voltage_calibration(py::buffer in, int start_window, bool isOldFirmware, int order, int nResidPoints, py::buffer packed_coeffs, bool isUsingResid, py::buffer packed_aveResid_volt, py::buffer packed_aveResid_adc)
{
  //make sure in is int16_t, get n

  auto in_info = in.request();
  if ( (in_info.format !=  py::format_descriptor<int16_t>::format()) || (in_info.ndim != 1) || (in_info.strides[0]!=0) )
  {
    std::cerr << "in must be of type int16_t, 1 dim, no stride" << std::endl;
    return py::none();
  }

  int N = in_info.size;

  auto packed_info = packed_coeffs.request();

  if (packed_info.format != py::format_descriptor<double>::format() || packed_info.ndim != 1 || packed_info.strides[0] != 0)
  {
    std::cerr << "packed_coeffs must be of type double, 1 dim, no stride" << std::endl;
    return py::none();
  }

  if (packed_info.size != (order+1) * mattak::k::num_lab4_samples)
  {
    std::cerr << "Packed coeffs not right size" << std::endl;
    return py::none();
  }

  auto ret = py::array_t<double>(N);

  if (!mattak::applyVoltageCalibration(N, (int16_t*) in.ptr(), (double*) ret.ptr(), start_window, isOldFirmware, order, nResidPoints, (double*) packed_coeffs.ptr(), isUsingResid, (double*) packed_aveResid_volt.ptr(), (double*) packed_aveResid_adc.ptr()))
  {
    return py::none();
  }

  return ret;
}

PYBIND11_MODULE(mattak, m)
{
  m.def("apply_voltage_calibration", &apply_voltage_calibration, "Apply voltage Calibration");
}

#endif
