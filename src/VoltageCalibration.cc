#include "mattak/VoltageCalibration.h"
#include <iostream>
#include <stdio.h>


#ifdef MATTAK_VECTORIZE
#include "vectorclass/vectorclass.h"
#endif

#ifdef LIBRNO_G_SUPPORT
#include "rno-g.h"
#endif


static double evalPars(double x, int order, const double * p)
{
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
  const double *voltTable = resid_volt;
  double *adcTable = new double[npoints];

  for (int i = 0; i < npoints; i++)
  {
    adcTable[i] = evalPars(voltTable[i], order, par) + resid_adc[i];
  }

  return adcTable;
}

static double adcToVolt(double in_adc, int npoints, const double * resid_volt, const double * resid_adc)
{
  const double *volt_array = resid_volt;
  const double *adc_array = resid_adc;
  double m;
  double out_volt = 0;

  // If in_adc is zero, out_volt is zero
  if (in_adc == 0) return out_volt;

  // If in_adc is out of range...
  if (in_adc < adc_array[0])
  {
    m = (volt_array[1] - volt_array[0])/(adc_array[1] - adc_array[0]);
    out_volt = volt_array[0] + (in_adc - adc_array[0]) * m;
    return out_volt;
  }
  if (in_adc > adc_array[npoints-1])
  {
    m = (volt_array[npoints-1] - volt_array[npoints-2])/(adc_array[npoints-1] - adc_array[npoints-2]);
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
        m = (volt_array[i+1] - volt_array[i])/(adc_array[i+1] - adc_array[i]);
        out_volt = volt_array[i] + (in_adc - adc_array[i]) * m;
        return out_volt;
      }
    }
  }

  return out_volt;
}

static bool inRange(unsigned low, unsigned high, unsigned x)
{
  return (low <= x && x <= high);
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
{
  setupFromTree(tree,branch_name,vref,fit_order,min,max,isUsingResid);
}

//almost copy-pasted from raw file reader...
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
  recalculateFits(fit_order, min, max, vref, isUsingResid);
}

mattak::VoltageCalibration::VoltageCalibration(const char * bias_scan_file, double vref, int fit_order, double min, double max, bool isUsingResid)
{

  /*
     Work out the extension of the file bias_scan_file.
     If the extension after the last period ('.') is '.root',
     assume the bias scan file has been loaded in ROOT format.
     Otherwise they must be in a file format that required librno-g to read.
  */

  const char * suffix = strrchr(bias_scan_file,'.');

  if (!strcmp(suffix,".root") && !(*(suffix + sizeof(".root")-1)))
  {
    TFile f(bias_scan_file);
    TTree * t = (TTree*) f.Get("pedestals");
    if (!t)
    {
      std::cerr << "Could not open tree pedestals in " << bias_scan_file << std::endl;
      return;
    }
    setupFromTree(t,"pedestals", vref, fit_order, min, max, isUsingResid);
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
  recalculateFits(fit_order, min, max, vref, isUsingResid);
  rno_g_close_handle(&h);
#endif

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
  gErrorIgnoreLevel = kFatal;

  hasBiasScanData = true;

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

  int dacType;
  int nResidSets[2]; // 2 DAC
  std::array<std::vector<double>, 2> residAve_volt;
  std::array<std::vector<double>, 2> residAve_adc;
  for (int j = 0; j < 2; j++)
  {
    nResidSets[j] = 0;
    for (int i = 0; i < scanSize(); i++)
    {
      residAve_volt[j].push_back(0);
      residAve_adc[j].push_back(0);
    }
  }


  TLinearFitter fit(1, formula[order]);
  for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++)
  {
    fit_coeffs[ichan].clear();
    if ( (mask & (1 << ichan))  == 0) continue;

    //check if channel seems broken by looking for zeroes

    fit_coeffs[ichan].resize((order+1) * mattak::k::num_lab4_samples, 0);

    int nbroken = 0;

    for (int i = 0; i < mattak::k::num_lab4_samples; i++)
    {
      graph[ichan][i] = new TGraph();
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
            if (fit_isUsingResid) graph[ichan][i]->SetPoint(graph[ichan][i]->GetN(),v,adc);
            else graph[ichan][i]->SetPoint(graph[ichan][i]->GetN(),adc,v);
            jmax = j;
          }
          last_adc = adc;

          if (!jmin) jmin = j;
        }
      }

      int npoints = graph[ichan][i]->GetN();
      double *data_adc;
      double *data_v;

      if (fit_isUsingResid)
      {
        data_adc = graph[ichan][i]->GetY();
        data_v = graph[ichan][i]->GetX();
      }
      else
      {
        data_adc = graph[ichan][i]->GetX();
        data_v = graph[ichan][i]->GetY();
      }

      for (int i = 0; i < npoints; i++)
      {
        if (fit_isUsingResid) fit.AddPoint(&data_v[i], data_adc[i]);
        else fit.AddPoint(&data_adc[i], data_v[i]);
      }

      if (vref) fit.FixParameter(0, 0);

      if (nzero > 1) nbroken++;
      else fit.Eval();

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

      if (fit_isUsingResid)
      {
        // Sum up all the ADC(V) residuals
        dacType = ichan >= mattak::k::num_radiant_channels / 2;
        nResidSets[dacType] ++;

        if (residAve_volt[dacType].size() == scanSize() || residAve_adc[dacType].size() == scanSize())
        {
          residAve_volt[dacType].resize(npoints);
          residAve_adc[dacType].resize(npoints);
        }

        for (int j = 0 ; j < npoints; j++)
        {
          double adc = data_adc[j];
          double v = data_v[j];
          double adcResid = adc - evalPars(v, fit_order, &fit_coeffs[ichan][i * (order+1)]);
          residAve_volt[dacType][j] += v;
          residAve_adc[dacType][j] += adcResid;
        }
      }

    }
    if (nbroken) printf("WARNING: Channel %d seems to have %d broken samples?\n", ichan, nbroken);
  }

  // Constants for the residual histograms
  const int nBinsY = 22;
  const int histLowY = -55;
  const int histHighY = 55;
  int nBinsX;
  double histLowX;
  double histHighX;

  // Calculate the average ADC(V) residuals and make a TGraph for each DAC
  for (int j = 0; j < 2; j++)
  {
    TString graphNameTitle = TString::Format("aveResid_dac%d", j+1);
    graph_residAve[j] = new TGraph();
    graph_residAve[j]->SetNameTitle(graphNameTitle);
    graph_residAve[j]->GetXaxis()->SetTitle("VBias [Volt]");
    graph_residAve[j]->GetYaxis()->SetTitle("ADC Residual");

    if (fit_isUsingResid)
    {
      int npoints_residGraph = residAve_volt[j].size();
      for (int ipoint = 0; ipoint < npoints_residGraph; ipoint++)
      {
        residAve_volt[j][ipoint] /= nResidSets[j];
        residAve_adc[j][ipoint] /= nResidSets[j];
        graph_residAve[j]->SetPoint(graph_residAve[j]->GetN(), residAve_volt[j][ipoint], residAve_adc[j][ipoint]);
      }

      resid_volt[j].resize(npoints_residGraph*2-1);
      resid_adc[j].resize(npoints_residGraph*2-1);
      nResidPoints[j] = resid_volt[j].size();

      //
      // Interpolating the average residuals
      //
      for (int i = 0; i < npoints_residGraph; i++)
      {
        resid_volt[j][i*2] = graph_residAve[j]->GetPointX(i);
        resid_adc[j][i*2] = graph_residAve[j]->GetPointY(i);
      }
      for (int i = 0; i < npoints_residGraph-1; i++)
      {
        resid_volt[j][i*2+1] = (resid_volt[j][i*2] + resid_volt[j][i*2+2])/2;
        resid_adc[j][i*2+1] = resid_adc[j][i*2] + (resid_volt[j][i*2+1] - resid_volt[j][i*2])*(resid_adc[j][i*2+2] - resid_adc[j][i*2])/(resid_volt[j][i*2+2] - resid_volt[j][i*2]);
      }

      // Residual histograms for the 2 DACs
      TString histNameTitle = TString::Format("residHist_dac%d", j+1);
      nBinsX = npoints_residGraph;
      histLowX = graph_residAve[j]->GetPointX(0);
      histHighX = graph_residAve[j]->GetPointX(npoints_residGraph-1);

      hist_resid[j] = new TH2S(histNameTitle, histNameTitle, nBinsX, histLowX, histHighX, nBinsY, histLowY, histHighY);
      hist_resid[j]->GetXaxis()->SetTitle("VBias [Volt]");
      hist_resid[j]->GetYaxis()->SetTitle("ADC Residual");
    }

  }


  //
  // Calculate max deviation and chi squared
  //
  std::cout << "\nCalculating max deviation and chi squared for fit quality validation..." << std::endl;

  std::vector<double> aveChisq;
  aveChisq.resize(mattak::k::num_radiant_channels);

  for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++)
  {
    std::cout << "Channel " << ichan;
    if ( (mask & (1 << ichan))  == 0) continue;

    aveChisq[ichan] = 0;

    int dacType = ichan >= mattak::k::num_radiant_channels / 2;

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

      int npoints = graph[ichan][i]->GetN();
      double *data_adc;
      double *data_v;
      double *residTablePerSamp_adc;

      if (fit_isUsingResid)
      {
        data_adc = graph[ichan][i]->GetY();
        data_v = graph[ichan][i]->GetX();
        residTablePerSamp_adc = adcTablePerSample(fit_order, nResidPoints[dacType], &fit_coeffs[ichan][i * (order+1)], &resid_volt[dacType][0], &resid_adc[dacType][0]);
      }
      else
      {
        data_adc = graph[ichan][i]->GetX();
        data_v = graph[ichan][i]->GetY();
      }

      // Calculate max deviation and chi squared
      for (int j = 0 ; j < npoints; j++)
      {
        double adc = data_adc[j];
        double v_meas = data_v[j];
        double v_pred;

        if (fit_isUsingResid) v_pred = adcToVolt(adc, nResidPoints[dacType], &resid_volt[dacType][0], &residTablePerSamp_adc[0]);
        else v_pred = evalPars(adc, fit_order, &fit_coeffs[ichan][i * (order+1)]);

        double delta = fabs(v_meas-v_pred);
        if (delta > fit_maxerr[ichan][i]) fit_maxerr[ichan][i] = delta;
        fit_chisq[ichan][i] += ( delta*delta/(0.002*0.002) );

        if (fit_isUsingResid)
        {
          double histX = v_meas;
          double histY = adc - (evalPars(v_meas, fit_order, &fit_coeffs[ichan][i * (order+1)]) + graph_residAve[dacType]->GetPointY(j));
          hist_resid[dacType]->Fill(histX, histY);
        }
      }

      if (fit_chisq[ichan][i]/fit_ndof[ichan][i] > 30.0)
      {
        channelHasBadFit = true;
        isBad_sampChisqPerDOF[ichan][i] = true;
        badFit.push_back(i);
      }

      aveChisq[ichan] = aveChisq[ichan] + (fit_chisq[ichan][i]/fit_ndof[ichan][i]);

      delete residTablePerSamp_adc;
    }

    // chi2 check for fit quality validation
    aveChisq[ichan] /= mattak::k::num_lab4_samples;
    if (aveChisq[ichan] > 6.0)
    {
      isBad_channelAveChisqPerDOF[ichan] = true;
      printf("\nBAD FITTING WARNING: The average chi2/DOF over all samples of CH%d is %f (> 6.0)!!!\n", ichan, aveChisq[ichan]);
    }

    if (aveChisq[ichan] <= 6.0 && channelHasBadFit)
    {
      for (int samp = 0; samp < badFit.size(); samp++)
      {
        int bad = badFit[samp];
        printf("\nBAD FITTING WARNING: chi2/DOF of sample %d in CH%d is %f (> 30.0)!!!\n", bad, ichan, fit_chisq[ichan][bad]/fit_ndof[ichan][bad]);
      }
    }

  }


  // Residual histograms and the box frame check for fit quality validation
  if (fit_isUsingResid)
  {
    TBox *smallBox[2];
    TBox *bigBox[2];

    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        isResidOutOfBoxFrame[i][j] = false;
      }

      nBinsX = hist_resid[i]->GetNbinsX();
      histLowX = hist_resid[i]->GetXaxis()->GetXmin();
      histHighX = hist_resid[i]->GetXaxis()->GetXmax();

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

      smallBox[i] = new TBox(smallBoxX1, smallBoxY1, smallBoxX2, smallBoxY2);
      smallBox[i]->SetFillStyle(0);
      smallBox[i]->SetLineStyle(9);
      smallBox[i]->SetLineWidth(3);
      smallBox[i]->SetLineColor(2);
      bigBox[i] = new TBox(bigBoxX1, bigBoxY1, bigBoxX2, bigBoxY2);
      bigBox[i]->SetFillStyle(0);
      bigBox[i]->SetLineStyle(9);
      bigBox[i]->SetLineWidth(3);
      bigBox[i]->SetLineColor(2);

      bool aboveSmallBoxY2 = false;
      bool belowSmallBoxY1 = false;
      bool aboveBigBoxY2 = false;
      bool belowBigBoxY1 = false;

      for (int binNumberX = 1; binNumberX <= graph_residAve[i]->GetN(); binNumberX++)
      {
        // Small Box Check
        if (binNumberX > smallBoxBinX1 && binNumberX <= smallBoxBinX2)
        {
          if (hist_resid[i]->GetBinContent(binNumberX, smallBoxBinY2) > 1 && !aboveSmallBoxY2)
          {
            aboveSmallBoxY2 = true;
            isResidOutOfBoxFrame[i][0] = true;
            printf("\nBAD FITTING WARNING: Some residuals in DAC-%d go beyond the first upper threshold (> 25 adu)!!!\n", i+1);
          }
          if (hist_resid[i]->GetBinContent(binNumberX, smallBoxBinY1) > 1 && !belowSmallBoxY1)
          {
            belowSmallBoxY1 = true;
            isResidOutOfBoxFrame[i][1] = true;
            printf("\nBAD FITTING WARNING: Some residuals in DAC-%d go below the first lower threshold (< -25 adu)!!!\n", i+1);
          }
        }

        // Big Box Check
        if (hist_resid[i]->GetBinContent(binNumberX, nBinsY) > 1 && !aboveBigBoxY2)
        {
          aboveBigBoxY2 = true;
          isResidOutOfBoxFrame[i][2] = true;
          printf("\nBAD FITTING WARNING: Some residuals in DAC-%d go beyond the second upper threshold (> 50 adu)!!!\n", i+1);
        }
        if (hist_resid[i]->GetBinContent(binNumberX, 1) > 1 && !belowBigBoxY1)
        {
          belowBigBoxY1 = true;
          isResidOutOfBoxFrame[i][3] = true;
          printf("\nBAD FITTING WARNING: Some residuals in DAC-%d go below the second lower threshold (< -50 adu)!!!\n", i+1);
        }

        if (aboveBigBoxY2 && belowBigBoxY1 && aboveSmallBoxY2 && belowSmallBoxY1) break;
      }

      hist_resid[i]->GetListOfFunctions()->Add(bigBox[i]);
      hist_resid[i]->GetListOfFunctions()->Add(smallBox[i]);
    }

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
  TGraph *g = new TGraph();
  g->SetName(Form("gsample_inverse_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time));
  g->SetTitle(Form("Station %d Ch %d sample %d [%d-%d]   %s", station_number, chan, samp, start_time, end_time, resid ? "(residuals)" : ""));
  g->GetXaxis()->SetTitle("VBias [Volt]");
  g->GetYaxis()->SetTitle(resid ? "ADC Residual" : "ADC");

  if (!hasBiasScanData) { printf("\nWARNING: Need to get data from a bias scan file in order to make graphs!\n"); return g; }

  if (!fit_isUsingResid) { printf("\nWARNING: Plots can only be made with function 'makeAdjustedInverseGraph()' when 'fit_isUsingResid' is TRUE!\n"); return g; }

  int dacType = chan >= mattak::k::num_radiant_channels / 2;

  int npoints = graph[chan][samp]->GetN();
  double *data_adc = graph[chan][samp]->GetY();
  double *data_v = graph[chan][samp]->GetX();

  TF1 * fn = new TF1(Form("fsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time), formula[fit_order], fit_min, fit_max, TF1::EAddToList::kNo);
  fn->SetParameters(getFitCoeffs(chan,samp));
  fn->SetRange(data_v[0],data_v[npoints-1]);
  fn->SetLineColor(2);

  for (int j = 0; j < npoints; j++)
  {
    double v = data_v[j];
    double adc = data_adc[j] - graph_residAve[dacType]->GetPointY(j);

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
  TGraph *g = new TGraph();
  g->SetName(Form("gsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time));
  g->SetTitle(Form("Station %d Ch %d sample %d [%d-%d], #chi^{2}= %g   %s", station_number, chan, samp, start_time, end_time, fit_chisq[chan][samp], resid ? "(residuals)" : ""));
  g->GetXaxis()->SetTitle("ADC");
  g->GetYaxis()->SetTitle(resid ? "(VBias - Predicted VBias) [Volt]" : "VBias [Volt]");

  if (!hasBiasScanData) { printf("\nWARNING: Need to get data from a bias scan file in order to make graphs!\n"); return g; }

  int dacType = chan >= mattak::k::num_radiant_channels / 2;

  int npoints = graph[chan][samp]->GetN();
  double *data_adc;
  double *data_v;
  double *residTablePerSamp_adc;

  TF1 *fn;

  if (fit_isUsingResid)
  {
    data_adc = graph[chan][samp]->GetY();
    data_v = graph[chan][samp]->GetX();
    residTablePerSamp_adc = adcTablePerSample(fit_order, nResidPoints[dacType], getFitCoeffs(chan,samp), getPackedAveResid_volt(chan), getPackedAveResid_adc(chan));
  }
  else
  {
    data_adc = graph[chan][samp]->GetX();
    data_v = graph[chan][samp]->GetY();

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
      if (fit_isUsingResid) v -= adcToVolt(adc, nResidPoints[dacType], getPackedAveResid_volt(chan), residTablePerSamp_adc);
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

  TTree chisqValidation_tree("chisqValidation_tree", "chisqValidation_tree");
  std::vector<bool> sampChisqPerDOF(mattak::k::num_lab4_samples);
  std::vector<bool> *p_sampChisqPerDOF = &sampChisqPerDOF;
  bool channelAveChisqPerDOF;
  chisqValidation_tree.Branch("sampChisqPerDOF", "std::vector<bool>", &p_sampChisqPerDOF);
  chisqValidation_tree.Branch("channelAveChisqPerDOF", &channelAveChisqPerDOF, "channelAveChisqPerDOF/O");
  chisqValidation_tree.SetDirectory(&f);

  const int nThresholdsPerDAC = 4; // 4 box frame thresholds for each DAC
  TTree residValidation_tree("residValidation_tree", "residValidation_tree");
  std::vector<bool> residOutOfBoxFrame(nThresholdsPerDAC);
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
    chisqValidation_tree.Fill();
  }
  fitCoeffs_tree.Write();
  chisqValidation_tree.Write();

  for (int i_dac = 0; i_dac < 2; i_dac++)
  {
    for (int i = 0; i < nThresholdsPerDAC; i++)
    {
      residOutOfBoxFrame[i] = isResidOutOfBoxFrame[i_dac][i];
    }
    residValidation_tree.Fill();
  }
  residValidation_tree.Write();

  getAveResidGraph_dac1()->Write();
  getAveResidGraph_dac2()->Write();

  std::cout << "\n\nAll voltage calibration constants saved in file: " << outFileName << std::endl;
  f.Close();
}

void mattak::VoltageCalibration::readFitCoeffsFromFile(const char * inFile)
{
  //
  // Get information about the bias scan and fit coefficients from the input root file
  //
  hasBiasScanData = false;

  TFile *inputFile = new TFile(inFile);

  TTree *general_tree = (TTree*)inputFile->Get("general_tree");
  general_tree->SetBranchAddress("fitOrder", &fit_order);
  general_tree->SetBranchAddress("stationNumber", &station_number);
  general_tree->SetBranchAddress("startTime", &start_time);
  general_tree->SetBranchAddress("endTime", &end_time);
  general_tree->SetBranchAddress("fit_isUsingResid", &fit_isUsingResid);
  general_tree->GetEntry(0);

  TTree *fitCoeffs_tree = (TTree*)inputFile->Get("coeffs_tree");
  std::vector<float> coeff(fit_order+1);
  std::vector<float> *p_coeff = &coeff;
  fitCoeffs_tree->SetBranchAddress("coeff", &p_coeff);

  TTree *chisqValidation_tree = (TTree*)inputFile->Get("chisqValidation_tree");
  std::vector<bool> sampChisqPerDOF(mattak::k::num_lab4_samples);
  std::vector<bool> *p_sampChisqPerDOF = &sampChisqPerDOF;
  bool channelAveChisqPerDOF;
  chisqValidation_tree->SetBranchAddress("sampChisqPerDOF", &p_sampChisqPerDOF);
  chisqValidation_tree->SetBranchAddress("channelAveChisqPerDOF", &channelAveChisqPerDOF);

  const int nThresholdsPerDAC = 4; // 4 box frame thresholds for each DAC
  TTree *residValidation_tree = (TTree*)inputFile->Get("residValidation_tree");
  std::vector<bool> residOutOfBoxFrame(nThresholdsPerDAC);
  std::vector<bool> *p_residOutOfBoxFrame = &residOutOfBoxFrame;
  residValidation_tree->SetBranchAddress("residOutOfBoxFrame", &p_residOutOfBoxFrame);

  if (fit_order < max_voltage_calibration_fit_order)
  {
    printf("\n%d-degree polynomials were used for the fitting...\n", fit_order);
    printf("SUGGESTION: Using order 9 (default) is highly suggested for getting better calibration results.\n");
  }

  if (!fit_isUsingResid) printf("\nNOTICE: 'fit_isUsingResid' is FALSE => The extra term residual function is not used!\n");

  for(int iChan = 0; iChan < mattak::k::num_radiant_channels; iChan++)
  {
    chisqValidation_tree->GetEntry(iChan);
    isBad_channelAveChisqPerDOF[iChan] = channelAveChisqPerDOF;
    if (isBad_channelAveChisqPerDOF[iChan]) printf("\nBAD FITTING WARNING: The average chi2/DOF over all samples of CH%d is greater than 6.0!!!\n", iChan);

    fit_coeffs[iChan].clear();
    fit_coeffs[iChan].resize((fit_order+1)*mattak::k::num_lab4_samples, 0);
    for(int iSamp = 0; iSamp < mattak::k::num_lab4_samples; iSamp++)
    {
      isBad_sampChisqPerDOF[iChan][iSamp] = sampChisqPerDOF[iSamp];
      if (!isBad_channelAveChisqPerDOF[iChan] && isBad_sampChisqPerDOF[iChan][iSamp])
      printf("\nBAD FITTING WARNING: chi2/DOF of sample %d in CH%d is greater than 30.0!!!\n", iSamp, iChan);

      fitCoeffs_tree->GetEntry(iChan*mattak::k::num_lab4_samples+iSamp);
      for(int iOrder = 0; iOrder < fit_order+1; iOrder++)
      {
        fit_coeffs[iChan][iSamp*(fit_order+1)+iOrder] = coeff[iOrder];
      }
    }
  }

  for (int j = 0; j < 2; j++)
  {
    residValidation_tree->GetEntry(j);
    for (int i = 0; i < nThresholdsPerDAC; i++)
    {
      isResidOutOfBoxFrame[j][i] = residOutOfBoxFrame[i];
      if (isResidOutOfBoxFrame[j][i])
      {
        if (i == 0) printf("\nBAD FITTING WARNING: Some residuals in DAC-%d are beyond the SMALL BOX FRAME upper threshold (25 adu)!!!\n", j+1);
        else if (i == 1) printf("\nBAD FITTING WARNING: Some residuals in DAC-%d are below the SMALL BOX FRAME lower threshold (-25 adu)!!!\n", j+1);
        else if (i == 2) printf("\nBAD FITTING WARNING: Some residuals in DAC-%d are beyond the BIG BOX FRAME upper threshold (50 adu)!!!\n", j+1);
        else printf("\nBAD FITTING WARNING: Some residuals in DAC-%d are below the BIG BOX FRAME lower threshold (-50 adu)!!!\n", j+1);
      }
    }

    // Average residuals for both types of DACs
    TString graphNameTitle = TString::Format("aveResid_dac%d", j+1);
    graph_residAve[j] = (TGraph*)inputFile->Get(graphNameTitle);

    if (fit_isUsingResid)
    {
      // Interpolating the average residuals
      int npoints_residGraph = graph_residAve[j]->GetN();
      resid_volt[j].resize(npoints_residGraph*2-1);
      resid_adc[j].resize(npoints_residGraph*2-1);
      nResidPoints[j] = resid_volt[j].size();

      for (int i = 0; i < npoints_residGraph; i++)
      {
        resid_volt[j][i*2] = graph_residAve[j]->GetPointX(i);
        resid_adc[j][i*2] = graph_residAve[j]->GetPointY(i);
      }
      for (int i = 0; i < npoints_residGraph-1; i++)
      {
        resid_volt[j][i*2+1] = (resid_volt[j][i*2] + resid_volt[j][i*2+2])/2;
        resid_adc[j][i*2+1] = resid_adc[j][i*2] + (resid_volt[j][i*2+1] - resid_volt[j][i*2])*(resid_adc[j][i*2+2] - resid_adc[j][i*2])/(resid_volt[j][i*2+2] - resid_volt[j][i*2]);
      }
    }
  }

  inputFile->Close();
}


#else
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

double * mattak::applyVoltageCalibration (int N, const int16_t * in, double * out, int start_window, bool isOldFirmware, int fit_order,
                            int nResidPoints, const double * packed_fit_params, bool isUsingResid, const double * packed_aveResid_volt, const double * packed_aveResid_adc)

{
  if (!out) out = new double[N];

  if (N % mattak::k::radiant_window_size)
  {
    std::cerr << "Not multiple of window size!" << std::endl;
    return 0;
  }

  if (fit_order > mattak::max_voltage_calibration_fit_order)
  {
    std::cerr << "Fit Order Too Large!" << std::endl;
    return 0;
  }

  bool is2ndHalfWindows;
  int nwindows = N / mattak::k::radiant_window_size;
  int isamp;
  int i = 0;
  int j;

  if (start_window > 15)
  {
    is2ndHalfWindows = true;
    j = start_window - mattak::k::radiant_windows_per_buffer;
  }
  else
  {
    is2ndHalfWindows = false;
    j = start_window;
  }


  for (int iwindow = 0; iwindow < nwindows; iwindow++)
  {
    if (isOldFirmware)
    {
      // Old Firmware
      // Wrap around if j % 16 = 0
      if (j % mattak::k::radiant_windows_per_buffer == 0) j = 0;
      if (inRange(0,4,j) || inRange(8,12,j)) isamp = (j+3) * mattak::k::radiant_window_size;
      if (inRange(5,7,j) || inRange(13,15,j)) isamp = (j-5) * mattak::k::radiant_window_size;
    }
    else
    {
      // New Firmware
      // Wrap around if j % 16 = 0
      if (j % mattak::k::radiant_windows_per_buffer == 0) j = 0;
      if (inRange(0,12,j)) isamp = (j+3) * mattak::k::radiant_window_size;
      if (inRange(13,15,j)) isamp = (j-13) * mattak::k::radiant_window_size;
    }

    if (is2ndHalfWindows) isamp += mattak::k::num_radiant_samples;

    j++;

#ifndef MATTAK_VECTORIZE
    for (int k = 0; k < mattak::k::radiant_window_size; k++)
    {
      const double *params = packed_fit_params + isamp * (fit_order+1);

      double *residTablePerSamp_adc;
      if (isUsingResid) residTablePerSamp_adc = adcTablePerSample(fit_order, nResidPoints, params, packed_aveResid_volt, packed_aveResid_adc);

      double adc = in[i];
      if (isUsingResid) out[i] = adcToVolt(adc, nResidPoints, packed_aveResid_volt, residTablePerSamp_adc);
      else out[i] = evalPars(adc, fit_order, params);

      isamp++;
      i++;

      delete residTablePerSamp_adc;
    }
#else

// Vectorized version, on intel x86x64 anyway. Optimized for AVX2.
#define VEC_SIZE 4
#define VECD vec4d
#define VECI vec4q
#define VEC_INCR VECI(0,1,2,3)
#define VEC_N mattak::k::radiant_window_size / (VEC_SIZE * VEC_UNROLL)

    VECD v[VEC_UNROLL];
    VECI vin[VEC_UNROLL];
    VECD x[VEC_UNROLL];
    VECI idx[VEC_UNROLL];


    for (int k = 0; k < vec_N; k++)
    {
      int iout = i; // for simplicity

      for (int u = 0; u < VEC_UNROLL; u++)
      {
        vin[u].load(in + i);
        x[u] = vin[u];
        idx[u]= VEC_INCR + isamp;
        v[u] = lookup < mattak::max_voltage_calibration_fit_order * mattak::k::num_lab4_samples > (idx, packed_fit_params);
        idx[u] *= (fit_order+1);
        isamp += VEC_SIZE;
        i += VEC_SIZE;
      }

      for (int j = fit_order-1; j >= 0; j--)
      {
        for (int u = 0; u < VEC_UNROLL; u++)
        {
          idx[u]++;
          v[u] = v[u] * x[u] + lookup < mattak::max_voltage_calibration_fit_order * mattak::k::num_lab4_samples > (idx, packed_fit_params);
        }
      }

      for (int u = 0; u < VEC_UNROLL; u++)
      {
        v[u].store(out+iout);
        iout += VEC_SIZE;
      }
    }
#endif
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
