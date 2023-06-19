#include "mattak/VoltageCalibration.h"
#include "TFile.h"
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
  int i = order -1;
  while (i>=0)
  {
    ans = x * ans + p[i--];
  }
  return ans;
}

static double adcToVolt(double in_adc, int order, const double * par, const double * resid_volt, const double * resid_adc)
{
  double out_volt;
  double adc_array[mattak::npoints_interpolatedAveResid];
  double volt_array[mattak::npoints_interpolatedAveResid];

  for (int i = 0; i < mattak::npoints_interpolatedAveResid; i++)
  {
    volt_array[i] = resid_volt[i];
    adc_array[i] = evalPars(volt_array[i], order, par) + resid_adc[i];

    if (in_adc == adc_array[i]) out_volt = volt_array[i]; // Lucky if this happens...
  }

  // Most likely we will get out_volt from interpolation
  for (int i = 0; i < mattak::npoints_interpolatedAveResid-1; i++)
  {
    if (in_adc > adc_array[i] && in_adc < adc_array[i+1])
    out_volt = volt_array[i] + (in_adc - adc_array[i])*(volt_array[i+1] - volt_array[i])/(adc_array[i+1] - adc_array[i]);
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
#include "TF1.h"
#include "TList.h"
#include "TString.h"

ClassImp(mattak::VoltageCalibration);


mattak::VoltageCalibration::VoltageCalibration(TTree * tree, const char * branch_name, double vref, int fit_order, double min, double max, bool getMaxErrAndChi2)
{
  setupFromTree(tree,branch_name,vref,fit_order,min,max,getMaxErrAndChi2);
}

//almost copy-pasted from raw file reader...
void mattak::VoltageCalibration::setupFromTree(TTree * tree, const char * branch_name, double vref, int fit_order, double min, double max, bool getMaxErrAndChi2)
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
  recalculateFits(fit_order, min, max, vref, getMaxErrAndChi2);
}

mattak::VoltageCalibration::VoltageCalibration(const char * bias_scan_file, double vref, int fit_order, double min, double max, bool getMaxErrAndChi2)
{
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
    setupFromTree(t,"pedestals", vref, fit_order, min, max, getMaxErrAndChi2);
    return;
  }
#ifndef LIBRNO_G_SUPPORT

  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) bias_scan_file;
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
  recalculateFits(fit_order, min, max, vref, getMaxErrAndChi2);
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


void mattak::VoltageCalibration::recalculateFits(int order, double min, double max, double vref, bool getMaxErrAndChi2, uint32_t mask, int turnover_threshold)
{
  gErrorIgnoreLevel = kFatal;

  fit_order = order < max_voltage_calibration_fit_order ? order : max_voltage_calibration_fit_order;
  order = fit_order;
  fit_min = min;
  this->turnover_threshold = turnover_threshold;
  fit_max = max;
  fit_vref = vref;
  fit_getMaxErrAndChi2 = getMaxErrAndChi2;


  if ( min > max)
  {
    fit_min = max;
    fit_max = min;
  }

  int icount_dac1 = 0;
  int icount_dac2 = 0;
  std::vector<double> sum_vol_dac1;
  std::vector<double> sum_adcResid_dac1;
  std::vector<double> sum_vol_dac2;
  std::vector<double> sum_adcResid_dac2;
  sum_vol_dac1.resize(npoints_general);
  sum_adcResid_dac1.resize(npoints_general);
  sum_vol_dac2.resize(npoints_general);
  sum_adcResid_dac2.resize(npoints_general);
  for(int ipoint = 0; ipoint < npoints_general; ipoint++)
  {
    sum_vol_dac1[ipoint] = 0;
    sum_adcResid_dac1[ipoint] = 0;
    sum_vol_dac2[ipoint] = 0;
    sum_adcResid_dac2[ipoint] = 0;
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
             else //have too interpolate
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
            graph[ichan][i]->SetPoint(graph[ichan][i]->GetN(),v,adc);
            jmax = j;
          }
          last_adc = adc;

          if (!jmin) jmin = j;
        }
      }

      int npoints = graph[ichan][i]->GetN();
      double *data_adc = graph[ichan][i]->GetY();
      double *data_v = graph[ichan][i]->GetX();

      for (int i = 0; i < npoints; i++)
      {
        fit.AddPoint(&data_v[i], data_adc[i]);
      }

      if (vref)
      {
        fit.FixParameter(0, 0);
      }
      if (nzero > 1)
      {
        nbroken++;
      }
      else
      {
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

      // Sum up all the ADC(V) residuals
      if (npoints == npoints_general)
      {
        if (ichan < 12) icount_dac1++;
        else icount_dac2++;

        for (int j = 0 ; j < npoints; j++)
        {
          double adc = data_adc[j];
          double v = data_v[j];
          double adcResid = adc - evalPars(v, fit_order, &fit_coeffs[ichan][i * (order+1)]);
          if (ichan < 12)
          {
            sum_vol_dac1[j] = sum_vol_dac1[j] + v;
            sum_adcResid_dac1[j] = sum_adcResid_dac1[j] + adcResid;
          }
          else
          {
            sum_vol_dac2[j] = sum_vol_dac2[j] + v;
            sum_adcResid_dac2[j] = sum_adcResid_dac2[j] + adcResid;
          }
        }
      }
      else std::cout << "WARNING: Number of data points is not " << npoints_general << "!!! / Channel" << ichan << " Sample" << i << std::endl;

    }

    if (nbroken) printf("WARNING: Channel %d seems to have %d broken samples?\n", ichan, nbroken);
  }

  // Calculate the average ADC(V) residuals and make a TGraph for each DAC
  double average_vol_dac1[npoints_general];
  double average_adcResid_dac1[npoints_general];
  double average_vol_dac2[npoints_general];
  double average_adcResid_dac2[npoints_general];
  for (int ipoint = 0; ipoint < npoints_general; ipoint++)
  {
    average_vol_dac1[ipoint] = sum_vol_dac1[ipoint]/icount_dac1;
    average_adcResid_dac1[ipoint] = sum_adcResid_dac1[ipoint]/icount_dac1;
    average_vol_dac2[ipoint] = sum_vol_dac2[ipoint]/icount_dac2;
    average_adcResid_dac2[ipoint] = sum_adcResid_dac2[ipoint]/icount_dac2;
  }
  graph_residAve[0] = new TGraph(npoints_general, average_vol_dac1, average_adcResid_dac1);
  graph_residAve[1] = new TGraph(npoints_general, average_vol_dac2, average_adcResid_dac2);
  graph_residAve[0]->SetNameTitle("aveResid_dac1");
  graph_residAve[1]->SetNameTitle("aveResid_dac2");
  for (int i = 0; i < 2; i++)
  {
    graph_residAve[i]->GetXaxis()->SetTitle("VBias [Volt]");
    graph_residAve[i]->GetYaxis()->SetTitle("ADC Residual");
  }

  //
  // Interpolating the average residuals
  //
  for (int j = 0; j < 2; j++)
  {
    for (int i = 0; i < npoints_general; i++)
    {
      resid_volt[j][i*2] = graph_residAve[j]->GetPointX(i);
      resid_adc[j][i*2] = graph_residAve[j]->GetPointY(i);
    }
    for (int i = 0; i < npoints_general-1; i++)
    {
      resid_volt[j][i*2+1] = (resid_volt[j][i*2] + resid_volt[j][i*2+2])/2;
      resid_adc[j][i*2+1] = resid_adc[j][i*2] + (resid_volt[j][i*2+1] - resid_volt[j][i*2])*(resid_adc[j][i*2+2] - resid_adc[j][i*2])/(resid_volt[j][i*2+2] - resid_volt[j][i*2]);
    }
  }

  if (fit_getMaxErrAndChi2)
  {
    std::cout << "\nCalculating max deviation and chi squared..." << std::endl;
    int dacType;

    for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++)
    {
      std::cout << "Channel " << ichan;
      if ( (mask & (1 << ichan))  == 0) continue;

      if (ichan < 12) dacType = 0;
      else dacType = 1;

      for (int i = 0; i < mattak::k::num_lab4_samples; i++)
      {
        if (!(i%128))
        {
          printf(".");
          fflush(stdout); // Print a dot every 128 samples processed
        }

        fit_chisq[ichan][i] = 0;
        fit_maxerr[ichan][i] = 0;

        int npoints = graph[ichan][i]->GetN();
        double *data_adc = graph[ichan][i]->GetY();
        double *data_v = graph[ichan][i]->GetX();

        // Calculate max deviation and chi squared
        for (int j = 0 ; j < npoints; j++)
        {
          double adc = data_adc[j];
          double v_meas = data_v[j];
          double v_pred = adcToVolt(adc, fit_order, &fit_coeffs[ichan][i * (order+1)], &resid_volt[dacType][0], &resid_adc[dacType][0]);
          double delta = fabs(v_meas-v_pred);
          if (delta > fit_maxerr[ichan][i]) fit_maxerr[ichan][i] = delta;
          fit_chisq[ichan][i] += ( delta*delta/(0.002*0.002) );
        }
      }
    }
  }

}


TH2S * mattak::VoltageCalibration::makeHist(int chan) const
{
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
  int npoints = graph[chan][samp]->GetN();
  double *data_adc = graph[chan][samp]->GetY();
  double *data_v = graph[chan][samp]->GetX();

  TGraph *g = new TGraph();
  g->SetName(Form("gsample_inverse_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time));
  g->SetTitle(Form("Station %d Ch %d sample %d [%d-%d]  %s", station_number, chan, samp, start_time, end_time, resid ? "(residuals)" : ""));
  g->GetXaxis()->SetTitle("VBias [Volt]");
  g->GetYaxis()->SetTitle(resid ? "ADC Residual" : "ADC");

  TF1 * fn = new TF1(Form("fsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time), formula[fit_order], fit_min, fit_max, TF1::EAddToList::kNo);
  fn->SetParameters(getFitCoeffs(chan,samp));
  fn->SetRange(data_v[0],data_v[npoints-1]);
  fn->SetLineColor(2);

  int dacType;
  if (chan < 12) dacType = 0;
  else dacType = 1;

  if (npoints != npoints_general) std::cout << "WARNING: Number of data points in this sample is not " << npoints_general << "!!! No adjusted plot!!!" << std::endl;
  else
  {
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
  }

  return g;
}

TGraph * mattak::VoltageCalibration::makeSampleGraph(int chan, int samp, bool resid) const
{
  int npoints = graph[chan][samp]->GetN();
  double *data_adc = graph[chan][samp]->GetY();
  double *data_v = graph[chan][samp]->GetX();

  TGraph *g = new TGraph();
  g->SetName(Form("gsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time));
  g->SetTitle(Form("Station %d Ch %d sample %d [%d-%d], #chi^{2}= %g  %s", station_number, chan, samp, start_time, end_time, fit_chisq[chan][samp], resid ? "(residuals)" : ""));
  g->GetXaxis()->SetTitle("ADC");
  g->GetYaxis()->SetTitle(resid ? "(VBias - Predicted VBias) [Volt]" : "VBias [Volt]");

  int dacType;
  if (chan < 12) dacType = 0;
  else dacType = 1;

  for (int j = 0; j < npoints; j++)
  {
    double v = data_v[j];
    double adc = data_adc[j];

    if (resid) v -= adcToVolt(adc, fit_order, getFitCoeffs(chan,samp), getPackedAveResid_volt(chan), getPackedAveResid_adc(chan));
    g->SetPoint(g->GetN(),adc,v);
  }

  return g;
}

void mattak::VoltageCalibration::saveFitCoeffsInFile()
{
  const TString outFileName = TString::Format("fitCoeffs_s%d_%d-%d.root", station_number, getStartTime(), getEndTime());
  TFile f(outFileName, "RECREATE", "", ROOT::RCompressionSetting::EDefaults::kUseSmallest);

  TTree *general_tree = new TTree("general_tree", "general_tree");
  general_tree->Branch("fitOrder", &fit_order, "fitOrder/I");
  general_tree->Branch("stationNumber", &station_number, "stationNumber/I");
  //general_tree->Branch("fitMin", &fit_min, "fitMin/D");
  //general_tree->Branch("fitMax", &fit_max, "fitMax/D");
  general_tree->Branch("startTime", &start_time, "startTime/i");
  general_tree->Branch("endTime", &end_time, "endTime/i");
  general_tree->SetDirectory(&f);
  general_tree->Fill();
  general_tree->Write();

  TTree fitCoeffs_tree("coeffs_tree", "coeffs_tree");
  std::vector<float> coeff(fit_order+1);
  std::vector<float> *p_coeff = &coeff;
  fitCoeffs_tree.Branch("coeff", "std::vector<float>", &p_coeff);
  fitCoeffs_tree.SetDirectory(&f);

  for (int iChan = 0; iChan < mattak::k::num_radiant_channels; iChan++)
  {
    for (int iSamp = 0; iSamp < mattak::k::num_lab4_samples; iSamp++)
    {
      for (int iOrder = 0; iOrder <= fit_order; iOrder++)
      {
        coeff[iOrder] = getFitCoeff(iChan, iSamp, iOrder);
      }
      fitCoeffs_tree.Fill();
    }
  }

  fitCoeffs_tree.Write();

  getAveResidGraph_dac1()->Write();
  getAveResidGraph_dac2()->Write();

  std::cout << "\nAll fit coefficients saved in file: " << outFileName << std::endl;
  f.Close();
}

void mattak::VoltageCalibration::readFitCoeffsFromFile(const char * inFile)
{
  //
  // Get information about the bias scan and fit coefficients from the input root file
  //
  TFile *inputFile = new TFile(inFile);

  TTree *general_tree = (TTree*)inputFile->Get("general_tree");
  general_tree->SetBranchAddress("fitOrder", &fit_order);
  general_tree->SetBranchAddress("stationNumber", &station_number);
  general_tree->SetBranchAddress("startTime", &start_time);
  general_tree->SetBranchAddress("endTime", &end_time);
  general_tree->GetEntry(0);

  TTree *fitCoeffs_tree = (TTree*)inputFile->Get("coeffs_tree");
  std::vector<float> coeff(fit_order+1);
  std::vector<float> *p_coeff = &coeff;
  fitCoeffs_tree->SetBranchAddress("coeff", &p_coeff);

  for(int iChan = 0; iChan < mattak::k::num_radiant_channels; iChan++)
  {
    fit_coeffs[iChan].clear();
    fit_coeffs[iChan].resize((fit_order+1)*mattak::k::num_lab4_samples, 0);
    for(int iSamp = 0; iSamp < mattak::k::num_lab4_samples; iSamp++)
    {
      fitCoeffs_tree->GetEntry(iChan*mattak::k::num_lab4_samples+iSamp);
      for(int iOrder = 0; iOrder < fit_order+1; iOrder++)
      {
        fit_coeffs[iChan][iSamp*(fit_order+1)+iOrder] = coeff[iOrder];
      }
    }
  }

  // Average residuals for both types of DACs
  graph_residAve[0] = (TGraph*)inputFile->Get("aveResid_dac1");
  graph_residAve[1] = (TGraph*)inputFile->Get("aveResid_dac2");

  //
  // Interpolating the average residuals
  //
  for (int j = 0; j < 2; j++)
  {
    for (int i = 0; i < npoints_general; i++)
    {
      resid_volt[j][i*2] = graph_residAve[j]->GetPointX(i);
      resid_adc[j][i*2] = graph_residAve[j]->GetPointY(i);
    }
    for (int i = 0; i < npoints_general-1; i++)
    {
      resid_volt[j][i*2+1] = (resid_volt[j][i*2] + resid_volt[j][i*2+2])/2;
      resid_adc[j][i*2+1] = resid_adc[j][i*2] + (resid_volt[j][i*2+1] - resid_volt[j][i*2])*(resid_adc[j][i*2+2] - resid_adc[j][i*2])/(resid_volt[j][i*2+2] - resid_volt[j][i*2]);
    }
  }


  inputFile->Close();
}


#else
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

double * mattak::applyVoltageCalibration (int N, const int16_t * in, double * out, int start_window, bool is2ndBoard, bool isOldFirmware,
                            int fit_order, const double * packed_fit_params, const double * packed_aveResid_volt, const double * packed_aveResid_adc)

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

  int nwindows = N / mattak::k::radiant_window_size;

  int isamp;
  start_window = start_window % mattak::k::radiant_windows_per_buffer;
  int i = 0;
  int j = start_window;

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

    if (is2ndBoard) isamp += mattak::k::num_radiant_samples;

    j++;

#ifndef MATTAK_VECTORIZE
    for (int k = 0; k < mattak::k::radiant_window_size; k++)
    {
      const double * params = packed_fit_params + isamp * (fit_order+1);

      double adc = in[i];
      out[i] = adcToVolt(adc, fit_order, params, packed_aveResid_volt, packed_aveResid_adc);

      isamp++;
      i++;
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

static py::array_t<double> apply_voltage_calibration(py::buffer in, int start_window, bool is2ndBoard, bool isOldFirmware, int order, py::buffer packed_coeffs, py::buffer packed_aveResid_volt, py::buffer packed_aveResid_adc)
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

  if (!mattak::applyVoltageCalibration(N, (int16_t*) in.ptr(), (double*) ret.ptr(), start_window, is2ndBoard, isOldFirmware, order, (double*) packed_coeffs.ptr(), (double*) packed_aveResid_volt.ptr(), (double*) packed_aveResid_adc.ptr()))
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
