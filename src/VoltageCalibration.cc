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



#ifndef MATTAK_NOROOT

#include "mattak/Pedestals.h"

#include "TLinearFitter.h"
#include "TF1.h"
#include "TList.h"
#include "TString.h"

ClassImp(mattak::VoltageCalibration);

//ROOT STUFF





mattak::VoltageCalibration::VoltageCalibration(TTree * tree, const char * branch_name, double vref, int fit_order, double min, double max)
{
  setupFromTree(tree,branch_name,vref,fit_order,min,max);


}
//almost copy-pasted from raw file reader...
void mattak::VoltageCalibration::setupFromTree(TTree * tree, const char * branch_name, double vref, int fit_order, double min, double max)
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
  recalculateFits(fit_order, min, max, vref);
}

mattak::VoltageCalibration::VoltageCalibration(const char * bias_scan_file, double vref, int fit_order, double min, double max)
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
    setupFromTree(t,"pedestals", vref, fit_order, min, max);
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
  recalculateFits(fit_order, min, max, vref);
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


void mattak::VoltageCalibration::recalculateFits(int order, double min, double max, double vref, uint32_t mask, int turnover_threshold)
{
  fit_order = order < max_voltage_calibration_fit_order ? order : max_voltage_calibration_fit_order;
  order = fit_order;
  fit_min = min;
  this->turnover_threshold = turnover_threshold;
  fit_max = max;
  fit_vref = vref;

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

      turnover_index[ichan][i] = npoints;

      TF1 * fn = new TF1("polFit", formula[fit_order]);
      fn->SetParameters(getFitCoeffs(ichan,i));
      fn->SetRange(data_v[0],data_v[npoints-1]);

      // Sum up all the ADC(V) residuals
      if (turnover_index[ichan][i] == npoints_general)
      {
        if (ichan < 12)
        {
          icount_dac1++;
          for (int j = 0 ; j < npoints; j++)
          {
            double adc = data_adc[j];
            double v = data_v[j];
            double adcResid = adc - fn->Eval(v);
            sum_vol_dac1[j] = sum_vol_dac1[j] + v;
            sum_adcResid_dac1[j] = sum_adcResid_dac1[j] + adcResid;
          }
        }
        else
        {
          icount_dac2++;
          for (int j = 0 ; j < npoints; j++)
          {
            double adc = data_adc[j];
            double v = data_v[j];
            double adcResid = adc - fn->Eval(v);
            sum_vol_dac2[j] = sum_vol_dac2[j] + v;
            sum_adcResid_dac2[j] = sum_adcResid_dac2[j] + adcResid;
          }
        }
      }

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

  std::cout << "\nCalculating max deviation and chi squared..." << std::endl;
  for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++)
  {
    std::cout << "Channel " << ichan;
    if ( (mask & (1 << ichan))  == 0) continue;
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

      TF1 * fn = new TF1("polFit", formula[fit_order]);
      fn->SetParameters(getFitCoeffs(ichan,i));
      fn->SetRange(data_v[0],data_v[npoints-1]);

      // Calculate max deviation and chi squared
      for (int j = 0 ; j < npoints; j++)
      {
        double adc_adjusted;
        if (ichan < 12) { adc_adjusted = data_adc[j] - graph_residAve[0]->GetPointY(j); }
        else { adc_adjusted = data_adc[j] - graph_residAve[1]->GetPointY(j); }
        double v_meas = data_v[j];
        double v_pred = fn->GetX(adc_adjusted);
        double delta = fabs(v_meas-v_pred);
        if (isnan(delta)) continue;
        else
        {
          if (delta > fit_maxerr[ichan][i]) fit_maxerr[ichan][i] = delta;
          fit_chisq[ichan][i] += ( delta*delta/(0.001*0.001) );
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
                      Form("Bias Scan, Station %d, Channel %d, Time [%d-%d], VRef=%g ; Sample ; Vbias [V] ; adu", station_number, chan, start_time, end_time,fit_vref),
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
  g->GetXaxis()->SetTitle("VBias");
  g->GetYaxis()->SetTitle(resid ? "ADC - Predicted ADC" : "ADC");

  TF1 * fn = new TF1(Form("fsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time), formula[fit_order], fit_min, fit_max, TF1::EAddToList::kNo);
  fn->SetParameters(getFitCoeffs(chan,samp));
  fn->SetRange(data_v[0],data_v[npoints-1]);
  fn->SetLineColor(2);

  for (unsigned j = 0; j < turnover_index[chan][samp]; j++)
  {
    double adc_adjusted;
    if (chan < 12) { adc_adjusted = data_adc[j] - graph_residAve[0]->GetPointY(j); }
    else { adc_adjusted = data_adc[j] - graph_residAve[1]->GetPointY(j); }
    double v = data_v[j];
    if (resid) adc_adjusted -= fn->Eval(v);
    g->SetPoint(g->GetN(),v,adc_adjusted);
  }

  if (!resid)
  {
    g->GetListOfFunctions()->Add(fn);
    fn->SetParent(g);
    fn->Save(data_v[0],data_v[npoints-1],0,0,0,0);
  }

  return g;
}

TGraph * mattak::VoltageCalibration::makeAdjustedSampleGraph(int chan, int samp, bool resid) const
{
  int npoints = graph[chan][samp]->GetN();
  double *data_adc = graph[chan][samp]->GetY();
  double *data_v = graph[chan][samp]->GetX();

  TGraph *g = new TGraph();
  g->SetName(Form("gsample_adjusted_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time));
  g->SetTitle(Form("Station %d Ch %d sample %d [%d-%d], #chi^{2}= %g  %s", station_number, chan, samp, start_time, end_time, fit_chisq[chan][samp], resid ? "(residuals)" : ""));
  g->GetXaxis()->SetTitle("ADC");
  g->GetYaxis()->SetTitle(resid ? "VBias - Predicted VBias" : "VBias");

  TF1 * fn = new TF1(Form("fsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time), formula[fit_order], fit_min, fit_max, TF1::EAddToList::kNo);
  fn->SetParameters(getFitCoeffs(chan,samp));
  fn->SetRange(data_v[0],data_v[npoints-1]);

  for (unsigned j = 0; j < turnover_index[chan][samp]; j++)
  {
    double adc_adjusted;
    if (chan < 12) { adc_adjusted = data_adc[j] - graph_residAve[0]->GetPointY(j); }
    else { adc_adjusted = data_adc[j] - graph_residAve[1]->GetPointY(j); }
    double v = data_v[j];
    if (isnan(fn->GetX(adc_adjusted))) continue;
    else
    {
      if (resid) v -= fn->GetX(adc_adjusted);
      g->SetPoint(g->GetN(),adc_adjusted,v);
    }
  }

  return g;
}

TGraph * mattak::VoltageCalibration::makeOriginalSampleGraph(int chan, int samp) const
{
  double *data_adc = graph[chan][samp]->GetY();
  double *data_v = graph[chan][samp]->GetX();

  TGraph *g = new TGraph();
  g->SetName(Form("gsample_s%d_c%d_s%d_%d_%d", station_number, chan, samp, start_time, end_time));
  g->SetTitle(Form("Station %d Ch %d sample %d [%d-%d]", station_number, chan, samp, start_time, end_time));
  g->GetXaxis()->SetTitle("ADC");
  g->GetYaxis()->SetTitle("VBias");

  for (unsigned j = 0; j < turnover_index[chan][samp]; j++)
  {
    double adc = data_adc[j];
    double v = data_v[j];
    g->SetPoint(g->GetN(),adc,v);
  }

  return g;
}



#else
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

          ///
//SHARED  /// This function is not used or needs to be modified as the new fitting approach is applied!
          ///

double * mattak::applyVoltageCalibration (int N, const int16_t * in, double * out, int start_window,
                                   const double * packed_fit_params, int fit_order, double fit_min, double fit_max)

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

  int start_samp = start_window * mattak::k::radiant_window_size;
  int nwindows = N / mattak::k::radiant_window_size;

  int i = 0;
  for (int iwindow = 0; iwindow < nwindows; iwindow++)
  {
    //wrap around if isamp % 2048 == 0
    if (iwindow % mattak::k::radiant_windows_per_buffer == 0 )
    {
      iwindow -= mattak::k::radiant_windows_per_buffer; ;
    }

    int isamp = iwindow * mattak::k::radiant_window_size;


#ifndef MATTAK_VECTORIZE
    for (int k  = 0; k < mattak::k::radiant_window_size; k++)
    {
      const double * params = packed_fit_params + isamp * (fit_order+1);

      //Evaluate via horner's method
      double x = in[i];
      out[i] = evalPars(x, fit_order, params);
      isamp++;
      i++;
    }
#else

//vectorized version, on intel x86x64 anyway. Optimized for AVX2.
#define VEC_SIZE 4
#define VECD vec4d
#define VECI vec4q
#define VEC_INCR VECI(0,1,2,3)
#define VEC_N mattak::k::radiant_window_size / (VEC_SIZE * VEC_UNROLL)

    VECD v[VEC_UNROLL];
    VECI vin[VEC_UNROLL];
    VECD x[VEC_UNROLL];
    VECI idx[VEC_UNROLL];


    for ( k = 0; k < vec_N ; k++)
    {
      int iout = i; // for simplicity

      for (int u = 0; u < VEC_UNROLL; u++)
      {
        vin[u].load(in + i);
        x[u] = vin[u];
        idx[u]= VEC_INCR + isamp:
        v[u] = lookup < mattak::max_voltage_calibration_fit_order * mattak::k::num_lab4_samples > ( idx, packed_fit_params);
        idx[u] *= (fit_order+1);
        isamp += VEC_SIZE;
        i+= VEC_SIZE;
      }

      for (int j = fit_order-1; j>=0; j--)
      {
        for (int u = 0; u < VEC_UNROLL; u++)
        {
          idx[u]++;
          v[u] = v[u] * x[u] + lookup < mattak::max_voltage_calibration_fit_order * mattak::k::num_lab4_samples > ( idx, packed_fit_params);
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


//PYTHON BINDING

#ifdef MATTAK_NOROOT

namespace py = pybind11;

static py::array_t<double> apply_voltage_calibration(py::buffer in, int start_window, py::buffer packed_coeffs, int order, double min, double max)
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

  if (!mattak::applyVoltageCalibration(N, (int16_t*) in.ptr(), (double*) ret.ptr(), start_window, (double*) packed_coeffs.ptr(), order, min, max));
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
