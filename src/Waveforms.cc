#include "mattak/Waveforms.h" 
#include "TPad.h" 
#include "TCanvas.h" 
#include "TGraph.h" 
#include "TPaveText.h" 
#include <iostream> 


ClassImp(mattak::Waveforms); 
ClassImp(mattak::IWaveforms); 
ClassImp(mattak::CalibratedWaveforms); 


mattak::Waveforms::Waveforms(const rno_g_waveform_t * wf ) 
  : Waveforms() 
{

#ifndef LIBRNO_G_SUPPORT
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) wf; 
#else

  this->run_number = wf->run_number; 
  this->event_number = wf->event_number; 
  this->buffer_length = wf->radiant_nsamples; 
  this->station_number = wf->station; 
  for (unsigned i = 0; i < mattak::k::num_radiant_channels; i++) 
  {
    memcpy(this->radiant_data[i], wf->radiant_waveforms[i], sizeof(int16_t) * mattak::k::num_radiant_samples); 
  }

#endif
}


mattak::CalibratedWaveforms::CalibratedWaveforms(const Waveforms & wf, const Header & hdr,  const VoltageCalibration & vc, bool isOldFirmware) 
{
  run_number = wf.run_number; 
  event_number = wf.event_number; 
  station_number = wf.station_number; 
  buffer_length = wf.buffer_length; 

  if (hdr.run_number != wf.run_number && hdr.event_number != wf.run_number && hdr.station_number != wf.station_number) 
  {
    std::cerr << "WARNING: Possible event-header mismatch" << std::endl; 
  }

  for (int i = 0; i < mattak::k::num_radiant_channels; i++) 
  {
    vc.apply( i, buffer_length, wf.radiant_data[i], hdr.trigger_info.radiant_info.start_windows[i][0], radiant_data[i], isOldFirmware); 
  }
}

template <typename T> 
static const char * getYaxisLabel () { return "amplitude [arb]" ; }

template<>
const char * getYaxisLabel<mattak::Waveforms> () { return "amplitude [adu]"; }

template<>
const char * getYaxisLabel<mattak::CalibratedWaveforms> () { return "amplitude [V]"; }


template <typename T> 
TGraph* graphImpl(const T & wf, int chan, bool ns) 
{
  TGraph * g  = new TGraph(wf.buffer_length); 
  g->SetBit(TGraph::kNotEditable); 
  g->SetBit(TGraph::kIsSortedX); 

  for (int isamp = 0; isamp < wf.buffer_length; isamp++) 
  {
    double x = ns ? isamp/3.2 : isamp; 
    g->GetX()[isamp] = x; 
    g->GetY()[isamp] = wf.radiant_data[chan][isamp]; 
  }


  g->GetXaxis()->SetTitle(ns ? "ns" : "sample"); 
  g->GetXaxis()->SetLimits(g->GetX()[0], g->GetX()[wf.buffer_length-1]); 
  g->GetXaxis()->SetRangeUser(g->GetX()[0], g->GetX()[wf.buffer_length-1]); 
  g->GetYaxis()->SetTitle(getYaxisLabel<T>()); 

  g->SetName(Form("g_s%d_r%d_e%d_ch%d", wf.station_number, wf.run_number, wf.event_number, chan));
  g->SetTitle(Form("Station %d, Run %d, Event %d, Ch %d", wf.station_number, wf.run_number, wf.event_number, chan));

  return g; 
}

template <typename T>
static TVirtualPad * drawImpl(const T & wf, const mattak::WaveformPlotOptions & opt, TVirtualPad * where) 
{
  int nplots = __builtin_popcount(opt.mask); 
  if (!nplots) return nullptr; 

  if (!where ) 
  {
    where = new TCanvas(Form("c_s%d_r%d_ev%d", wf.station_number, wf.run_number, wf.event_number), Form("Station %d, Run %d, Event %d", wf.station_number, wf.run_number, wf.event_number), opt.width, opt.height); 
  }


  //figure out how to divide the canvases

  if (nplots > 1)
  {

    int nrows = nplots < 4 ? 1: 
                nplots < 9 ? 2: 
                nplots < 12? 3: 
                4; 

    int ncols = ceil (nplots / nrows); 

    where->Divide(ncols,nrows,0.001,0.001); 
  }




  double gmin = DBL_MAX; 
  double gmax = -DBL_MAX; 

  int ipad = nplots == 1 ? 0 : 1; 
  if (opt.global_scale) 
  {
    for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++) 
    {
      if ( (1 <<ichan) & ~opt.mask) continue; 
      for (int isamp = 0; isamp < wf.buffer_length; isamp++) 
      {
        if (wf.radiant_data[ichan][isamp] > gmax) 
          gmax = wf.radiant_data[ichan][isamp]; 

        if (wf.radiant_data[ichan][isamp] < gmin) 
          gmin = wf.radiant_data[ichan][isamp]; 
      }
    }
  }


  for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++) 
  {

    if ((1<<ichan) & ~opt.mask) continue; 

    where->cd(ipad++); 

//    TGraph * g  = new TGraph(wf.buffer_length); 
    TGraph * g =graphImpl(wf, ichan, opt.ns); 

    g->SetBit(TGraph::kCanDelete); 

    double this_min = DBL_MAX; 
    double this_max  = -DBL_MAX; 
    if (!opt.global_scale || opt.stats) 
    {
      for (int isamp = 0; isamp < wf.buffer_length; isamp++) 
      {
         if (g->GetY()[isamp] > this_max) this_max = g->GetY()[isamp]; 
         if (g->GetY()[isamp] < this_min) this_min = g->GetY()[isamp]; 
      }
    }

    if (!opt.autoscale) 
    {
      g->GetYaxis()->SetRangeUser(opt.min,opt.max); 
    }
    else
    {
      if (!opt.global_scale) 
      {
        gmin = this_min; 
        gmax = this_max; 
      }

      if (opt.symmetric_scale) 
      {
        double abs_max = gmax > -gmin ? gmax : -gmin; 
        abs_max *= opt.scale_fact; 
        g->GetYaxis()->SetRangeUser(-abs_max, abs_max); 
      }
      else
      {
        g->GetYaxis()->SetRangeUser(gmin*opt.scale_fact, gmax*opt.scale_fact); 
      }
    }

    g->SetLineColor(kAzure+2); 
    g->Draw("al"); 

    g->GetXaxis()->SetRangeUser(g->GetX()[0], g->GetX()[wf.buffer_length-1]); 
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetRightMargin(0.01); 

    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetLabelSize(0.04);

    g->GetYaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetLabelSize(0.04);
    g->GetYaxis()->SetTitleOffset(0.9);
    if (opt.stats) 
    {
      TPaveText * pt = new TPaveText(0.75,0.7,0.99,0.9,"NB NDC"); 
      pt->SetFillStyle(0); 
      pt->SetLineWidth(0); 
      pt->SetBit(TObject::kCanDelete); 
      pt->AddText(Form("#mu: %g",g->GetMean(2)));
      pt->AddText(Form("#sigma: %g",g->GetRMS(2)));
      pt->AddText(Form("V_{pp}: %g",this_max-this_min));
      pt->Draw(); 
    }
  }
  return where; 
}


TVirtualPad * mattak::Waveforms::drawWaveforms(const WaveformPlotOptions & opt, TVirtualPad * where)  const
{
  return drawImpl(*this, opt, where); 
}

TVirtualPad * mattak::CalibratedWaveforms::drawWaveforms(const WaveformPlotOptions & opt, TVirtualPad * where)  const
{
  return drawImpl(*this, opt, where); 
}


TGraph * mattak::Waveforms::makeGraph(int chan, bool ns)  const
{
  return graphImpl(*this, chan,ns);
}

TGraph * mattak::CalibratedWaveforms::makeGraph(int chan, bool ns)  const
{
  return graphImpl(*this, chan,ns);
}

