#include "mattak/Waveforms.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
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
  this->radiant_sampling_rate = wf->radiant_sampling_rate;

  for (unsigned i = 0; i < mattak::k::num_radiant_channels; i++)
  {
    if(wf->digitizer_readout_delay[i]!=0) this->digitizer_readout_delay_ns[i]=float(wf->digitizer_readout_delay[i])*128./float(wf->radiant_sampling_rate)*1000;
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

  for (int ch = 0; ch < mattak::k::num_radiant_channels; ch++)
  {
    vc.apply(ch, buffer_length, wf.radiant_data[ch], hdr.trigger_info.radiant_info.start_windows[ch][0], radiant_data[ch], isOldFirmware);
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

  double rate = wf.radiant_sampling_rate / 1000.;
  for (int isamp = 0; isamp < wf.buffer_length; isamp++)
  {
    double x = ns ? isamp/rate : isamp;
    g->GetX()[isamp] = x;
    g->GetY()[isamp] = wf.radiant_data[chan][isamp];
  }


  g->GetXaxis()->SetTitle(ns ? "time [ns]" : "sample number");
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

  int nrows = 1;
  int ncols = 1;

  double left_margin =opt.left_margin < 0.02 ? 0.02 : opt.left_margin;
  double first_col_ratio = 1;

  if (nplots > 1)
  {

    nrows =opt.rows ?:
                nplots < 4 ? 1:
                nplots < 9 ? 2:
                nplots < 12? 3:
                4;

    ncols = ceil (nplots / (float(nrows)));

    if (!opt.share_xaxis && !opt.share_yaxis)
    {
      where->Divide(ncols,nrows,0.0001,0.0001);
    }
    else
    {
      where->cd();
      double eps = 0.0001;
      first_col_ratio = 1. /( 1- (left_margin - 0.02));
      double first_row_ratio = 1. /( 1- (0.09));

      double total_x = first_col_ratio + ncols-1;
      double total_y = first_row_ratio + nrows-1;
      double first_col_width = first_col_ratio / total_x;
      double first_row_width = first_row_ratio / total_y;
      double row_width = 1. / total_y;
      double col_width = 1. / total_x;

      for (int irow = 0; irow < nrows; irow++)
      {
        double y = (irow > 0 ) * (first_row_width + (irow-1) * row_width);
        double ynext =  first_row_width + (irow) * row_width;
        for (int icol = 0; icol < ncols; icol++)
        {
          double x = (icol > 0 ) * (first_col_width + (icol-1) * col_width);
          double xnext =  first_col_width + (icol) * col_width;

          TPad * p = new TPad(Form("%s_%d_%d", where->GetName(),irow,icol), Form("%d%d", irow,icol), x+eps,y+eps, xnext-eps,ynext-eps);
          p->SetNumber(1+icol+ ncols * irow);
          p->Draw();
        }
      }
      where->Modified();
    }

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


  int chan_counter = 0;
  for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++)
  {

    if ((1<<ichan) & ~opt.mask) continue;

    where->cd(ipad++);

//    TGraph * g  = new TGraph(wf.buffer_length);
    TGraph * g =graphImpl(wf, ichan, opt.ns);

    if (opt.titles_map.count(ichan)) g->SetTitle(opt.titles_map.at(ichan));

    if (!opt.show_title || opt.share_xaxis) g->SetTitle("");

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

    int min_sample = opt.min_sample;
    int max_sample = opt.max_sample;

    while (min_sample < 0) min_sample += g->GetN();
    while (max_sample < 0) max_sample += g->GetN();
    if (min_sample >= g->GetN()) min_sample = 0;
    if (max_sample >= g->GetN()) max_sample = g->GetN()-1;

    double min_t = g->GetX()[min_sample];
    double max_t = g->GetX()[max_sample];
    g->GetXaxis()->SetRangeUser(min_t,max_t);

    g->SetLineColor(opt.line_colors_map.count(ichan) ? opt.line_colors_map.at(ichan) : opt.line_color);
    g->SetLineWidth(opt.line_width);
    g->SetLineStyle(opt.line_style);

    int row = chan_counter / ncols;
    int col = chan_counter % ncols;
    double text_scale_factor = col == 0 ? 1. / first_col_ratio : 1 ;

    TString the_title = g->GetTitle();
    g->SetTitle(""); // we'll draw it ourselves...

    g->Draw("al");

    gPad->SetGridx();
    gPad->SetGridy();


    TLatex ltx;
    ltx.SetTextFont(42);
    ltx.SetTextSize(opt.title_size * text_scale_factor);


    if (!opt.show_title || opt.share_xaxis) gPad->SetTopMargin(0.01);
    else
    {
      ltx.SetTextAlign(22);
      ltx.DrawLatexNDC(0.5, 0.95, the_title.Data());
    }


    if (opt.annotations_map.count(ichan))
    {
      ltx.SetTextAlign(11);
      ltx.DrawLatexNDC(opt.share_yaxis && col > 0 ? 0.1 : 0.1 + opt.left_margin, 0.8, opt.annotations_map.at(ichan));
    }



    gPad->SetRightMargin(0.02);


    if (!opt.share_yaxis || col == 0)
    {
      gPad->SetLeftMargin(opt.left_margin);
      g->GetYaxis()->SetTitleSize(opt.ytitle_size * text_scale_factor);
      g->GetYaxis()->SetLabelSize(opt.ylabel_size * text_scale_factor);
      g->GetYaxis()->SetLabelOffset(opt.ylabel_offset / text_scale_factor);
      g->GetYaxis()->CenterTitle(opt.ytitle_center * text_scale_factor);
      g->GetYaxis()->SetTitleOffset(opt.ytitle_offset / text_scale_factor);
    }
    else
    {
      gPad->SetLeftMargin(0.02);
      g->GetYaxis()->SetTitle("");
      g->GetYaxis()->SetLabelSize(0);
    }


    if (!opt.share_xaxis || row == nrows - 1)
    {
      g->GetXaxis()->SetTitleSize(opt.xtitle_size * text_scale_factor);
      g->GetXaxis()->SetTitleOffset(opt.xtitle_offset / text_scale_factor);
      g->GetXaxis()->SetLabelOffset(opt.xlabel_offset / text_scale_factor);
      g->GetXaxis()->SetLabelSize(opt.xlabel_size * text_scale_factor);
      g->GetXaxis()->CenterTitle(opt.xtitle_center * text_scale_factor);
    }
    else
    {
      g->GetXaxis()->SetTitle("");
      g->GetXaxis()->SetLabelSize(0);
      gPad->SetBottomMargin(0.01);
    }

    g->GetYaxis()->SetNdivisions(opt.yndivisions);
    g->GetXaxis()->SetNdivisions(opt.xndivisions);

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
    chan_counter++;
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
