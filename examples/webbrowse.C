// web browser example
// This supports only one user at a time 

R__LOAD_LIBRARY(libRootFftwWrapper.so) 
#include "FFTtools.h" 

mattak::Dataset * d = 0; 
mattak::Dataset * d2 = 0; 
TCanvas * cweb = 0;
TCanvas * cthresh = 0;
TCanvas * crates = 0;
TPad * pinfo = 0; 
TPad * pwf = 0; 
int current; 
THttpServer * srv = 0; 
TPaveText * txt = 0; 

bool fft = false;  
double * x = 0; 
std::complex<double> *X = 0 ;
double * mag = 0; 

int mode(std::string m) 
{
  if (m == "fft")
  {
    fft = true; 
    return 1; 
  }
  else fft = false;
  return 0; 
}

int setsize(int w, int h) 
{
  cweb->SetCanvasSize(w,h); 
  cweb->Update(); 
  cthresh->SetCanvasSize(w,h); 
  cthresh->Update(); 
  crates->SetCanvasSize(w,h*0.5); 
  crates->Update(); 
  return 0; 
}

int go(int i) 
{
  if (d->N() <= 0) return -1; 

  current= i; 
  if (current< 0) current= 0;
  if (current>= d->N()) current = d->N()-1; 
  d->setEntry(current); 
  cweb->cd(); 
  pinfo->cd(); 
  pinfo->Clear(); 
  if (txt) delete txt; 
  txt = new TPaveText(0.01,0.01,0.99,0.99); 
  txt->AddText(Form("S %d, R: %d, ev: %d",d->header()->station_number, d->header()->run_number, d->header()->event_number)); 
  TTimeStamp trigtime(int(d->header()->trigger_time), int(1e9*(d->header()->trigger_time-int(d->header()->trigger_time))));
  TTimeStamp readtime(int(d->header()->readout_time), int(1e9*(d->header()->readout_time-int(d->header()->readout_time))));
  TString str_trig = trigtime.AsString();
  TString str_read = readtime.AsString();
  txt->AddText(Form("Trigger time: %s, Readout time: %s", str_trig.Data(), str_read.Data())); 
  txt->AddText(Form("Force? %d, LT? %d, PPS? %d, RAD? %s", d->header()->trigger_info.force_trigger, d->header()->trigger_info.lt_trigger, d->header()->trigger_info.pps_trigger, 
        d->header()->trigger_info.radiant_trigger ?   
        ( d->header()->trigger_info.which_radiant_trigger == 0 ? "RF0" : d->header()->trigger_info.which_radiant_trigger == 1 ?"RF1" : "RFX") : "no"));
  txt->Draw(); 
  pwf->Clear(); 


  if (fft) 
  {
    //we'll be a bit sleazy here
    pwf->Divide(6,4,0.001,0.001); 

    for (int chan = 0; chan < 24; chan++) 
    {


      TGraph * g = d->raw()->makeGraph(chan); 
      TGraph * P = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(g); 
      TString gname = g->GetName();
      gname[0]='P'; 
      P->GetXaxis()->SetTitle("Freq [MHz]"); 
      P->GetXaxis()->SetRangeUser(0,1599); 
      P->GetYaxis()->SetRangeUser(-20,60); 
      P->GetYaxis()->SetTitle("Power [dBArb]"); 
      P->SetTitle(""); 
      P->SetTitle(g->GetTitle()); 
      P->SetName(gname); 
      P->SetLineColor(kRed+2); 
      P->SetBit(TGraph::kCanDelete);
      P->SetBit(TGraph::kIsSortedX); 
      P->SetBit(TGraph::kNotEditable); 
      TVirtualPad * p = pwf->cd(chan+1); 
      p->SetGridx();
      p->SetGridy();
      P->Draw("al"); 
    }
  }
  else 
  {
    d->raw()->drawWaveforms(mattak::WaveformPlotOptions(),pwf);
  }

  pwf->Update(); 
  gSystem->ProcessEvents(); 
  srv->ProcessRequests(); 
  return current; 
}


void draw_ds()
{
  TTree *st = d->daqStatusTree(); 
  cthresh->cd(); 

  TMultiGraph *grad_th = new TMultiGraph; 
  grad_th->SetTitle("RADIANT Thresholds (surface channels only); readout time ; threshold [V]");
  grad_th->SetBit(TObject::kCanDelete,true); 

  TMultiGraph *grad_sc = new TMultiGraph; 
  grad_sc->SetTitle("RADIANT Scalers (surface channels only); readout time ; rate [Hz]");
  grad_sc->SetBit(TObject::kCanDelete,true); 

  for (int ch = 12; ch <=20; ch++) 
  {
    int N = st->Draw(Form("radiant_thresholds[%d]/16777215*2.5:radiant_scalers[%d]:readout_time_radiant", ch, ch),"","goff"); 
    TGraph * gth = new TGraph(N, st->GetV3(), st->GetV1()); 
    gth->SetTitle(Form("RAD ch%d; readout time; threshold [V]", ch));
    gth->SetName(Form("radch%d", ch));
    gth->GetXaxis()->SetTimeDisplay(1); 
    TGraph * gsc = new TGraph(N, st->GetV3(), st->GetV2()); 
    gsc->SetTitle(Form("RAD ch%d; readout time; rate [Hz]", ch));
    gth->SetName(Form("radch%d", ch));
    gsc->GetXaxis()->SetTimeDisplay(1); 
    grad_th->Add(gth);
    grad_sc->Add(gsc); 
  }

  TMultiGraph * gflo_th = new TMultiGraph; 
  gflo_th->SetTitle("FLOWER Thresholds; readout time; threshold [arb]");
  gflo_th->SetBit(TObject::kCanDelete,true); 

  TMultiGraph * gflo_sc = new TMultiGraph; 
  gflo_sc->SetTitle("FLOWER Scalers; readout time; rate [Hz]");
  gflo_sc->SetBit(TObject::kCanDelete,true); 

  int N = st->Draw("lt_1Hz_scalers.trig_coinc:lt_1Hz_scalers.servo_coinc:readout_time_lt:lt_1Hz_gated_scalers.trig_coinc:lt_1Hz_gated_scalers.servo_coinc","","goff"); 

  TGraph * gcoinc = new TGraph(N, st->GetV3(), st->GetV1()); 
  gcoinc->SetTitle("Trig Coinc.; readout time; scaler [Hz]"); 
  gcoinc->SetName("lttrigcoinc");
  gflo_sc->Add(gcoinc); 

  TGraph * gcoinc_servo = new TGraph(N, st->GetV3(), st->GetV2()); 
  gcoinc_servo->SetTitle("Servo Coinc.; readout time; scaler [Hz]"); 
  gcoinc_servo->SetName("ltservocoinc");
  gflo_sc->Add(gcoinc_servo); 

  TGraph * gcoinc_gated = new TGraph(N, st->GetV3(), st->GetV4()); 
  gcoinc_gated->SetTitle("Trig Coinc. (gated); readout time; scaler [Hz]"); 
  gcoinc_gated->SetName("lttrigcoincgate");
  gflo_sc->Add(gcoinc_gated); 

  TGraph * gcoinc_servo_gated = new TGraph(N, st->GetV3(), st->GetVal(4)); 
  gcoinc_servo_gated->SetTitle("Servo Coinc. (gated); readout time; scaler [Hz]"); 
  gcoinc_servo_gated->SetName("ltservocoincgate");
  gflo_sc->Add(gcoinc_servo_gated); 


  for (int ch = 0; ch < 4; ch++)
  {
    N = st->Draw(Form("lt_1Hz_scalers.trig_per_chan[%d]:lt_1Hz_scalers.servo_per_chan[%d]:lt_trigger_thresholds[%d]:lt_servo_thresholds[%d]:readout_time_lt", ch,ch,ch,ch),"","goff"); 

    TGraph * gth = new TGraph(N, st->GetVal(4), st->GetV3()); 
    gth->SetTitle(Form("LT ch%d trig.; readout time; threshold [arb]", ch));
    gth->SetName(Form("ltch%dtrig",ch));
    gflo_th->Add(gth);

    TGraph * gth_servo = new TGraph(N, st->GetVal(4), st->GetV4()); 
    gth_servo->SetTitle(Form("LT ch%d servo.; readout time; threshold [arb]", ch));
    gth_servo->SetName(Form("ltch%dservo",ch));
    gflo_th->Add(gth_servo);

    TGraph * gsc = new TGraph(N, st->GetVal(4), st->GetV1()); 
    gsc->SetTitle(Form("LT ch%d trig; readout time; Scaler [Hz]", ch));
    gsc->SetName(Form("ltch%dtrig",ch));
    gflo_sc->Add(gsc); 

    TGraph * gsc_servo = new TGraph(N, st->GetVal(4), st->GetV2()); 
    gsc_servo->SetTitle(Form("LT ch%d servo; readout time; Scaler [Hz]", ch));
    gsc_servo->SetName(Form("ltch%dservo",ch));
    gflo_sc->Add(gsc_servo); 
  }




  cthresh->Clear(); 
  cthresh->Divide(2,2); 
  cthresh->cd(1); 
  grad_th->GetXaxis()->SetTimeDisplay(1);
  grad_th->GetXaxis()->SetTimeOffset(0,"GMT");
  gStyle->SetPalette(kRainBow); 
  grad_th->Draw("a pmc plc"); 
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->BuildLegend();

  cthresh->cd(2); 
  grad_sc->GetXaxis()->SetTimeDisplay(1);
  grad_sc->GetXaxis()->SetTimeOffset(0,"GMT");
  grad_sc->Draw("a pmc plc"); 

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->BuildLegend();
  cthresh->cd(3); 
  gflo_th->GetXaxis()->SetTimeDisplay(1);
  gflo_th->Draw("a pmc plc"); 
  gflo_th->GetXaxis()->SetTimeOffset(0,"GMT");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->BuildLegend();
  cthresh->cd(4); 
  gflo_sc->Draw("a pmc plc"); 
  gflo_sc->GetXaxis()->SetTimeDisplay(1);
  gflo_sc->GetXaxis()->SetTimeOffset(0,"GMT");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->BuildLegend();
  gPad->Update(); 


}

void draw_rates() 
{
  TTree *hd = d2->headTree(); 
  crates->Clear(); 
  crates->cd(); 
  //we need to figure out start time, end time and binning
 
  double min = hd->GetMinimum("trigger_time"); 
  double max = hd->GetMaximum("trigger_time"); 
  //we want approximately 10 second bins, rounded to the nearest minute 
  int nbins = ceil((max-min)/10); 

  TH1 * total = new TH1F("total","Trigger Rate;time; Rate [Hz]",nbins,min,max);
  total->SetLineColor(kBlack);
  total->SetMarkerColor(kBlack);
  total->SetLineWidth(2); 
  total->GetXaxis()->SetTimeDisplay(1); 
  total->GetXaxis()->SetTimeOffset(0,"GMT");
  TH1 * force = new TH1F("force","Force",nbins,min,max);
  force->SetLineColor(kRed-2); 
  force->SetMarkerColor(kRed-2); 
  TH1 * rad = new TH1F("radiant","Radiant",nbins,min,max);
  rad->SetLineColor(kBlue-2); 
  rad->SetMarkerColor(kBlue-2); 
  TH1 * lt = new TH1F("lt","LT",nbins,min,max);
  lt->SetLineColor(kMagenta+2); 
  lt->SetMarkerColor(kMagenta+2); 
  TH1 * pps = new TH1F("pps","PPS",nbins,min,max);
  pps->SetLineColor(kGray); 
  pps->SetMarkerColor(kGray); 

  total->SetStats(false); 
  hd->Draw("trigger_time >> total","","goff"); 
  total->Scale(1./total->GetBinWidth(1));
  total->SetBit(TObject::kCanDelete);
  total->Draw("L"); 
  hd->Draw("trigger_time >> force","trigger_info.force_trigger","goff"); 
  force->Scale(1./total->GetBinWidth(1));
  force->SetBit(TObject::kCanDelete);
  force->Draw("Lsame");
  hd->Draw("trigger_time >> lt","trigger_info.lt_trigger","goff"); 
  lt->Scale(1./total->GetBinWidth(1));
  lt->SetBit(TObject::kCanDelete);
  lt->Draw("Lsame"); 
  hd->Draw("trigger_time >> pps","trigger_info.pps_trigger","goff"); 
  pps->Scale(1./total->GetBinWidth(1));
  pps->SetBit(TObject::kCanDelete);
  pps->Draw("Lsame");
  hd->Draw("trigger_time >> radiant","trigger_info.radiant_trigger","goff"); 
  rad->Scale(1./total->GetBinWidth(1));
  rad->SetBit(TObject::kCanDelete);
  rad->Draw("Lsame"); 
  gPad->BuildLegend(); 
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->Update(); 
}


int set_run(int station, int run)
{
  d->loadRun(station,run); 
  d2->loadRun(station,run,false); 
  draw_ds();
  draw_rates(); 



  return go(0);
}


int  first() 
{
  return go (0);
}

int last() 
{
  return go(d->N()-1); 
}

int next() 
{
  return go(current+1); 
}

int previous() 
{
  return go(current-1); 
}


void init(int port, const char * data_dir) 
{
  d = new mattak::Dataset(data_dir); 
  d2 = new mattak::Dataset(data_dir); 
  cthresh = new TCanvas("cthresh","cthresh", 2400,1600); 
  crates = new TCanvas("crates","crates", 2400,1600); 
  cweb = new TCanvas("cweb","cweb", 2400,1600); 
  pinfo = new TPad("pinfo","pinfo",0.01,0.81, 0.99,0.99);
  pwf = new TPad("pwf","pwf",0.01,0.01, 0.99,0.79);

  pinfo->Draw(); 
  pwf->Draw(); 

  srv = new THttpServer (Form("http:%d;rw", port)); 
  srv->RegisterCommand("/Get","go(%arg1%)"); 
  srv->RegisterCommand("/SetStationRun","set_run(%arg1%,%arg2%)"); 
  srv->RegisterCommand("/Resize","setsize(%arg1%,%arg2%)"); 
  srv->RegisterCommand("/SetMode","mode(\"%arg1%\")"); 
  srv->RegisterCommand("/First","first()"); 
  srv->RegisterCommand("/Reload","go(current)"); 
  srv->RegisterCommand("/Last","last()"); 
  srv->RegisterCommand("/Next","next()"); 
  srv->RegisterCommand("/Previous","previous()"); 
  srv->SetDefaultPage("examples/webbrowse.htm"); 
  srv->SetTimer(); 
}


int webbrowse(int station, int run, const char * data_dir = NULL, int port = 12345) 
{
  if (data_dir && !*data_dir) data_dir=0; 
  gROOT->SetBatch(true); 
  init(port, data_dir); 
  set_run(station,run); 
//  while (!gSystem->ProcessEvents()); 
  return 0; 
}
