mattak::Dataset * d = 0; 
TCanvas * cweb = 0;
TPad * pinfo = 0; 
TPad * pwf = 0; 
int current; 
THttpServer * srv = 0; 
TPaveText * txt = 0; 


int setsize(int w, int h) 
{
  cweb->SetCanvasSize(w,h); 
  cweb->Update(); 
  return 0; 
}

int go(int i) 
{

  current= i; 
  if (current< 0) current= 0;
  if (current>= d->N()) current = d->N()-1; 
  d->setEntry(current); 
  pinfo->cd(); 
  pinfo->Clear(); 
  if (txt) delete txt; 
  txt = new TPaveText(0.01,0.01,0.99,0.99); 
  txt->AddText(Form("S %d, R: %d, ev: %d",d->header()->station_number, d->header()->run_number, d->header()->event_number)); 
  TTimeStamp trigtime(int(d->header()->trigger_time), int(1e9*(d->header()->trigger_time-int(d->header()->trigger_time))));
  TTimeStamp readtime(int(d->header()->trigger_time), int(1e9*(d->header()->trigger_time-int(d->header()->trigger_time))));
  TString str_trig = trigtime.AsString();
  TString str_read = readtime.AsString();
  txt->AddText(Form("Trigger time: %s, Readout time: %s", str_trig.Data(), str_read.Data())); 
  txt->AddText(Form("Force? %d, LT? %d, PPS? %d, RAD? %s", d->header()->trigger_info.force_trigger, d->header()->trigger_info.lt_trigger, d->header()->trigger_info.pps_trigger, 
        d->header()->trigger_info.radiant_trigger ?   
        ( d->header()->trigger_info.which_radiant_trigger == 0 ? "RF0" : d->header()->trigger_info.which_radiant_trigger == 1 ?"RF1" : "RFX") : "no"));
  txt->Draw(); 
  pwf->Clear(); 
  d->raw()->drawWaveforms(mattak::WaveformPlotOptions(),pwf);
  pwf->Update(); 
  gSystem->ProcessEvents(); 
  srv->ProcessRequests(); 
  return current; 
}

void set_run(int station, int run)
{
  d->loadRun(station,run); 
  go(0); 
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
  cweb = new TCanvas("cweb","cweb", 2400,1600); 
  pinfo = new TPad("pinfo","pinfo",0.01,0.81, 0.99,0.99);
  pwf = new TPad("pwf","pwf",0.01,0.01, 0.99,0.79);

  pinfo->Draw(); 
  pwf->Draw(); 

  srv = new THttpServer (Form("http:%d;rw", port)); 
  srv->RegisterCommand("/Get","go(%arg1%)"); 
  srv->RegisterCommand("/SetStationRun","set_run(%arg1%,%arg2%)"); 
  srv->RegisterCommand("/Resize","setsize(%arg1%,%arg2%)"); 
  srv->RegisterCommand("/First","first()"); 
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
