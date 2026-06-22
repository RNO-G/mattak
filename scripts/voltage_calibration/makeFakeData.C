
/* This makes some fake data */ 

void makeFakeData(int station = 12345, int run = 1, int nsecs = 3600, double event_rate = 1, double daq_status_interval = 10, double sensors_interval = 10) 
{

  const char * data_dir = getenv("RNO_G_DATA"); 

  if (!data_dir)
  {
    std::cerr << "Please define RNO_G_DATA" << std::endl; 
    return ; 
  }

  TTimeStamp now; 
  unsigned year, month, day; 
  now.GetDate(true,0,&year,&month,&day); 
  unsigned hr, min, sec; 
  now.GetTime(true,0,&hr,&min,&sec); 





  gSystem->mkdir(Form("%s/station_%d/runs/run_%05d",data_dir, station, run), true); 
  gSystem->mkdir(Form("%s/station_%d/sensors/%u/%02u/%02u", data_dir, station, year, month, day), true); 

  TFile wf_file(Form("%s/station_%d/runs/run_%05d/waveforms.root", data_dir, station, run),"RECREATE"); 
  auto wf_tree = new TTree("waveforms","waveforms"); 
  auto wf= new mattak::Waveforms;
  wf_tree->Branch("waveforms", &wf); 


  TFile hd_file(Form("%s/station_%d/runs/run_%05d/header.root", data_dir, station, run),"RECREATE"); 
  auto hd_tree = new TTree("header","header"); 
  auto hd = new mattak::Header;
  hd_tree->Branch("header",&hd); 

  TFile st_file(Form("%s/station_%d/runs/run_%05d/daqstatus.root", data_dir, station, run),"RECREATE"); 
  auto st_tree = new TTree("daqstatus","daqstatus"); 
  auto st = new mattak::DAQStatus; 
  st_tree->Branch("daqstatus",&st); 

  TFile hk_file(Form("%s/station_%d/sensors/%d/%02u/%02u/%02u%02u%02u.root", data_dir, station, year, month, day, hr, min, sec),"RECREATE"); 
  auto hk_tree = new TTree("sensors","sensors"); 
  auto hk = new mattak::Sensors; 
  hk_tree->Branch("sensors",&hk); 


  double t0 = now.AsDouble(); 

  double t = t0; 
  gRandom->SetSeed(0); 

  double last_sensors_t = 0;
  double last_status_t = 0;

  hd->run_number= run; 
  hd->station_number = station;
  wf->run_number= run; 
  wf->station_number = station;


  hd->buffer_length = 2048;
  wf->buffer_length = 2048;
  hd->pretrigger_samples = 100;
  hd->trigger_info.rf_trigger =true;
  hd->trigger_info.radiant_trigger =true;

  int sysclk_rate = 100e6; // whatever
  int sysclk_pps_offset = gRandom->Uniform(0,sysclk_rate); 
  hd->sysclk = sysclk_pps_offset; 

  while (t < t0 + nsecs) 
  {

    if (t > last_sensors_t + sensors_interval)
    {
      last_sensors_t = t; 
      hk->readout_time = t + gRandom->Gaus(0,0.2); 

      hk->PV_voltage =  TMath::Min(0., 48 * sin ( t / (3600*24) )) + gRandom->Gaus(0,0.1) ; // whatever
      hk->PV_current =  TMath::Min(0., 1000*sin ( t / (3600*24) )) ; + gRandom->Gaus(0,100); // whatever
      hk->battery_voltage = 24 + gRandom->Gaus(0,0.1); 
      hk->battery_current = 1000 + gRandom->Gaus(0,100); 


      hk->T_power_board = gRandom->Gaus(20,2);
      hk->T_case = gRandom->Gaus(20,2);
      hk->T_amps = gRandom->Gaus(20,2);
      hk->T_controller_board = gRandom->Gaus(20,2);

      hk->radiant_current = gRandom->Gaus(2000,200); 
      hk->lt_current = 0; 
      hk->sbc_current = gRandom->Gaus(200,10); 

      hk->downhole_currents[0] = gRandom->Gaus(200,10); 
      hk->downhole_currents[1] = gRandom->Gaus(200,10); 
      hk->downhole_currents[2] = gRandom->Gaus(200,10); 
      
      hk->amp_currents[0] = gRandom->Gaus(200,10); 
      hk->amp_currents[1] = gRandom->Gaus(200,10); 
      hk->amp_currents[2] = gRandom->Gaus(200,10); 
      hk->amp_currents[3] = gRandom->Gaus(200,10); 
      hk->amp_currents[4] = gRandom->Gaus(200,10); 
      hk->amp_currents[5] = gRandom->Gaus(200,10); 

      hk->disk_space_SD_card = 100; 
      hk->disk_space_internal = 1.1; 
      hk->free_mem = 123.45; 
      hk->loadavg[0] =1;
      hk->loadavg[1] =1;
      hk->loadavg[2] =1;
      hk->cpu_uptime = t-t0-1000;
      hk->mcu_uptime = t-t0-1000;

      hk->lte_rsrq = -5; 
      hk->lte_rsrp = -90; 
      hk->lte_rssi = -80; 
      hk->lte_snr = -5; 

      hk_tree->Fill(); 
    }

    if (t > last_status_t + daq_status_interval)
    {
      last_status_t = t; 
      st->readout_time = t+gRandom->Exp(0.1); 

      //TODO: fill rest 
      //
      st_tree->Fill(); 
    }


    hd->event_number++; 
    wf->event_number++; 
    hd->trigger_number++; 
    hd->readout_time = t+gRandom->Exp(0.1); 
    hd->pps_num = int(t); 

    //TODO: figure out last_pps and last_last_pps 
    hd->trigger_info.surface_trigger = gRandom->Rndm() < 0.1;


    for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++)
    {
      for (int isamp = 0; isamp < mattak::k::num_radiant_samples; isamp++)
      {
        wf->radiant_data[ichan][isamp] = gRandom->Gaus(2048,20); 
      }
    }

    hd_tree->Fill(); 
    wf_tree->Fill(); 
    double dt = gRandom->Exp(event_rate); 
    hd->sysclk += dt * sysclk_rate; 
    t+= dt; 
  }; 


  wf_file.Write();
  hd_file.Write(); 
  st_file.Write(); 
  hk_file.Write(); 


}; 

