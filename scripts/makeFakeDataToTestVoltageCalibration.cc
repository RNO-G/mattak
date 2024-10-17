
/* This makes some fake data */
#include "TTimeStamp.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TFile.h"

#include "mattak/Waveforms.h"
#include "mattak/Header.h"

/*
  To compile:

  gcc makeFakeDataToTestVoltageCalibration.cc \
    -o makeFakeDataToTestVoltageCalibration.out -I $PWD/../install/include \
    -lmattak -L $PWD/../install/lib  $(root-config --libs) $(root-config --cflags) -Wall -lstdc++
*/

int main()
{

  int station = 12345;
  int run = 3;
  const char * data_dir = "./";

  TTimeStamp now;
  unsigned year, month, day;
  now.GetDate(true,0,&year,&month,&day);
  unsigned hr, min, sec;
  now.GetTime(true,0,&hr,&min,&sec);

  gSystem->mkdir(Form("%s/station_%d/runs/run_%05d",data_dir, station, run), true);

  TFile wf_file(Form("%s/station_%d/runs/run_%05d/waveforms.root", data_dir, station, run),"RECREATE");
  auto wf_tree = new TTree("waveforms","waveforms");
  auto wf = new mattak::Waveforms;
  wf_tree->Branch("waveforms", &wf);


  TFile hd_file(Form("%s/station_%d/runs/run_%05d/headers.root", data_dir, station, run),"RECREATE");
  auto hd_tree = new TTree("header","header");
  auto hd = new mattak::Header;
  hd_tree->Branch("header",&hd);


  gRandom->SetSeed(0);

  hd->run_number= run;
  hd->station_number = station;
  wf->run_number= run;
  wf->station_number = station;


  hd->buffer_length = 2048;
  wf->buffer_length = 2048;
  hd->pretrigger_samples = 100;
  hd->trigger_info.rf_trigger = true;
  hd->trigger_info.radiant_trigger = true;

  int sysclk_rate = 100e6; // whatever
  int sysclk_pps_offset = gRandom->Uniform(0, sysclk_rate);
  hd->sysclk = sysclk_pps_offset;
  double t = 0;
  double event_rate = 1;
  for (int x = -1000; x < 1000; x += 50)
  {
    // Create each trace with window 0 and 16 to test all samples
    for (int window = 0; window < 2; window++) {

      hd->event_number++;
      wf->event_number++;
      hd->trigger_number++;
      hd->readout_time = t+gRandom->Exp(0.1);
      hd->pps_num = int(t);

      for (int ichan = 0; ichan < mattak::k::num_radiant_channels; ichan++)
      {

        if (window == 0) {
          hd->trigger_info.radiant_info.start_windows[ichan][0] = 0;
          hd->trigger_info.radiant_info.start_windows[ichan][1] = 0;
        } else {
          hd->trigger_info.radiant_info.start_windows[ichan][0] = 16;
          hd->trigger_info.radiant_info.start_windows[ichan][1] = 16;
        }

        for (int isamp = 0; isamp < mattak::k::num_radiant_samples; isamp++)
        {
          wf->radiant_data[ichan][isamp] = x;
        }
      }

      hd_tree->Fill();
      wf_tree->Fill();
      double dt = gRandom->Exp(event_rate);
      hd->sysclk += dt * sysclk_rate;
      t += dt;

    }

  };

  wf_file.Write();
  hd_file.Write();
};
