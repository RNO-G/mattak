/**
 * Offline replay of the rno-g-acq readout trigger-rate cap (see check_write_event() in
 * rno-g-ice-software/src/rno-g-acq.c).
 *
 * The live DAQ tracks RADIANT's two RF triggers (RF0/RF1) and the LT board's trigger
 * independently: each source gets its own trailing window of trigger times, and once more
 * than max_trigger_rate * trigger_rate_window of them fall within the trailing
 * trigger_rate_window seconds, further events from that source get their readout skipped
 * until the rate drops again. Force/PPS triggers, and RADIANT triggers that can't be
 * confidently attributed to RF0 or RF1 (RFX), are never capped.
 *
 * This walks a run in acquisition order and reproduces that same per-source decision for
 * each event, using the header's readout_time (closest offline analog to the CLOCK_MONOTONIC
 * timestamp the live code takes right after reading the event out) rather than trigger_time.
 *
 * Note this doesn't reproduce the live code's RATE_WINDOW_MAX_ENTRIES=4096 memory bound
 * (a fixed-size ring buffer per source on the embedded side); this offline version just
 * uses a plain deque, so pick absurd rate/window combinations and you'll see a difference.
 *
 * Writes one entry per event to an output ROOT tree, and prints a per-run summary.
 */

#include <iostream>
#include <deque>
#include <stdlib.h>
#include <string.h>
#include "mattak/Dataset.h"
#include "TFile.h"
#include "TTree.h"

enum RateSource
{
  RF0 = 0,
  RF1 = 1,
  LT = 2,
  OTHER = 3 // force/pps triggers, or RFX (ambiguous RADIANT0/1); never rate-capped
};

static const char * source_names[] = {"RF0", "RF1", "LT", "OTHER"};

static RateSource classify(const mattak::Header * h)
{
  if (h->trigger_info.radiant_trigger)
  {
    if (h->trigger_info.which_radiant_trigger == 0) return RF0;
    if (h->trigger_info.which_radiant_trigger == 1) return RF1;
    return OTHER; // RFX
  }
  if (h->trigger_info.lt_trigger) return LT;
  return OTHER; // force/pps
}

// Mirrors check_write_event()'s sliding-window logic: append, prune anything older than
// window, then skip if more than max_rate*window entries remain.
struct RateTracker
{
  std::deque<double> times;

  bool check(double t, double window, double max_rate)
  {
    times.push_back(t);
    while (!times.empty() && t - times.front() > window) times.pop_front();

    int max_count = (int) (max_rate * window);
    if (max_count < 1) max_count = 1;

    return (int) times.size() > max_count; // true == would be skipped
  }
};

static void usage()
{
  std::cerr << "Usage: rno-g-label-rate-capped station run[:endrun] output.root [max_rate_hz=2] [window_s=60] [data_dir=$RNO_G_DATA]" << std::endl;
}

int main(int nargs, char ** args)
{
  if (nargs < 4)
  {
    usage();
    return 1;
  }

  int station = atoi(args[1]);

  int start_run = -1, end_run = -1;
  if (strchr(args[2], ':')) sscanf(args[2], "%d:%d", &start_run, &end_run);
  else start_run = end_run = atoi(args[2]);

  if (start_run < 0)
  {
    usage();
    return 1;
  }

  double max_rate = nargs > 4 ? atof(args[4]) : 2;
  double window = nargs > 5 ? atof(args[5]) : 60;

  mattak::DatasetOptions opt;
  opt.partial_skip_incomplete = false;
  if (nargs > 6) opt.base_data_dir = args[6];

  TFile of(args[3], "RECREATE");
  TTree * t = new TTree("rate_cap", "Trigger rate cap labels");

  int station_number = 0;
  unsigned run_number = 0;
  unsigned event_number = 0;
  double readout_time = 0;
  int source = 0;      // RateSource: 0=RF0, 1=RF1, 2=LT, 3=OTHER (never capped)
  int would_skip = 0;
  int current_count = 0;
  double current_rate = 0;

  t->Branch("station", &station_number);
  t->Branch("run", &run_number);
  t->Branch("event_number", &event_number);
  t->Branch("readout_time", &readout_time);
  t->Branch("source", &source);
  t->Branch("would_skip", &would_skip);
  t->Branch("current_count", &current_count);
  t->Branch("current_rate", &current_rate);

  for (int run = start_run; run <= end_run; run++)
  {
    mattak::Dataset d(station, run, opt);

    if (d.N() <= 0)
    {
      std::cout << "Station " << station << " Run " << run << ": could not find dataset or dataset empty" << std::endl;
      continue;
    }

    RateTracker trackers[3]; // RF0, RF1, LT; reset for each run
    int nskipped[4] = {0, 0, 0, 0};
    int ntotal[4] = {0, 0, 0, 0};

    for (int i = 0; i < d.N(); i++)
    {
      d.setEntry(i);
      mattak::Header * h = d.header();

      station_number = h->station_number;
      run_number = h->run_number;
      event_number = h->event_number;
      readout_time = h->readout_time;

      RateSource src = classify(h);
      source = (int) src;
      ntotal[src]++;

      if (src == OTHER)
      {
        would_skip = 0;
        current_count = 0;
        current_rate = 0;
      }
      else
      {
        would_skip = trackers[src].check(readout_time, window, max_rate) ? 1 : 0;
        current_count = (int) trackers[src].times.size();
        current_rate = window > 0 ? current_count / window : 0;
        if (would_skip) nskipped[src]++;
      }

      t->Fill();
    }

    std::cout << "Station " << station << " Run " << run << " (max_rate=" << max_rate
              << " Hz, window=" << window << " s):" << std::endl;
    for (int s = RF0; s <= LT; s++)
    {
      std::cout << "    " << source_names[s] << ": " << nskipped[s] << " / " << ntotal[s] << " would be skipped" << std::endl;
    }
    std::cout << "    " << source_names[OTHER] << " (never capped): " << ntotal[OTHER] << std::endl;
  }

  of.cd();
  t->Write();

  return 0;
}
