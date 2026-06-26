"""
Calculate the cal-pulser trigger efficiency of RNO-G runs from the header files alone.

The cal pulser fires once per second, locked to the PPS, so over a run of duration T
seconds we expect T cal pulses. The trigger efficiency is the fraction of those that
were actually recorded:

    efficiency = n_calpulser / T .

Cal-pulser events are tagged as in plot_calpulser_amplitudes.py: an event is a cal
pulser when it triggers within ``--threshold`` seconds of the last PPS, i.e.
(t_trigger - t_pps) < threshold, where the time since the last PPS is derived from the
sysclk counter:
    dt = ((sysclk - sysclk_last_pps) mod 2**32) / f ,
with f (ticks/s) taken from the two most recent PPS sysclks (~100 MHz).

Runs from several stations may be passed at once; results are grouped by station and a
plot of trigger efficiency vs time is produced with one line per station.

Everything is read through the mattak.Dataset header interface (skip_incomplete=False);
only headers.root is needed.

Usage
-----
    python calpulser_trigger_efficiency.py /data/.../station23/run1 /data/.../station21/run2

Each positional argument may be a run directory or any .root file inside it.
"""
import argparse
import logging
import os
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import mattak.Dataset


def run_efficiency(rundir, threshold):
    """
    Cal-pulser trigger efficiency of a single run, read from its header file.

    Returns a dict with: station, run, n_events, n_cal (events within ``threshold`` of
    the PPS), duration (s), efficiency (n_cal / duration, or None if the duration is
    unavailable) and time (first-event trigger time, Unix seconds, for the time axis).
    """
    dataset = mattak.Dataset.Dataset(
        station=0, run=0, data_path=rundir, skip_incomplete=False)

    # warn if the run config did not enable the cal pulser (enable_cal=0)
    try:
        if not dataset.is_calibration_run():
            logging.warning("%s is not a cal-pulser run (enable_cal=0 in run config); "
                            "selected events may be noise/RF near the PPS", rundir)
    except ValueError as e:
        logging.warning("Could not determine cal-pulser status for %s: %s", rundir, e)

    dataset.setEntries((0, dataset.N()))
    infos = dataset.eventInfo()

    sysclk = np.array([i.sysclk for i in infos], dtype=np.int64)
    pps0 = np.array([i.sysclkLastPPS[0] for i in infos], dtype=np.int64)
    pps1 = np.array([i.sysclkLastPPS[1] for i in infos], dtype=np.int64)

    # ticks per second from the two most recent PPS sysclks (uint32, handle wrap)
    f = (pps0 - pps1) % (2 ** 32)
    # time since the last PPS (handle uint32 wrap-around)
    dt = ((sysclk - pps0) % (2 ** 32)) / f

    n_cal = int(np.sum(dt < threshold))
    duration = dataset.duration()
    efficiency = n_cal / duration if duration and duration > 0 else None

    return {
        "station": int(dataset.station),
        "run": int(dataset.run),
        "n_events": int(dt.size),
        "n_cal": n_cal,
        "duration": duration,
        "efficiency": efficiency,
        "time": float(infos[0].triggerTime) if len(infos) else None,
    }


def resolve_rundir(path):
    """A path may be a run directory or any file inside it; return the directory."""
    return path if os.path.isdir(path) else os.path.dirname(path) or "."


def print_table(results):
    """Print a per-run efficiency table (grouped by station) plus a combined total."""
    header = f"{'station/run':>14}  {'events':>8}  {'cal':>7}  {'dur [s]':>9}  {'eff':>7}"
    print(header)
    print("-" * len(header))

    tot_cal = 0
    tot_dur = 0.0
    # group by station, then order each station's runs by time
    for station in sorted({r["station"] for r in results}):
        runs = sorted((r for r in results if r["station"] == station),
                      key=lambda r: (r["time"] is None, r["time"]))
        for r in runs:
            eff = f"{r['efficiency']:.3f}" if r["efficiency"] is not None else "n/a"
            dur = f"{r['duration']:.0f}" if r["duration"] else "n/a"
            print(f"{'s%d/r%d' % (r['station'], r['run']):>14}  {r['n_events']:>8}  "
                  f"{r['n_cal']:>7}  {dur:>9}  {eff:>7}")
            tot_cal += r["n_cal"]
            if r["duration"]:
                tot_dur += r["duration"]

    print("-" * len(header))
    tot_eff = f"{tot_cal / tot_dur:.3f}" if tot_dur > 0 else "n/a"
    print(f"{'TOTAL':>14}  {'':>8}  {tot_cal:>7}  {tot_dur:>9.0f}  {tot_eff:>7}")


def plot_efficiency_vs_time(results, output):
    """Plot cal-pulser trigger efficiency vs time, with one line per station."""
    by_station = defaultdict(list)
    for r in results:
        if r["efficiency"] is None or r["time"] is None:
            continue
        by_station[r["station"]].append((r["time"], r["efficiency"]))

    if not by_station:
        logging.warning("No runs with both a valid efficiency and time; skipping plot.")
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    for station in sorted(by_station):
        pts = sorted(by_station[station])
        times = [datetime.fromtimestamp(t, tz=timezone.utc) for t, _ in pts]
        effs = [e for _, e in pts]
        ax.plot(times, effs, marker="o", label=f"{station}")

    ax.set_xlabel("Time (UTC)")
    ax.set_ylabel("Cal-pulser trigger efficiency")
    ax.set_title("Cal-pulser trigger efficiency vs time")
    ax.grid()
    ax.legend(title="Station", ncols=3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H:%M"))
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(output)
    print(f"Saved figure to {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+",
                        help="run directory or any .root file inside it (one per run)")
    parser.add_argument("--threshold", type=float, default=3e-6,
                        help="max (t_trigger - t_pps) in seconds to tag a cal pulser (default: 3e-6)")
    parser.add_argument("--output", type=str, default="calpulser_efficiency_vs_time.png",
                        help="output figure path (default: calpulser_efficiency_vs_time.png)")
    args = parser.parse_args()

    results = [run_efficiency(resolve_rundir(p), args.threshold) for p in args.paths]
    print_table(results)
    plot_efficiency_vs_time(results, args.output)


if __name__ == "__main__":
    main()
