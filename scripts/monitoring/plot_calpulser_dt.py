"""
Plot the distribution of dt = (t_trigger - t_pps), the time since the last PPS, over
all events of all given runs, to help optimize the cal-pulser time threshold.

dt is derived from the sysclk counter via the mattak.Dataset header interface
(skip_incomplete=False), exactly as in plot_calpulser_amplitudes.py:
    dt = ((sysclk - sysclk_last_pps) mod 2**32) / f ,
with f (ticks/s) from the two most recent PPS sysclks (~100 MHz).

Cal-pulser events pile up at dt ~ 0 (locked to the PPS); random/RF triggers are
uniform in [0, 1) s. The left panel zooms into small dt to expose the cal-pulser
peak; the right panel shows the full second on a log scale.

Usage
-----
    python plot_calpulser_dt.py /data/.../station23/run1 /data/.../station23/run2

Each positional argument may be a run directory or any .root file inside it.
"""
import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

import mattak.Dataset


# candidate thresholds (seconds) summarized on the terminal to help pick a cut
CANDIDATE_THRESHOLDS = [1e-6, 2e-6, 3e-6, 5e-6, 1e-5, 2e-5]


def run_cal_label(dataset):
    """ Cal-pulser "channel/type" label from the run config, or None if not a cal run. """
    if not dataset.get_config("calib", "enable_cal"):
        return None
    parts = [dataset.get_config("calib", key) for key in ("channel", "type")]
    parts = [p for p in parts if p not in (None, "none")]
    return "/".join(parts) if parts else None


def event_dt(rundir):
    """
    Return (dt, station, run, cal_label): dt = time since last PPS (seconds) for every
    event in the run, read through the mattak.Dataset header interface, plus the
    cal-pulser channel/type label (e.g. "coax/pulser", or None if not a cal run).
    """
    dataset = mattak.Dataset.Dataset(
        station=0, run=0, data_path=rundir, skip_incomplete=False)
    dataset.setEntries((0, dataset.N()))
    infos = dataset.eventInfo()

    sysclk = np.array([i.sysclk for i in infos], dtype=np.int64)
    pps0 = np.array([i.sysclkLastPPS[0] for i in infos], dtype=np.int64)
    pps1 = np.array([i.sysclkLastPPS[1] for i in infos], dtype=np.int64)

    # ticks per second from the two most recent PPS sysclks (uint32, handle wrap)
    f = (pps0 - pps1) % (2 ** 32)
    # time since the last PPS (handle uint32 wrap-around)
    dt = ((sysclk - pps0) % (2 ** 32)) / f
    return dt, int(dataset.station), int(dataset.run), run_cal_label(dataset)


def resolve_rundir(path):
    """A path may be a run directory or any file inside it; return the directory."""
    return path if os.path.isdir(path) else os.path.dirname(path) or "."


def collect_dt(paths):
    """Concatenate dt over all runs; return (dt_all, runs, cal_labels).

    ``cal_labels`` maps station -> set of cal-pulser channel/type labels seen.
    """
    dts = []
    runs = []
    cal_labels = {}
    for path in paths:
        rundir = resolve_rundir(path)
        dt, station, run, cal_label = event_dt(rundir)
        dts.append(dt)
        runs.append((station, run))
        cal_labels.setdefault(station, set()).add(cal_label)
        print(f"station{station}/run{run}: {dt.size} events")
    return np.concatenate(dts), runs, cal_labels


def print_threshold_summary(dt_all):
    """Print, for each candidate threshold, the count/fraction of events below it."""
    n = dt_all.size
    print(f"\nThreshold summary ({n} events total):")
    print(f"  {'threshold':>12}  {'n below':>10}  {'fraction':>10}")
    for thr in CANDIDATE_THRESHOLDS:
        below = int(np.sum(dt_all < thr))
        print(f"  {thr:12.1e}  {below:10d}  {below / n:10.4f}")


def format_station_title(cal_labels):
    """Build a per-station title like "station23 (coax/pulser)" from a
    {station: set(cal-pulser labels)} mapping, joining multiple stations with ", "."""
    parts = []
    for station in sorted(cal_labels):
        labels = sorted(lbl for lbl in cal_labels[station] if lbl is not None)
        suffix = f" ({', '.join(labels)})" if labels else ""
        parts.append(f"station{station}{suffix}")
    return ", ".join(parts)


def plot_dt(dt_all, runs, output, zoom_us, bins, threshold, station_title=None):
    """Two-panel dt distribution: zoom near 0 (left) and full second log-scale (right)."""
    fig, (ax_zoom, ax_full) = plt.subplots(1, 2, figsize=(14, 5))

    # left: zoom into small dt, in microseconds
    dt_us = dt_all * 1e6
    ax_zoom.hist(dt_us, bins=bins, range=(0, zoom_us))
    ax_zoom.axvline(threshold * 1e6, color="r", ls="--",
                    label=f"threshold = {threshold * 1e6:g} us")
    ax_zoom.set_xlabel(r"$dt = t_\mathrm{trigger} - t_\mathrm{pps}$  ($\mu$s)")
    ax_zoom.set_ylabel("Events")
    ax_zoom.set_yscale("log")
    ax_zoom.set_title(f"Zoom: dt < {zoom_us:g} us")
    ax_zoom.legend()
    ax_zoom.grid()

    # right: full second, log-scale, to show the uniform random-trigger background
    ax_full.hist(dt_all, bins=bins, range=(0, 1))
    ax_full.axvline(threshold, color="r", ls="--")
    ax_full.set_xlabel(r"$dt = t_\mathrm{trigger} - t_\mathrm{pps}$  (s)")
    ax_full.set_ylabel("Events")
    ax_full.set_yscale("log")
    ax_full.set_title("Full second")
    ax_full.grid()

    tag = ", ".join(f"s{s}r{r}" for s, r in runs)
    suptitle = f"dt distribution ({dt_all.size} events; {tag})"
    if station_title and len(station_title) < 100:
        suptitle += f"\n{station_title}"
    fig.suptitle(suptitle)
    fig.tight_layout()
    fig.savefig(output)
    print(f"\nSaved figure to {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+",
                        help="run directory or any .root file inside it (one per run)")
    parser.add_argument("--zoom-us", type=float, default=20.0,
                        help="upper edge of the zoomed panel, in microseconds (default: 20)")
    parser.add_argument("--bins", type=int, default=100, help="number of histogram bins")
    parser.add_argument("--threshold", type=float, default=3e-6,
                        help="candidate cal-pulser threshold (s) drawn as a marker (default: 3e-6)")
    parser.add_argument("--output", type=str, default=None,
                        help="output figure path (default: calpulser_dt_<runs>.png)")
    args = parser.parse_args()

    dt_all, runs, cal_labels = collect_dt(args.paths)
    print_threshold_summary(dt_all)

    output = args.output
    if output is None:
        tag = "_".join(f"s{s}r{r}" for s, r in runs)
        output = f"calpulser_dt_{tag}.png"

    plot_dt(dt_all, runs, output, args.zoom_us, args.bins, args.threshold,
            station_title=format_station_title(cal_labels))


if __name__ == "__main__":
    main()
