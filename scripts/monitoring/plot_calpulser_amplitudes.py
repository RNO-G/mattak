"""
Select cal-pulser events from RNO-G runs and histogram their maximum amplitude.

For every (monitoring.root, run-directory) pair:
  1. Read the per-event ``max_abs_amplitude`` (24 channels) from monitoring.root.
  2. Read the run header via the ``mattak.Dataset`` interface (skip_incomplete=False).
  3. Tag an event as a cal pulser when it triggers within ``--threshold`` seconds
     of the last PPS, i.e. (t_trigger - t_pps) < threshold, where the time since
     the last PPS is derived from the sysclk counter:
         dt = ((sysclk - sysclk_last_pps) mod 2**32) / f ,
     with f (ticks/s) taken from the two most recent PPS sysclks (~100 MHz).
  4. Match cal-pulser events to the monitoring amplitudes by event number and
     histogram the max. abs. amplitude for the deep (in-ice) channels.

Usage
-----
    python plot_calpulser_amplitudes.py /data/.../station23/run1 /data/.../station23/run2

Each positional argument may be a run directory or either of the two files it
contains (monitoring.root / headers.root); both are resolved from it.
"""
import argparse
import logging
import os
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

import uproot
import ROOT
import mattak.backends.pyroot.mattakloader  # noqa: F401  (loads the EventSummary dictionary)

import mattak.Dataset


# Deep (in-ice) channels: the three power-string deep channels live in 0-11 and 21-23.
DEEP_CHANNELS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 22, 23]

# 12-bit ADC: codes 0..4095 span the 2.5 V dynamic range.
ADC_MAX_CODE = 4095


def run_cal_label(dataset):
    """ Cal-pulser "channel/type" label from the run config, or None if not a cal run. """
    if not dataset.get_config("calib", "enable_cal"):
        return None
    parts = [dataset.get_config("calib", key) for key in ("channel", "type")]
    parts = [p for p in parts if p not in (None, "none")]
    return "/".join(parts) if parts else None


def read_monitoring_summary(path):
    """
    Return (amplitudes, rms): two {event_number: np.ndarray[24]} dicts from a
    monitoring.root file, holding max_abs_amplitude and the per-channel RMS.
    """
    f = ROOT.TFile(path, "READ")
    event_tree = f.Get("events")

    amplitudes = {}
    rms = {}
    for entry in event_tree:
        summary = entry.EventSummary
        ev = int(summary.event_number)
        amplitudes[ev] = np.array(summary.max_abs_amplitude, dtype=float)
        rms[ev] = np.array(summary.rms, dtype=float)

    f.Close()
    return amplitudes, rms


def pedestal_near_rail(rundir):
    """
    Return the per-channel distance (in ADC counts) from the centered baseline to the
    nearest ADC rail, i.e. min(pedestal, 4095 - pedestal). The mattak raw waveforms are
    centered near 0, but the true baseline is not at mid-scale, so the rails are
    asymmetric; the nearer rail is where a growing signal clips first.

    Returns an array of shape (24,) or None if pedestal.root is not available.
    """
    path = os.path.join(rundir, "pedestal.root")
    if not os.path.exists(path):
        return None

    tree = uproot.open(path)["pedestals"]
    ped = tree["pedestals[24][4096]"].array(library="np")[0]  # (24, 4096) ADC codes
    mean_ped = ped.mean(axis=1)
    return np.minimum(mean_ped, ADC_MAX_CODE - mean_ped)


def calpulser_event_numbers(rundir, threshold):
    """
    Identify cal-pulser and force-trigger events using the mattak.Dataset interface
    (skip_incomplete=False). Cal-pulser events are those triggering within ``threshold``
    seconds of the last PPS.

    Returns (cal_evs, force_evs, station, run, duration, cal_label): the cal-pulser and
    force-trigger event-number sets, the station, run, run duration (seconds) and the
    cal-pulser channel/type label (e.g. "coax/pulser", or None if not a cal run).
    """
    dataset = mattak.Dataset.Dataset(
        station=0, run=0, data_path=rundir, backend="pyroot",
        skip_incomplete=False, voltage_calibration=None)

    # warn if the run config did not enable the cal pulser (enable_cal=0)
    try:
        if not dataset.is_calibration_run():
            logging.warning("%s is not a cal-pulser run (enable_cal=0 in run config); "
                            "selected events may be noise/RF near the PPS", rundir)
    except ValueError as e:
        logging.warning("Could not determine cal-pulser status for %s: %s", rundir, e)

    dataset.setEntries((0, dataset.N()))
    infos = dataset.eventInfo()

    event_number = np.array([i.eventNumber for i in infos])
    trigger_type = np.array([i.triggerType for i in infos])
    sysclk = np.array([i.sysclk for i in infos], dtype=np.int64)
    pps0 = np.array([i.sysclkLastPPS[0] for i in infos], dtype=np.int64)
    pps1 = np.array([i.sysclkLastPPS[1] for i in infos], dtype=np.int64)

    # ticks per second from the two most recent PPS sysclks (uint32, handle wrap)
    f = (pps0 - pps1) % (2 ** 32)
    # time since the last PPS (handle uint32 wrap-around)
    dt = ((sysclk - pps0) % (2 ** 32)) / f

    cal_evs = set(event_number[dt < threshold].tolist())
    force_evs = set(event_number[trigger_type == "FORCE"].tolist())
    cal_label = run_cal_label(dataset)
    return (cal_evs, force_evs,
            int(dataset.station), int(dataset.run), dataset.duration(), cal_label)


def saturation_limits_for_run(rundir, override, margin):
    """
    Per-channel max_abs_amplitude limit (ADC counts) above which a channel saturates.

    With ``override`` the same value is used for every channel. Otherwise each channel's
    limit is its own pedestal near-rail distance (from pedestal.root) reduced by
    ``margin`` to allow for per-sample pedestal spread. Returns (limits, source) where
    ``limits`` is a dict {channel: adc} over the deep channels, or (None, reason) if
    unavailable.
    """
    if override is not None:
        return {ch: float(override) for ch in DEEP_CHANNELS}, "user-specified"

    near_rail = pedestal_near_rail(rundir)
    if near_rail is None:
        return None, "pedestal.root not found"

    limits = {ch: (1.0 - margin) * float(near_rail[ch]) for ch in DEEP_CHANNELS}
    return limits, f"pedestal-derived (per-channel near-rail x {1 - margin:.2f})"


def combine_limits(per_run_limits):
    """Combine per-run, per-channel saturation limits into one conservative (min) limit per channel."""
    valid = [lim for lim in per_run_limits if lim is not None]
    if not valid:
        return None
    return {ch: min(lim[ch] for lim in valid) for ch in DEEP_CHANNELS}


def collect_calpulser_amplitudes(monitoring_files, rundirs, threshold,
                                 saturation_override=None, saturation_margin=0.02):
    """Gather deep-channel max-amplitude arrays for all cal-pulser events across runs.

    Also prints, per run, the cal-pulser trigger efficiency (detected / run-duration
    seconds, since the pulser fires once per second) and warns about ADC saturation.
    """
    # per channel -> list of amplitudes (cal pulsers) / vrms (force triggers)
    amps = defaultdict(list)
    force_rms = defaultdict(list)
    n_cal_total = 0
    runs = []  # (station, run) of every processed run
    per_run_limits = []  # per-channel saturation limits of every run

    cal_labels = {}  # station -> set of cal-pulser channel/type labels seen
    for mon_path, rundir in zip(monitoring_files, rundirs):
        amplitudes, rms = read_monitoring_summary(mon_path)
        cal_evs, force_evs, station, run, duration, cal_label = calpulser_event_numbers(rundir, threshold)
        runs.append((station, run))
        cal_labels.setdefault(station, set()).add(cal_label)

        # force-trigger vrms (noise level), matched to the monitoring file
        force_matched = force_evs & rms.keys()
        for ev in force_matched:
            r = rms[ev]
            for ch in DEEP_CHANNELS:
                force_rms[ch].append(r[ch])

        matched = cal_evs & amplitudes.keys()
        missing = cal_evs - amplitudes.keys()
        n_cal_total += len(cal_evs)
        print(f"station{station}/run{run}: {len(cal_evs)} cal-pulser events, "
              f"{len(matched)} matched in monitoring file"
              + (f" ({len(missing)} not in monitoring file)" if missing else "")
              + f"; {len(force_matched)} force triggers for vrms")

        # trigger efficiency: one expected cal pulse per second of run duration
        if duration and duration > 0:
            print(f"    trigger efficiency: {len(cal_evs)} / {duration:.0f} s "
                  f"= {len(cal_evs) / duration:.3f}")
        else:
            logging.warning("    could not determine run duration; skipping efficiency")

        # collect per-channel amplitudes for this run (for histogram + saturation check)
        run_amps = {ch: np.array([amplitudes[ev][ch] for ev in matched]) for ch in DEEP_CHANNELS}
        for ch in DEEP_CHANNELS:
            amps[ch].extend(run_amps[ch].tolist())

        limits, source = saturation_limits_for_run(rundir, saturation_override, saturation_margin)
        per_run_limits.append(limits)
        warn_saturation(run_amps, limits, source)

    print(f"Total cal-pulser events selected: {n_cal_total}")

    # mean force-trigger vrms per channel (noise floor), pooled over all runs
    vrms_levels = {ch: float(np.mean(force_rms[ch])) for ch in DEEP_CHANNELS if force_rms[ch]}

    return ({ch: np.array(v) for ch, v in amps.items()}, runs,
            combine_limits(per_run_limits), vrms_levels or None,
            format_station_title(cal_labels))


def format_station_title(cal_labels):
    """Build a per-station title like "station23 (coax/pulser)" from a
    {station: set(cal-pulser labels)} mapping, joining multiple stations with ", "."""
    parts = []
    for station in sorted(cal_labels):
        labels = sorted(lbl for lbl in cal_labels[station] if lbl is not None)
        suffix = f" ({', '.join(labels)})" if labels else ""
        parts.append(f"station{station}{suffix}")
    return ", ".join(parts)


def warn_saturation(run_amps, limits, source):
    """Warn if cal-pulser events reach the per-channel ADC saturation limit on deep channels."""
    if limits is None:
        logging.warning("    saturation check skipped (%s)", source)
        return

    flagged = {}
    for ch in DEEP_CHANNELS:
        data = run_amps[ch]
        if data.size:
            n = int(np.sum(data >= limits[ch]))
            if n:
                flagged[ch] = (n, data.size, limits[ch])

    if flagged:
        details = ", ".join(f"ch{ch}: {n}/{tot} (>= {lim:.0f} ADC, {100 * n / tot:.1f}%)"
                            for ch, (n, tot, lim) in flagged.items())
        logging.warning("    SATURATION (%s): %s", source, details)
    else:
        print(f"    no saturation (per-channel limits, {source})")


def plot_histograms(amps, output, bins=50, sat_limits=None, vrms_levels=None, title=None):
    """Grid of max-abs-amplitude histograms, one panel per deep channel.

    If ``sat_limits`` (dict {channel: adc}) is given, each channel's saturation limit
    is drawn as a red vertical line. If ``vrms_levels`` (dict {channel: adc}) is given,
    each channel's mean force-trigger vrms (noise floor) is drawn as a green line.
    If ``title`` is given, it is used as the figure suptitle (station + cal-pulser type).
    """
    n = len(DEEP_CHANNELS)
    ncols = 5
    nrows = int(np.ceil(n / ncols))

    fig, axs = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2 * nrows),
                            sharex=True, sharey=True, squeeze=False)

    for idx, ch in enumerate(DEEP_CHANNELS):
        ax = axs[idx // ncols][idx % ncols]
        data = amps.get(ch, np.array([]))
        if data.size:
            ax.hist(data, bins=bins)
        if vrms_levels is not None and ch in vrms_levels:
            ax.axvline(vrms_levels[ch], color="g", ls=":", lw=1,
                       label=f"$V_{{rms}}$ = {vrms_levels[ch]:.1f} ADC")
        if sat_limits is not None and ch in sat_limits:
            ax.axvline(sat_limits[ch], color="r", ls="--", lw=1,
                       label=f"Sat. = {sat_limits[ch]:.0f} ADC")
        ax.legend(title=f"Channel {ch} (N={data.size})",
                  fontsize="x-small", title_fontsize="small", ncols=1)
        ax.grid()

    # hide unused panels
    for idx in range(n, nrows * ncols):
        axs[idx // ncols][idx % ncols].set_visible(False)


    fig.supxlabel("Max. abs. amplitude")
    fig.supylabel("Cal-pulser events", x=0)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output)
    print(f"Saved histogram to {output}")


def resolve_paths(path, monitoring_name="monitoring.root"):
    """
    Resolve a user-supplied path into (monitoring_file, run_directory).

    ``path`` may be a run directory or either of the two files it contains
    (monitoring.root / headers.root); both live in the same directory.
    """
    rundir = path if os.path.isdir(path) else os.path.dirname(path) or "."
    monitoring_file = os.path.join(rundir, monitoring_name)
    return monitoring_file, rundir


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+",
                        help="run directory or monitoring.root/headers.root file (one per run)")
    parser.add_argument("--monitoring-name", default="monitoring.root",
                        help="name of the monitoring file inside the run directory")
    parser.add_argument("--threshold", type=float, default=3e-6,
                        help="max (t_trigger - t_pps) in seconds to tag a cal pulser (default: 3e-6)")
    parser.add_argument("--bins", type=int, default=50, help="number of histogram bins")
    parser.add_argument("--saturation-threshold", type=float, default=None,
                        help="fixed max_abs_amplitude (ADC counts) flagging saturation; "
                             "default is derived from pedestal.root (smallest deep near-rail)")
    parser.add_argument("--saturation-margin", type=float, default=0.02,
                        help="fractional margin below the near-rail for the pedestal-derived "
                             "default saturation threshold (default: 0.02)")
    parser.add_argument("--output", type=str, default=None,
                        help="output figure path (default: calpulser_max_amplitude_<station><runs>.png)")
    args = parser.parse_args()

    monitoring_files, rundirs = zip(*(resolve_paths(p, args.monitoring_name) for p in args.paths))

    amps, runs, sat_limits, vrms_levels, station_title = collect_calpulser_amplitudes(
        monitoring_files, rundirs, args.threshold,
        saturation_override=args.saturation_threshold,
        saturation_margin=args.saturation_margin)

    output = args.output
    if output is None:
        tag = "_".join(f"s{station}r{run}" for station, run in runs)
        output = f"calpulser_max_amplitude_{tag}.png"

    plot_histograms(amps, output, bins=args.bins,
                    sat_limits=sat_limits, vrms_levels=vrms_levels, title=station_title)


if __name__ == "__main__":
    main()
