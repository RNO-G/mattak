# `scripts/`

Helper scripts and tools that live alongside the mattak library. They are
grouped by purpose into the sub-directories described below.

A subset of these is installed onto the `PATH` (into `${prefix}/bin`) by the
`INSTALL_SCRIPTS` list in [`../CMakeLists.txt`](../CMakeLists.txt); those are
marked **(installed)** below. The rest are run in place from the source tree
(e.g. via the systemd units in [`services/`](services), or by hand).

## `daq/`

Scripts driving the data-acquisition conversion pipeline (raw station data →
ROOT files). These are the production scripts used by the DAQ.

- `rno-g-autoconverter` — long-running converter loop, run per station by the
  systemd service in [`services/`](services). Watches the ingress directory and
  rootifies new runs, builds combined files and event lists.
- `rno-g-convert-run` **(installed)** — convert a single raw run directory into
  the per-run ROOT files (`waveforms.root`, `headers.root`, `daqstatus.root`,
  `pedestal.root`, `runinfo.root`) and trigger calibration/monitoring steps.
- `rno-g-convert-station` **(installed)** — convert all runs of a station in
  parallel (wraps `rno-g-convert-run` via GNU `parallel`).
- `rno-g-reconvert` — re-run the conversion for a single station/run.
- `rno-g-make-eventlists` **(installed)** — produce per-run event lists so files
  can be remade without reconverting everything.
- `rno-g-check-daqstatus` — find runs missing `daqstatus.root` and convert just
  the daqstatus for them.

## `metadata/`

Tools for collecting and inspecting run/event metadata.

- `collect_runinfo.py` — gather per-run information into a summary.
- `collect_event_info_uproot.py` — collect per-event information via `uproot`
  (uses `collect_runinfo`).
- `read_daqstatus.py` — read/plot fields from `daqstatus.root`.
- `plot_flower_meta_data.py` — plot FLOWER (trigger board) metadata.

## `monitoring/`

Run-level monitoring data and diagnostic plots.

- `mon_write.py` **(installed)** — write `monitoring.root` for a run (invoked by
  `rno-g-convert-run`).
- `mon_read.py`, `mon_read2.py` — example readers for `monitoring.root`.
- `plot_calpulser_amplitudes.py` — select cal-pulser events and histogram their
  maximum amplitude.
- `plot_calpulser_dt.py` — plot the dt = (t_trigger − t_pps) distribution to help
  tune the cal-pulser time threshold.

## `voltage_calibration/`

Bias-scan voltage calibration: tooling, plots, fake-data generators, and
templates.

- `convert_bias_scan_to_calibration.py` **(installed)** — convert a bias scan
  into calibration constants (called by `rno-g-convert-run`).
- `rno-g-run-missing-voltage-calibration` — find runs without calibration and
  generate it.
- `rno-g-rename-calibration-files` — rename `volCalConsts_*` files from timestamp
  to run-number naming.
- `plot_voltage_calibration.py`, `plot_voltage_calibration_waveforms.py` —
  diagnostic plots.
- `makeFakeData.C`, `makeFakeDataToTestVoltageCalibration*.cc` — generators for
  test data used to validate the voltage calibration.
- `template4_wo_hardwareResponse_max_amp_*.txt` — input templates.

## `services/`

systemd units and helpers for running the DAQ converter as a service.

- `rno-g-autoconverter@.service`, `rno-g-autoconverter.target` — per-station
  service template and target. Installed with `make install-systemd` (see the
  top-level [`Makefile`](../Makefile)).
- `tmux_follow_services.sh` — open a tmux session tailing each station's journal.
- `README.md` — service start/stop/log commands.

## `plots/`

Output directory for generated plots (PNGs) from the plotting scripts above.
