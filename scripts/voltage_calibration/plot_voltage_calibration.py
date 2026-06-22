from mattak.backends.uproot.voltage_calibration import VoltageCalibration
import sys
import numpy as np
import os
from matplotlib import pyplot as plt
import mattak.Dataset
import argparse

"""
This script needs to read in pseudo data produced with scripts/makeFakeDataToTestVoltageCalibration.C
"""


def get_data(dset, vmin=-1.3, vmax=0.7):
    dset.setEntries((0, dset.N()))
    wfs = dset.wfs(calibrated=True) * 1000 # convert to mV

    wfs[wfs > vmax * 1000] = vmax * 1000
    wfs[wfs < vmin * 1000] = vmin * 1000

    wfs = np.swapaxes(wfs, 0, 1)
    wfs = wfs.reshape(24, wfs.shape[1] // 2, 4096)
    wfs = np.swapaxes(wfs, 0, 1)

    out = np.swapaxes([np.mean(wfs, axis=-1), np.std(wfs, axis=-1)], 0, 1)

    return out, wfs


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="")

    argparser.add_argument('data_path', type=str, default=None)
    argparser.add_argument('--voltage_calibration', '-vc', type=str, default=None)
    argparser.add_argument('--bias_scan', type=str, default=None)
    argparser.add_argument('--calibrate', action="store_true")
    args = argparser.parse_args()

    dset_pyroot = mattak.Dataset.Dataset(
        0, 0, args.data_path,
        backend="pyroot", verbose=True, read_daq_status=False,
        voltage_calibration=args.voltage_calibration)

    vc = VoltageCalibration(args.voltage_calibration, caching_mode="interpolate", table_step_size=0)

    dset_uproot = mattak.Dataset.Dataset(
        0, 0, args.data_path,
        backend="uproot", verbose=True, read_daq_status=False,
        voltage_calibration=vc)

    vc2 = VoltageCalibration(args.bias_scan, caching_mode="interpolate", table_step_size=0)

    dset_uproot2 = mattak.Dataset.Dataset(
        0, 0, args.data_path,
        backend="uproot", verbose=True, read_daq_status=False,
        voltage_calibration=vc2)

    label = os.path.basename(args.voltage_calibration).replace("volCalConsts_pol9_", "").replace(".root", "")

    xs = np.arange(-1000, 1000, 50)

    linear_calib = 2500 / 4095  # mV / ADC

    out_uproot, wfs_uproot = get_data(dset_uproot)
    out_uproot2, wfs_uproot2 = get_data(dset_uproot2)
    out_pyroot, wfs_pyroot = get_data(dset_pyroot, vc.fit_min, vc.fit_max)
    print(wfs_pyroot.shape, out_pyroot.shape)

    fig, (ax, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10, 6))
    ax.plot(xs, xs * linear_calib, "k--", lw=1, label=f"{linear_calib:.3f} mV / ADC")

    for ch in range(24):

        ax.errorbar(xs, out_pyroot[:, 0, ch], out_pyroot[:, 1, ch], elinewidth=2, alpha=0.5, ls="", color="C0")
        ax2.errorbar(xs, out_pyroot[:, 0, ch] - xs * linear_calib, out_pyroot[:, 1, ch], elinewidth=2, alpha=0.5, ls="", color="C0")

        # ax.errorbar(xs, out_uproot[:, 0, ch], out_uproot[:, 1, ch], elinewidth=2, alpha=0.5, ls="", label="uproot - calibration file")
        # ax2.errorbar(xs, out_uproot[:, 0, ch] - xs * linear_calib, out_uproot[:, 1, ch], elinewidth=2, alpha=0.5, ls="")

        # ax.errorbar(xs, out_uproot2[:, 0, ch], out_uproot2[:, 1, ch], elinewidth=2, alpha=0.5, ls="", label="uproot - bias scan file")
        # ax2.errorbar(xs, out_uproot2[:, 0, ch] - xs * linear_calib, out_uproot2[:, 1, ch], elinewidth=2, alpha=0.5, ls="")

        ax3.plot(xs, out_uproot[:, 0, ch] - out_pyroot[:, 0, ch], marker="s", markersize=2, alpha=0.5, ls="", color="C1")
        ax3.plot(xs, out_uproot2[:, 0, ch] - out_pyroot[:, 0, ch], marker="s", markersize=2, alpha=0.5, ls="", color="C2")

    ax3.plot(xs, out_uproot[:, 0, ch] - out_pyroot[:, 0, ch], marker="s", markersize=2, alpha=0.5, ls="", color="C1", label="uproot - calibration file")
    ax3.plot(xs, out_uproot2[:, 0, ch] - out_pyroot[:, 0, ch], marker="s", markersize=2, alpha=0.5, ls="", color="C2", label="uproot - bias scan file")
    ax.errorbar(np.nan, np.nan, 0, elinewidth=2, alpha=0.5, ls="", color="C0", label="pyroot")

    ax3.set_ylabel("uproot - pyroot / mV")
    ax3.set_xlabel("input ADC")
    ax3.grid()

    ax2.set_ylabel("residual from\nlinear / mV")
    ax2.grid()

    ax.set_ylabel("output voltage / mV")
    ax.grid()
    ax.legend(ncols=2)
    ax3.legend()

    fig.tight_layout()
    plt.savefig(f"compare_calibrations_all_channels_{label}.png", transparent=False, dpi=300)


    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 6))

    for ch in range(24):
        ax.plot(xs * linear_calib, out_pyroot[:, 0, ch] - xs * linear_calib, marker="s", markersize=2, alpha=0.3, ls="", color="C0")


    diff = wfs_uproot - wfs_pyroot
    diff = diff.reshape(len(diff), -1)
    ax2.errorbar(xs * linear_calib, np.mean(diff, axis=-1), np.std(diff, axis=-1), marker="s", markersize=2, alpha=0.3, ls="", color="C1", label="uproot - calibration file")

    diff = wfs_uproot2 - wfs_pyroot
    diff = diff.reshape(len(diff), -1)
    ax2.errorbar(xs * linear_calib, np.mean(diff, axis=-1), np.std(diff, axis=-1), marker="s", markersize=2, alpha=0.3, ls="", color="C2", label="uproot - interpolate")


    ax2.set_ylabel("uproot - pyroot / mV")
    ax2.set_xlabel("ADC in mV")
    ax2.grid()

    ax.set_ylabel("residual from\nlinear / mV")
    ax.grid()

    fig.tight_layout()
    plt.savefig(f"compare_calibrations_channels_samples_{label}.png", transparent=False, dpi=300)
