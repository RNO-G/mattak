import numpy as np
from matplotlib import pyplot as plt
import mattak.Dataset
import argparse
import datetime as dt
from collections import defaultdict
import matplotlib.dates as mdates


'''
This scripts expects proper (can be fake) waveforms.
'''


def calculate_power(traces):
    return np.sum(traces ** 2, axis=-1)


def plot_comparison_violin(axs, wfs1, wfs2):

    amax = np.amax(wfs1, axis=-1) - np.amax(wfs2, axis=-1)
    amin = np.amin(wfs1, axis=-1) - np.amin(wfs2, axis=-1)
    power = calculate_power(wfs1)
    power2 = calculate_power(wfs2)
    power_ratio = (power - power2) / power

    rms = np.sqrt(np.mean(wfs1 ** 2, axis=-1)) - np.sqrt(np.mean(wfs2 ** 2, axis=-1))

    label = ""
    if wfs2.ndim == 1:
        label = "2.5V / 4095"

    _plot_plot_comparison_violin(axs, amax, amin, rms, power_ratio, label)

def _plot_plot_comparison_violin(axs, amax, amin, rms, power_ratio, label=""):


    parts2 = axs[0].violinplot(
        amax, np.arange(24), showextrema=True, showmedians=True,
        vert=False, side="high", widths=1.8)

    axs[0].axvline(0, color="k", ls="--", lw=1, label=label)
    axs[0].set_xlabel("max amplitude / mV")

    parts2 = axs[1].violinplot(
        amin, np.arange(24), showextrema=True, showmedians=True,
        vert=False, side="high", widths=1.8)

    axs[1].axvline(0, color="k", ls="--", lw=1, label=label)
    axs[1].set_xlabel("min amplitude / mV")

    parts2 = axs[2].violinplot(
        rms, np.arange(24), showextrema=True, showmedians=True,
        vert=False, side="high", points=300, widths=1.8)

    axs[2].axvline(0, color="k", ls="--", lw=1, label=label)
    axs[2].set_xlabel(r"$\Delta$RMS / mV")

    parts2 = axs[3].violinplot(
        power_ratio, np.arange(24), showextrema=True, showmedians=True,
        vert=False, side="high", points=300, widths=1.8)

    axs[3].axvline(0, color="k", ls="--", lw=1, label=label)
    axs[3].set_xlabel(r"power ratio")

    if label != "":
        axs[0].legend()

    axs[0].set_ylabel("channel id")


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="")

    argparser.add_argument('data_path', type=str, default=None)
    argparser.add_argument('--voltage_calibration', '-vc', nargs="+", type=str, default=None)
    argparser.add_argument('--label', '-l', type=str, default="")
    argparser.add_argument('--backend', '-b', type=str, default="pyroot")
    argparser.add_argument('--lin', action="store_true")

    args = argparser.parse_args()

    dset1 = mattak.Dataset.Dataset(
        0, 0, args.data_path,
        backend=args.backend, verbose=True, read_daq_status=False,
        voltage_calibration=args.voltage_calibration[0])

    if args.label != "" and not args.label.startswith("_"):
        args.label = f"_{args.label}"

    # for the uncalibrate data every waveform is the same (besides position)
    dset1.setEntries(1)
    orig_adc_trace = np.array(dset1.wfs(), dtype=float)
    adc_ref = int(np.amax(orig_adc_trace))
    orig_adc_trace = orig_adc_trace * 2500 / 4095

    if not args.lin:
        dset1.setEntries((0, dset1.N()))
        wfs1 = dset1.wfs(calibrated=True)
        wfs1 = wfs1 * 1000 # convert to mV

    st = dset1.vc.getStationNumber()
    t1 = dt.datetime.fromtimestamp(dset1.vc.getStartTime())
    t1_str = t1.strftime('%Y.%m.%d')

    if len(args.voltage_calibration) == 1:
        # compare with linear "pseudo" calibration
        fig, axs = plt.subplots(1, 4, sharey=True, figsize=(10, 8))
        plot_comparison_violin(axs, wfs1, orig_adc_trace)

        fig.suptitle(r"$A_\mathrm{max}^\mathrm{ref}$"
                     rf" = {adc_ref} ADC $\approx$ {adc_ref * 2500 / 4096:.2f}mV")
        fig.tight_layout()
        fig.savefig(
            f"voltage_calib_waveform_st{st}_{t1_str}_ampl{adc_ref}{args.label}.png")

    elif len(args.voltage_calibration) == 2:
        # compare to calibrations
        dset2 = mattak.Dataset.Dataset(
            0, 0, args.data_path,
            backend=args.backend, verbose=True, read_daq_status=False,
            voltage_calibration=args.voltage_calibration[1])

        dset2.setEntries((0, dset2.N()))
        wfs2 = dset2.wfs(calibrated=True)
        wfs2 = wfs2 * 1000 # convert to mV

        fig, axs = plt.subplots(1, 4, sharey=True, figsize=(10, 8))
        plot_comparison_violin(axs, wfs1, wfs2)

        t2_str = dt.datetime.fromtimestamp(dset2.vc.getStartTime()).strftime('%Y.%m.%d')

        fig.suptitle(
            r"$A_\mathrm{max}^\mathrm{ref}$"
            rf" = {adc_ref} ADC $\approx$ {adc_ref * 2500 / 4096:.2f}mV; {t1_str} vs {t2_str}")

        fig.tight_layout()
        tstr = t1_str + "-" + t2_str
        fig.savefig(
            f"diff_voltage_calib_waveform_st{st}_{tstr}_ampl{adc_ref}{args.label}.png")

    elif args.lin:
        data = defaultdict(list)
        for vc in args.voltage_calibration:

            dset = mattak.Dataset.Dataset(
                0, 0, args.data_path,
                backend="pyroot", verbose=True, read_daq_status=False,
                voltage_calibration=vc)

            dset.setEntries((0, dset.N()))
            wfs = dset.wfs(calibrated=True)
            wfs = wfs * 1000 # convert to mV

            data["amax"].append(np.amax(wfs, axis=-1))
            data["amin"].append(np.amin(wfs, axis=-1))
            data["rms"].append(np.sqrt(np.mean(wfs ** 2, axis=-1)))
            data["power"].append(calculate_power(wfs))

            data["times"].append(
                dt.datetime.fromtimestamp(dset.vc.getStartTime())
            )


        data = {key: np.array(value) for key, value in data.items()}
        times = data["times"]

        fig, axs = plt.subplots(1, 4, sharey=True, figsize=(10, 8))

        amax = np.vstack(data["amax"]) - np.amax(orig_adc_trace)
        print(amax.shape)
        amin = np.vstack(data["amin"]) - np.amin(orig_adc_trace)
        rms = np.vstack(data["rms"]) - np.sqrt(np.mean(orig_adc_trace ** 2))
        power = np.vstack(data["power"])
        power_ratio = (power - calculate_power(orig_adc_trace)) / power

        label = "2.5V / 4095"
        _plot_plot_comparison_violin(axs, amax, amin, rms, power_ratio, label=label)

        t2_str = times[-1].strftime('%Y.%m.%d')
        fig.suptitle(
            r"$A_\mathrm{max}^\mathrm{ref}$"
            rf" = {adc_ref} ADC $\approx$ {adc_ref * 2500 / 4096:.2f}mV; {t1_str} - {t2_str}")

        tstr = t1_str + "-" + t2_str + f"-{len(times)}"

        fig.tight_layout()
        fig.savefig(
            f"diff_from_lin_voltage_calib_waveform_st{st}_{tstr}_ampl{adc_ref}{args.label}.png")

    else:

        data = defaultdict(list)
        for vc in args.voltage_calibration[1:]:

            dset = mattak.Dataset.Dataset(
                0, 0, args.data_path,
                backend="pyroot", verbose=True, read_daq_status=False,
                voltage_calibration=vc)

            dset.setEntries((0, dset.N()))
            wfs = dset.wfs(calibrated=True)
            wfs = wfs * 1000 # convert to mV

            amax = np.amax(wfs1, axis=-1) - np.amax(wfs, axis=-1)
            amin = np.amin(wfs1, axis=-1) - np.amin(wfs, axis=-1)

            power1 = calculate_power(wfs1)
            power = calculate_power(wfs)
            power_ratio = (power1 - power) / power1

            data["power_ratios_mean"].append(np.mean(power_ratio, axis=0))
            data["amaxs_mean"].append(np.mean(amax, axis=0))
            data["amins_mean"].append(np.mean(amin, axis=0))

            data["power_ratios_std"].append(np.std(power_ratio, axis=0))
            data["amaxs_std"].append(np.std(amax, axis=0))
            data["amins_std"].append(np.std(amin, axis=0))

            data["times"].append(
                dt.datetime.fromtimestamp(dset.vc.getStartTime())
            )

        data = {key: np.array(value) for key, value in data.items()}
        times = data["times"]

        years = np.unique([t.year for t in times])

        def _plot(ax_):
            ax2_ = ax_.twinx()
            for y, yerr in zip(data["amaxs_mean"].T, data["amaxs_std"].T):
                ax_.errorbar(times, y, yerr, color='C0', alpha=0.3, ls="", marker="s")

            for y, yerr in zip(data["amins_mean"].T, data["amins_std"].T):
                ax_.errorbar(times, y, yerr, color='C1', alpha=0.3, ls="", marker="s")

            for y, yerr in zip(data["power_ratios_mean"].T, data["power_ratios_std"].T):
                ax2_.errorbar(times, y, yerr, color='C2', alpha=0.3, ls="", marker="s")

            ax_.errorbar(np.nan, np.nan, 0,
                        color='C0', alpha=0.3, ls="", marker="s", label=r"$A_\mathrm{max}$")

            ax_.errorbar(np.nan, np.nan, 0,
                        color='C1', alpha=0.3, ls="", marker="s", label=r"$A_\mathrm{min}$")

            ax_.errorbar(np.nan, np.nan, 0,
                        color='C2', alpha=0.3, ls="", marker="s", label=r"power ratio")

            return ax2_

        if len(years) == 1:

            fig, ax = plt.subplots()

            ax2 = _plot(ax)
            ax.axvline(t1, color="k", ls="--", label="reference calibration")
            ax.legend()

            ax.set_ylabel(r"$\Delta A$ / mV")
            ax2.set_ylabel(r"power ratio")
            ax.set_title(
                r"$A_\mathrm{max}^\mathrm{ref}$"
                rf" = {adc_ref} ADC $\approx$ {adc_ref * 2500 / 4096:.2f}mV")

            ax.tick_params(axis="x", rotation=25)

        else:
            fig, axs = plt.subplots(1, 2, sharey=True, width_ratios=[1, 3])
            fig.subplots_adjust(wspace=0.05, left=0.1, right=0.85, top=0.9, bottom=0.13)  # adjust space between Axes

            ax21 = _plot(axs[0])
            ax22 = _plot(axs[1])

            ax21.set_axis_off()
            y11, y12 = ax21.get_ylim()
            y21, y22 = ax22.get_ylim()
            ax21.set_ylim(np.amin([y11, y12]), np.amax([y21, y22]))
            ax21.set_ylim(np.amin([y11, y12]), np.amax([y21, y22]))

            axs[0].set_xlim(dt.datetime(2022, 9, 26), dt.datetime(2022, 10, 5))
            axs[1].set_xlim(dt.datetime(2023, 4, 23), None)

            # hide the spines between ax and ax2
            axs[0].spines.right.set_visible(False)
            axs[1].spines.left.set_visible(False)
            ax22.spines.left.set_visible(False)
            axs[1].tick_params(left=False)  # don't put tick on left axis

            d = .5  # proportion of vertical to horizontal extent of the slanted line
            kwargs = dict(marker=[(-1, -d), (1, d)], markersize=8,
                        linestyle="none", color='k', mec='k', mew=1, clip_on=False)
            axs[0].plot([1, 1], [0, 1], transform=axs[0].transAxes, **kwargs)
            axs[1].plot([0, 0], [0, 1], transform=axs[1].transAxes, **kwargs)

            axs[0].set_xticks([dt.datetime(2022, 9, 27), dt.datetime(2022, 10, 1)])

            for ax in axs:
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))
                ax.tick_params(axis="x", rotation=25)

            axs[0].set_ylabel(r"$\Delta A$ / mV")
            ax22.set_ylabel(r"power ratio")

            axs[1].axvline(t1, color="k", ls="--", label="reference calibration")
            axs[1].legend()

            fig.suptitle(
                r"$A_\mathrm{max}^\mathrm{ref}$"
                rf" = {adc_ref} ADC $\approx$ {adc_ref * 2500 / 4096:.2f}mV")

        # fig.tight_layout()
        tstr = t1_str + "-" + times[-1].strftime('%Y.%m.%d') + f"-{len(times)}"
        fig.savefig(
            f"diff_voltage_calib_waveform_st{st}_{tstr}_ampl{adc_ref}{args.label}.png")
