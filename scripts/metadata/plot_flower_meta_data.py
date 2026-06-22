import numpy as np
import matplotlib.pyplot as plt
import sys
import datetime
import matplotlib.dates as mdates
from utilities import plotting
import collections

def convert_thresholds_to_voltages(thresholds, gains, nevents):
    gains = np.array([np.repeat(ele, nevents).astype(int) for ele in gains.T]).T
    gains[gains > 13] = 13

    amplification_factors = np.array([1, 1.25, 2, 2.5, 4, 5, 8, 10, 12.5, 16, 20, 25, 32, 50])

    amplification_factors = np.array([amplification_factors[gain] for gain in gains])

    thresholds_volt = thresholds * 2 / 255 / amplification_factors

    return thresholds_volt


def get_frist_event_time_per_run(trigger_time, number_of_events):
    run_time = []
    idx = 0
    for n_events in number_of_events:
        jdx = 0
        while np.isinf(triggerTime[idx + jdx]):
            jdx += 1

        time = triggerTime[idx + jdx]

        run_time.append(triggerTime[idx + jdx])
        idx += n_events

    return np.array(run_time)


def plots_1(data):
    run_numbers = data["run_number"]
    number_of_events = data["number_of_events"]
    lowTrigThrs = data["lowTrigThrs"]
    triggerTime = data["triggerTime"]
    triggerType = data["triggerType"]
    flower_gain_codes = data["flower_gain_codes"]

    time_gain_dt = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in time_gain])

    mask = np.all([triggerType == "LT", triggerTime > datetime.datetime(2023, 1, 1).timestamp(), ~np.isinf(triggerTime)], axis=0)
    triggerTime = triggerTime[mask]
    lowTrigThrs = lowTrigThrs[mask]
    lowTrigThrs_volt = lowTrigThrs_volt[mask]

    triggerTime_dt = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in triggerTime])

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    axs[0].plot(triggerTime_dt[::100], lowTrigThrs[::100], "^")

    mask2 = time_gain_dt >  datetime.datetime(2023, 1, 1)
    axs[1].plot(time_gain_dt[mask2], flower_gain_codes[mask2], "v")

    axs[0].set_ylabel("trigger thresholds / ADC")

    axs[1].set_ylabel("flower gain code")

    axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))
    # ax.xaxis.set_major_locator(mdates.DayLocator())

    fig.tight_layout()
    fig.savefig("flower_threshold_adc_and_gain.png")

    fig2, axs2 = plt.subplots(1, 2)
    for ele in lowTrigThrs.T:
        over_under = np.sum(np.any([ele < 15, ele > 40], axis=0))
        plotting.pretty_trans_hist(axs2[0], ele, bins=np.arange(15, 41), label=over_under)
        # axs2[0].hist(ele, bins=np.arange(15, 41))

    lowTrigThrs_mean = np.mean(lowTrigThrs, axis=1)
    print(np.mean(lowTrigThrs_mean), np.std(lowTrigThrs_mean))
    over_under = np.sum(np.any([lowTrigThrs_mean < 15, lowTrigThrs_mean > 40], axis=0))
    plotting.pretty_trans_hist(axs2[0], lowTrigThrs_mean, bins=np.arange(15, 41), color="k", label=over_under)

    # axs2[0].hist(np.nan, bins=np.arange(15, 41), color="k", label=f"over/under threshold: {over_under}")
    axs2[0].set_xlabel("trigger thresholds / ADC")

    for ele in flower_gain_codes.T:
        plotting.pretty_trans_hist(axs2[1], ele, bins=np.arange(flower_gain_codes.min(), flower_gain_codes.max() + 1))
        # axs2[1].hist(ele[mask2], bins=np.arange(flower_gain_codes.min(), flower_gain_codes.max() + 1))
    axs2[1].set_xlabel("flower gain code")

    axs2[0].legend(title="over/under threshold", ncol=2)
    fig2.tight_layout()
    fig2.savefig("flower_threshold_adc_and_gain_hist.png")


    fig3, ax3 = plt.subplots()
    ax3.plot(triggerTime_dt[::100], lowTrigThrs_volt[::100] * 1000, "^")
    ax3.set_ylabel("trigger thresholds / mV")

    fig3.savefig("flower_threshold_volt.png")
    ax3.set_ylim(5, 10)
    fig3.savefig("flower_threshold_volt_zoom.png")

    fig4, ax4 = plt.subplots()
    for ele in lowTrigThrs_volt.T:
        ele *= 1000
        over_under = np.sum(np.any([ele < 5, ele > 10], axis=0))
        plotting.pretty_trans_hist(ax4, ele, bins=np.linspace(5, 10), label=over_under)

    ax4.legend(title="over/under threshold", ncol=2)
    ax4.set_xlabel("trigger thresholds / mV")
    fig4.savefig("flower_threshold_volt_hist.png")


def apply_mask(data, mask, run_level=False):

    if run_level:
        mask2 = np.repeat(mask, data["number_of_events"])

    for key in data:
        if run_level and len(data[key]) != len(mask):
            data[key] = data[key][mask2]
        else:
            data[key] = data[key][mask]

    return data


def per_run_data(arr, number_of_events):
    idx = 0
    data = []
    for n_events in number_of_events:
        data.append(arr[idx:idx + n_events])
        idx += n_events

    return data


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_flower_meta_data.py <flower_meta_data.npz>")
        sys.exit(1)

    data = collections.defaultdict(list)
    for path in sys.argv[1:]:
        tmp_data = np.load(path)
        for key in tmp_data:
            data[key].extend(tmp_data[key].tolist())

    data = {key: np.array(data[key]) for key in data}

    run_numbers = data["run_number"]
    number_of_events = data["number_of_events"]
    lowTrigThrs = data["lowTrigThrs"]
    triggerTime = data["triggerTime"]
    flower_gain_codes = data["flower_gain_codes"]

    if isinstance(flower_gain_codes[0][0], str):
        flower_gain_codes = np.array([[int(x) for x in code.split(" ")] for code in flower_gain_codes])


    lowTrigThrs_volt = convert_thresholds_to_voltages(lowTrigThrs, flower_gain_codes, number_of_events)
    run_time = get_frist_event_time_per_run(triggerTime, number_of_events)

    data["run_time"] = run_time
    data["lowTrigThrs_volt"] = lowTrigThrs_volt

    lowTrigThrs_run = per_run_data(lowTrigThrs, number_of_events)
    lowTrigThrs_run_mean = np.array([np.mean(ele) for ele in lowTrigThrs_run])

    fig, ax = plt.subplots()

    run_time_dt = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in run_time])

    mask = data["duration"] > 7000

    lt_trigger_rate = data["lt_trigger_rate"]
    ax.plot(run_time_dt[mask], lt_trigger_rate[mask], "o")
    # ax.plot(run_time_dt[~mask], lt_trigger_rate[~mask], "o")
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))

    fig, ax = plt.subplots()
    ax.plot(lowTrigThrs_run_mean[mask], lt_trigger_rate[mask], "o")
    ax.set_xlabel("mean trigger threshold / ADC")
    ax.set_ylabel("LT trigger rate / Hz")


    fig, ax = plt.subplots()

    mask = data["duration"] > 7000

    lt_trigger_rate = data["lt_trigger_rate"]
    plotting.pretty_trans_hist(ax, lt_trigger_rate[mask], np.linspace(0, 1, 20), color="k")

    ax.set_xlabel("LT trigger rate / Hz")
    fig.tight_layout()
    plt.show()