import mattak.Dataset
import time
import numpy
import argparse


def compare_event_info(ev1, ev2):

    if not numpy.all(list(ev1.__dict__.keys()) == list(ev2.__dict__.keys())):
        return False

    for key in ev1.__dict__:
        try:
            if not ev1.__dict__[key] == ev2.__dict__[key]:
                return False
        except ValueError:
            if not numpy.all(ev1.__dict__[key] == ev2.__dict__[key]):
                print(key, ev1.__dict__[key], ev2.__dict__[key])
                return False

    return True


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Test mattak.")
    argparser.add_argument('--station', type=int, default=0)
    argparser.add_argument('--run', type=int, default=0)
    argparser.add_argument('--data_dir', type=str, default=None)
    argparser.add_argument('--voltage_calibration', '-vc', type=str, default=None)
    argparser.add_argument('--calibrate', action="store_true")
    args = argparser.parse_args()

    calibrated = args.calibrate

    dset_pyroot = mattak.Dataset.Dataset(
        args.station, args.run, data_path=args.data_dir,
        backend="pyroot", verbose=True, read_daq_status=True,
        voltage_calibration=args.voltage_calibration)

    dset_uproot = mattak.Dataset.Dataset(
        args.station, args.run, data_path=args.data_dir,
        backend="uproot", verbose=True, read_daq_status=True,
        voltage_calibration=args.voltage_calibration)

    dset_pyroot.setEntries((0, 100))
    dset_uproot.setEntries((0, 100))

    wfs = dset_pyroot.wfs(calibrated=calibrated)
    print(wfs)

    if calibrated:
        # 0.0001 -> 0.1 mV (1 ADC -> 0.61mV,  0.6 / sqrt(12) -> 1.76)
        # at maximum of 700mV max tolerance is 0.1 mV + 0.001 * 700 mV = 0.17 mV
        # at minmim of -1300mV max tolerance is 0.1mV + 0.001 * 1300 mV = 0.23 mV
        assert numpy.allclose(wfs, dset_uproot.wfs(calibrated=calibrated), rtol=0.001, atol=0.0001), "Both backends return different waveforms."
    else:
        assert numpy.all(wfs == dset_uproot.wfs(calibrated=calibrated)), "Backends return different waveforms."


    eventInfos_pyroot = dset_pyroot.eventInfo()
    eventInfos_uproot = dset_uproot.eventInfo()
    for ev_py, ev_up in zip(eventInfos_pyroot, eventInfos_uproot):
        assert compare_event_info(ev_py, ev_up), "Backends return different event infos"

    print("Successful")