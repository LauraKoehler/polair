"""
titile: device.py
author: Laura Köhler
institution: Alfred-Wegener-Institut, Bremerhaven, Germany
contact: laura.koehler@awi.de
date: 2026-05-08
content: processing command for different devices
comment: part of polair package
"""

from . import _helpers as h
from . import _corr_fcts as corr 
import numpy as np
import xarray as xr
import pandas as pd
import os

def configure_device_parser(parser):
    parser.add_argument(
        "-f",
        "--flight",
        help="Research flight number (integer)",
        default=None,
        required=True,
    )
    parser.add_argument(
        "-c",
        "--config",
        metavar="CONFIG_FILE",
        help="configuration file (yaml)",
        default=None,
        required=True,
    )
    parser.add_argument(
        "-i",
        "--instrument",
        help="instrument to be processed, options are mcpc, partector, partector_dms, kt19, radiation",
        default=None,
        required=True,
    )
    parser.add_argument(
        "-p",
        "--platform",
        help='platform where instrument is installed, options are polar6 or tbird',
        required=False,
        default="polar6",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        metavar="DEBUG",
        help='Set the level of verbosity [DEBUG, INFO," " WARNING, ERROR]',
        required=False,
        default="INFO",
    )

    parser.set_defaults(func=run)
    
def run(args):
    flight = int(args.flight)
    config_file = args.config
    config = h.import_dictionary(config_file)
    logfile = h.create_logfile(config)
    dev = args.instrument
    pf = args.platform
    if pf == "polar6":
        dev_pf = dev
    elif pf == "tbird":
        dev_pf = f"{dev}_{pf}"
    dev_name = dev
    if dev[-4:] == "_dms":
        dev_name = dev[:-4]
    if pf == "tbird":
        dev_name = f"{dev_name}_{pf}"
    out_vars = h.import_dictionary(config["paths"]["processed_variables"])
    out_vars = out_vars[dev]
    if dev in ["kt19", "radiation"]:
        indir = config["flights"][flight]["data_dir"]
    else:
        indir = config["flights"][flight][dev_name]
    outdir = config["paths"]["outdirs"][dev_name]
    flight_date = str(config["flights"][flight]["date"]).replace("-","")
    campaign = config["campaign"]["name"]
    fn_out = outdir+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_{dev_name}.nc"
    if pf == "polar6":
        start = np.datetime64(config["flights"][flight]["start"])
        stop = np.datetime64(config["flights"][flight]["stop"])
    elif pf == "tbird":
        start = np.datetime64(config["flights"][flight]["start_tbird"])
        stop = np.datetime64(config["flights"][flight]["stop_tbird"])
    h.add2logfile(logfile, f"{dev_name}: {config["campaign"]["name"]}, flight {flight}, {flight_date}")

    try:
        time_offset = config["flights"][flight]["time_offsets"][dev_name]
    except:
        print("No time offset for this instrument mentioned in config.")
        time_offset = 0

    if dev == "radiation":
        ds = corr.get_radiation(config, flight, out_vars)
    else:
        ds = h.import_device_data(indir, dev, time_offset)
    # Resampling to 1 sec time resolution is only done if resample = True. Otherwise, the original time stamps are kept.
    ds = h.resample2sec(ds, resample = False)

    var_list = list(out_vars.keys())
    
    for v in var_list:
        v_old = out_vars[v]["old"]
        var_data = ds[[v_old]]
        var_data = var_data.sel(time = slice(start,stop)).rename({v_old: v})
        var_data = h.convert_unit(var_data, out_vars, v)
        try:
            out_ds = xr.merge([out_ds, var_data])
        except:
            out_ds = var_data
    
        out_ds[v].attrs = {}
        out_ds = h.add_attrs_var(out_ds, v, out_vars)

    try:
        if pf == "polar6":
            indir_pf = config["paths"]["outdirs"]["noseboom"]
            fn_pf = indir_pf+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_noseboom_100Hz.nc"
        elif pf == "tbird":
            indir_pf = config["paths"]["outdirs"]["tbird"]
            fn_pf = indir_pf+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_tbird_100Hz.nc"
        data_pf = xr.open_dataset(fn_pf)
        data_pf_int = data_pf.interp({"time":out_ds.time.values}, method = "nearest")
        out_ds = out_ds.assign_coords({"lat": data_pf_int.lat, "lon": data_pf_int.lon, "alt": data_pf_int.alt})
        if dev in ["radiation"]:
            out_ds = out_ds.assign({"pitch": data_pf_int.pitch, "roll": data_pf_int["roll"]})

    except:
        print("No processed turbulence file with position data found for this flight")

    # Reduction to standard temperature and pressure for devices specified in config file
    devs4stp = config["campaign"]["stp"]
    if dev in devs4stp:
        out_ds = corr.stp_conditions(out_ds)
    
    if dev in ["mcpc", "partector"]:
        # Cut first and last 2 minutes since there are often very high counts due to the own emissions
        out_ds = out_ds.sel(time = slice(start + np.timedelta64(2,"m"), stop - np.timedelta64(2,"m")))
    if dev in ["mcpc", "partector", "partector_dms"]:
        out_ds = corr.mask_out_peaks(out_ds)
        out_ds = corr.check_flow(out_ds)

    out_ds = h.get_global_attributes(out_ds, config, dev_name, flight)
    out_ds = h.add_segment_coordinate(out_ds, config, flight)
    
    out_ds.to_netcdf(fn_out)