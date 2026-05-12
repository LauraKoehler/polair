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
        help="instrument to be processed, options are mcpc",
        default=None,
        required=True,
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
    pf = args.instrument
    out_vars = h.import_dictionary(config["paths"]["processed_variables"])
    out_vars = out_vars[pf]
    indir = config["flights"][flight][pf]
    outdir = config["paths"]["outdirs"][pf]
    flight_date = str(config["flights"][flight]["date"]).replace("-","")
    campaign = config["campaign"]["name"]
    fn_out = outdir+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_{pf}.nc"
    start = np.datetime64(config["flights"][flight]["start"])
    stop = np.datetime64(config["flights"][flight]["stop"])
    h.add2logfile(logfile, f"{pf}: {config["campaign"]["name"]}, flight {flight}, {flight_date}")

    try:
        time_offset = config["flights"][flight]["time_offsets"][pf]
    except:
        print("No time offset for this instrument mentioned in config.")
        time_offset = 0

    ds = h.import_device_data(indir, pf, time_offset)

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
        indir_nb = config["paths"]["outdirs"]["noseboom"]
        fn_nb = indir_nb+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_noseboom_100Hz.nc"
        data_nb = xr.open_dataset(fn_nb)
        data_nb_int = data_nb.interp({"time":out_ds.time.values}, method = "nearest")
        out_ds = out_ds.assign_coords({"lat": data_nb_int.lat, "lon": data_nb_int.lon, "alt": data_nb_int.alt})

    except:
        print("No processed noseboom file with position data found for this flight")

    # Cut first and last 2 minutes since there are often very high counts due to the own emissions
    if pf in ["mcpc"]:
        out_ds = out_ds.sel(time = slice(start + np.timedelta64(2,"m"), stop - np.timedelta64(2,"m")))
        out_ds = corr.mask_out_peaks(out_ds)
        out_ds = corr.check_flow(out_ds)

    out_ds = h.get_global_attributes(out_ds, config, pf, flight)
    out_ds = h.add_segment_coordinate(out_ds, config, flight)

    out_ds.to_netcdf(fn_out)