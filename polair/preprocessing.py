"""
titile: preprocessing.py
author: Laura Köhler
institution: Alfred-Wegener-Institut, Bremerhaven, Germany
contact: laura.koehler@awi.de
date: 2026-04-17
content: preprocessing command to convert all data to physical, calibrated values
comment: part of polair package
"""

from . import _helpers as h
from . import _calibration as calibration 
import numpy as np
import xarray as xr
import pandas as pd

def configure_preprocessing_parser(parser):
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
        "-v",
        "--verbose",
        metavar="DEBUG",
        help='Set the level of verbosity [DEBUG, INFO," " WARNING, ERROR]',
        required=False,
        default="INFO",
    )

    parser.set_defaults(func=run)
    
def run(args):
    # Preparation like importing config file etc.
    flight = int(args.flight)
    config_file = args.config
    config = h.import_dictionary(config_file)
    logfile = h.create_logfile(config)
    vars = h.import_dictionary(config["paths"]["variables"])
    cal_file = h.import_dictionary(config["paths"]["calibration"])
    fn_prefix = f"{config["flights"][flight]["data_dir"]}/{config["flights"][flight]["prefix"]}"
    outdir = config["paths"]["outdir"]
    flight_date = str(config["flights"][flight]["date"]).replace("-","")
    campaign = config["campaign"]["name"]
    fn_out = outdir+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_calibrated_raw_data.nc"
    h.add2logfile(logfile, f"Preprocessing: {config["campaign"]["name"]}, flight {flight}, {flight_date}")
    
    # Calibration of the raw data and interpolation to 100 Hz
    var_list = list(vars.keys())
    if not config["flights"][flight]["tbird"]:
        var_list = [
            v for v in var_list
            if vars[v]["group"] not in ["tbird", "inat"]
        ]
    for v in var_list:
        old_name = vars[v]["old"]
        fn = f"{config["flights"][flight]["data_dir"]}/{config["flights"][flight]["prefix"]}{old_name}.dat"
        df = pd.read_csv(fn, header  = 4, sep = r'\s+', names = ["date", "time", f"{v}"])
        df = h.get_timestamps(df)
        df = calibration.cal(v,cal_file,df, fn_prefix, vars)
        h.check_sampling(df, v, vars, logfile)
        h.find_gaps(df, v, vars, logfile)
        ds = h.interpolate_time(df, v, vars)
        ds = h.convert_unit(ds, vars, v)
        try:
            data_100Hz = xr.merge([data_100Hz, ds])
        except:
            data_100Hz = ds

    g_ratio = h.g_welmec(data_100Hz.lat_gprmc, data_100Hz.h_gpgga)/9.81
    for v in var_list:
        if str(vars[v]["units_old"])[:4] == "9.81":
            data_100Hz[v] = g_ratio * data_100Hz[v]
        data_100Hz = h.add_attrs_var(data_100Hz, v, vars)

    data_100Hz = h.add_global_attrs(data_100Hz, config, flight)
    data_100Hz.to_netcdf(fn_out)