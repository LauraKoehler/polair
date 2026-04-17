"""
titile: noseboom.py
author: Laura Köhler
institution: Alfred-Wegener-Institut, Bremerhaven, Germany
contact: laura.koehler@awi.de
date: 2026-04-17
content: noseboom processing command
comment: part of polair package
"""

from . import _helpers as h
from . import _corr_fcts as corr 
import numpy as np
import xarray as xr
import pandas as pd

def configure_noseboom_parser(parser):
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
    flight = int(args.flight)
    config_file = args.config
    config = h.import_dictionary(config_file)
    logfile = h.create_logfile(config)
    vars = h.import_dictionary(config["paths"]["variables"])
    out_vars = h.import_dictionary(config["paths"]["processed_variables"])
    fhp_params = h.import_dictionary(config["paths"]["fiveholeprobe"])
    outdir = config["paths"]["outdir"]
    flight_date = str(config["flights"][flight]["date"]).replace("-","")
    campaign = config["campaign"]["name"]
    fn_in = outdir+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_calibrated_raw_data.nc"
    fn_out = outdir+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_noseboom_100Hz.nc"
    start = config["flights"][flight]["start"]
    stop = config["flights"][flight]["stop"]
    h.add2logfile(logfile, f"Noseboom: {config["campaign"]["name"]}, flight {flight}, {flight_date}")

    data = xr.open_dataset(fn_in)

    # Calibration of the pressures. For this, the calibration segments are used and the calibration factors are saved in fhp_params. The calibration has to be performed manually. Use the jupyter notebook Get_segments.ipynb.
    pf = "noseboom"

    for v in ["qb", "qc", "ps", "alpha", "beta"]:
        da = corr.alignement_correction(data, fhp_params, v, pf).rename(v)
        try:
            data_corr[v] = da
        except:
            data_corr = da.to_dataset()
    
    # Some corrections:
    data_corr["qc"] = data_corr["qc"].where(data_corr["qc"]>0, other = 0)
    data_corr["alpha"] = data_corr["alpha"].where(data_corr["qc"]>500, other = 0)
    data_corr["beta"] = data_corr["beta"].where(data_corr["qc"]>500, other = 0)
    
    # adiabatic corrections for sensors in Rosemount/Goodrich housings
    etaE = 1.00025  # recovery factor for deiced sensors
    for temp in ["TejB", "ThuB", "Te_N", "ThuN", "TejN"]:
        name = f"{temp}_corr"
        if temp == "TejB":
            da = corr.sat_correction(data, data_corr, temp, recovery=etaE).rename(name)
        else:
            da = corr.sat_correction(data, data_corr, temp).rename(name)
        try:
            data_corr[name] = da
        except:
            data_corr = da.to_dataset()

    # Correct relative humidity using the Magnus formula for the saturation pressure.
    rh_corr = corr.humidity_correction(data.rFHuN, data.ThuN, data_corr.Te_N_corr)
    data_corr["rFHuN_corr"] = rh_corr

    # Correction of INS with GPS
    w_ins = corr.get_w_ins(data, start, stop)
    h_ins = corr.get_h_ins(w_ins["w_ins"])
    data = xr.merge([data, w_ins, h_ins])
    
    for v in ["lat", "lon", "gs", "h_ins", "w_ins", "vew", "vns"]:
        gps_corr = corr.correct_ins_with_gps(data, v)
        data_corr = xr.merge([data_corr, gps_corr])

    ttrk_corr = corr.correct_ttrk_ins_with_gps(data, data_corr, "ttrk")
    data_corr = xr.merge([data_corr, ttrk_corr])

    # Calculate true airspeed
    data_corr["tas"] = corr.get_true_air_speed(data_corr, pf)

    # Determine wind components
    for v in ["u", "v", "vertwind"]:
        wind_comp = corr.get_wind_component(data, data_corr, v, pf)
        data_corr[v] = wind_comp

    # Cleaning up and prepare the output dataset
    out_vars = out_vars[pf]
    var_list = list(out_vars.keys())
    for v in var_list:
        v_old = out_vars[v]["old"]
        if out_vars[v]["platform"] in ["noseboom"]:
            try:
                var_data = data[[v_old]]
            except:
                var_data = data_corr[[v_old]]
            var_data = var_data.sel(time = slice(start,stop)).rename({v_old: v})
            try:
                out_ds = xr.merge([out_ds, var_data])
            except:
                out_ds = var_data
            out_ds[v].attrs = {}
            out_ds = h.add_attrs_var(out_ds, v, out_vars)

    out_ds.to_netcdf(fn_out)