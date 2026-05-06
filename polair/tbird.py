"""
titile: tbird.py
author: Laura Köhler
institution: Alfred-Wegener-Institut, Bremerhaven, Germany
contact: laura.koehler@awi.de
date: 2026-04-17
content: T-Bird processing command
comment: part of polair package
"""

from . import _helpers as h
from . import _corr_fcts as corr 
import numpy as np
import xarray as xr
import pandas as pd

def configure_tbird_parser(parser):
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
    indir = config["paths"]["outdirs"]["raw"]
    outdir = config["paths"]["outdirs"]["tbird"]
    flight_date = str(config["flights"][flight]["date"]).replace("-","")
    campaign = config["campaign"]["name"]
    fn_in = indir+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_calibrated_raw_data.nc"
    fn_out = outdir+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_tbird_100Hz.nc"
    start = config["flights"][flight]["start_tbird"]
    stop = config["flights"][flight]["stop_tbird"]
    h.add2logfile(logfile, f"T-Bird: {config["campaign"]["name"]}, flight {flight}, {flight_date}")

    data = xr.open_dataset(fn_in)

    # Calibration of the pressures. For this, the calibration segments are used and the calibration factors are saved in fhp_params. The calibration has to be performed manually. Use the jupyter notebook Get_segments.ipynb.
    pf = "tbird"

    twist_angle = 0
    for v in ["qb", "qc", "ps", "alpha", "beta"]:
        da = corr.alignement_correction(data, fhp_params, v, pf, twist_angle).rename(v)
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
    for temp in ["Te_T"]:
        name = f"{temp}_corr"
        da = corr.sat_correction(data, data_corr, temp).rename(name)
        data_corr[name] = da

    # In case of switched antennas (stated in config), reverse angles
    switched_antennas = config["campaign"]["tbird_reversed_antennas"]

    for a in ["thdg_inat", "roll_inat", "pitch_inat", "ttrk_inat"]:
        name = f"{a}_corr"
        da = corr.reverse_antennas(data, a, switched_antennas).rename(name)
        data_corr[name] = da
    
    data_corr["tas"] = corr.get_true_air_speed(data_corr, pf)

    # True track is GPS corrected using lat and lon from INAT which come from GPS. This is used as true heading in the wind computation.
    ttrk_corr = corr.correct_ttrk_inat_with_gps(data, data_corr)
    data_corr["ttrk_inat_corr"] = ttrk_corr.ttrk_inat_corr

    data = data.sel(time = slice(start,stop))
    data_corr = data_corr.sel(time = slice(start,stop))
    
    for v in ["u", "v", "vertwind"]:
        wind_comp = corr.get_wind_component(data, data_corr, v, pf)
        data_corr[v] = wind_comp

    # Cleaning up and prepare the output dataset
    out_vars = out_vars[pf]
    var_list = list(out_vars.keys())
    for v in var_list:
        v_old = out_vars[v]["old"]
        if out_vars[v]["platform"] in ["tbird"]:
            try:
                var_data = data[[v_old]]
            except:
                var_data = data_corr[[v_old]]
            var_data = var_data.sel(time = slice(start,stop)).rename({v_old: v})
            try:
                out_ds = xr.merge([out_ds, var_data])
            except:
                out_ds = var_data
            out_ds = h.add_attrs_var(out_ds, v, out_vars)

    out_ds = h.get_global_attributes(out_ds, config, pf, flight)
    out_ds = h.add_segment_coordinate(out_ds, config, flight)

    out_ds.to_netcdf(fn_out)