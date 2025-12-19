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
    outdir = config["paths"]["outdir"]
    flight_date = str(config["flights"][flight]["date"]).replace("-","")
    campaign = config["campaign"]["name"]
    fn_in = outdir+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_calibrated_raw_data.nc"
    fn_out = outdir+"/"+campaign+"_"+flight_date+f"_RF{flight:02}_tbird_100Hz.nc"
    start = config["flights"][flight]["start_tbird"]
    stop = config["flights"][flight]["stop_tbird"]

    data = xr.open_dataset(fn_in)

    # adiabatic corrections for sensors in Rosemount/Goodrich housings
    etaE = 1.00025  # recovery factor for deiced sensors
    for temp in ["Te_T", "TejB", "ThuB", "Te_N", "ThuN", "TejN"]:
        name = f"{temp}_corr"
        if temp == "TejB":
            da = corr.sat_correction(data, temp, recovery=etaE).rename(name)
        else:
            da = corr.sat_correction(data, temp).rename(name)
        try:
            data_corr[name] = da
        except:
            data_corr = da.to_dataset()

    for a in ["thdg_inat", "ttrk_inat", "mtrk_inat", "roll_inat", "pitch_inat"]:
        name = f"{a}_corr"
        da = corr.reverse_antennas(data, a).rename(name)
        data_corr[name] = da

    data_corr["tang_track_inat"] = corr.bearing_from_latlon(data.lat_inat, data.lon_inat)

    pf = "tbird"

    for v in ["qb", "qc", "ps", "alpha", "beta"]:
        data_corr[v] = corr.alignement_correction(data, fhp_params, v, pf)
    
    # Some corrections:
    data_corr["qc"] = data_corr["qc"].where(data_corr["qc"]>0, other = 0)
    data_corr["alpha"] = data_corr["alpha"].where(data_corr["qc"]>500, other = 0)
    data_corr["beta"] = data_corr["beta"].where(data_corr["qc"]>500, other = 0)
    
    data_corr["tas"] = corr.get_true_air_speed(data_corr, pf)
    
    data_corr["lat_corr"] = data.lat_inat # This is only to avoid more complications in the definition of the function deltah
    
    deltah = corr.deltah(data_corr)
    data_corr["h_baro"] = deltah.cumsum("time", skipna=True)
    data_corr["w_baro"] = deltah/0.01
    
    for v in ["u", "v", "vertwind"]:
        wind_comp = corr.get_wind_component(data, data_corr, v, pf)
        data_corr[v] = wind_comp

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

    out_ds.to_netcdf(fn_out)