"""
Processing command for dropsondes

This module is part of the polair package and provides functionality for
processing the dropsondes launched during the entire campaign.

| title: dropsondes.py
| author: Laura Köhler
| institution: Alfred-Wegener-Institut, Bremerhaven, Germany
| contact: laura.koehler@awi.de
| date: 2026-10-02

The dropsonde command processes dropsonde data for the entire campaign.

.. code-block:: bash

    polair dropsondes -c <config file>


"""

from . import _helpers as h
from . import _corr_fcts as corr 
import numpy as np
import xarray as xr
import pandas as pd
import os

def configure_dropsonde_parser(parser):
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
    dev = "dropsondes"
    config_file = args.config
    config = h.import_dictionary(config_file)
    logfile = h.create_logfile(config)

    outdir = config["paths"]["outdirs"][dev]
    campaign = config["campaign"]["name"]
    fn_out = outdir+"/"+campaign+f"_{dev}.nc"

    out_vars = h.import_dictionary(config["paths"]["processed_variables"])
    out_vars = out_vars[dev]

    ds = h.import_dropsondes(config)

    var_list = list(out_vars.keys())
    out_ds = None
    
    for v in var_list:
        v_old = out_vars[v]["old"]
        var_data = ds[[v_old]]
        var_data = var_data.rename({v_old: v})
        var_data[v] = var_data[v].where(var_data[v]>-900, other = np.nan)
        var_data = h.convert_unit(var_data, out_vars, v)
        out_ds = var_data if out_ds is None else xr.merge([out_ds, var_data])
    
        out_ds[v].attrs = {}
        out_ds = h.add_attrs_var(out_ds, v, out_vars)

    out_ds = corr.mask_out_dropsonde_peaks(out_ds, out_vars)
    
    out_ds = h.get_global_attributes_without_flight(out_ds, config, dev)
    
    out_ds.to_netcdf(fn_out)