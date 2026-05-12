"""
titile: _helpers.py
author: Laura Köhler
institution: Alfred-Wegener-Institut, Bremerhaven, Germany
contact: laura.koehler@awi.de
date: 2026-04-17
content: defintions for DMS data processing
comment: part of polair package
"""

import numpy as np
import pandas as pd
import xarray as xr
import yaml
import re
from pathlib import Path
from datetime import datetime
import pint_xarray
import os

def import_dictionary(yaml_file):
    """
    Import config file for respective campaign and other dictionaries
    
    The config file contains:
    - campaign information
    - Paths to xml-file, data

    Parameters:
    - yaml_file: .yaml file 
        Dictionary containing basic information about the campaign

    Returns:
    - dict
        The config file as dictionary
    """
    try:
        with open(yaml_file, "r") as f:
            config = yaml.safe_load(f)
        return config
    except:
        print(f"{yaml_file}: yaml file not found or import error!!!")

def add2logfile(logfile, text):
    """
    Writes text in the logfile

    Parameters:
    - logfile: .txt file
        Logfile for the processing, location defined in the config file
    - text: str
        Text which should be added to the logfile
    """
    with logfile.open("a") as f:
        f.write(text+"\n")

def create_logfile(config):
    """
    Check if processing log file already exists, if not, create it, print date and time
    The config file needs to contain the path information of the log file.

    Paramters:
    - config: dict
        Configuration dictionary

    Returns:
    - logfile: str
        Path to logfile
    """
    fn = config["paths"]["processing_log_file"]+"turbulence_processing_log.txt"
    logfile = Path(fn)
    # Create the file if it doesn't exist
    if not logfile.exists():
        logfile.touch()
    # Append current date and time
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    add2logfile(logfile, f"{now}: processing starting")
    return logfile

def get_variable_names(xml_file):
    """
    Read a DMS order XML file and extract the variable names created by DMS.
    
    Background:
    - DMS filenames depend on the device configuration during specific flights or campaigns.
    - This function reads an XML order file (e.g., 'RAD_flightname.xml') and extracts
      variable names based on <channel> tags that include 'deviceShortName' and 'channelShortName'.

    Paramters:
    - xml_file: str
        Location of the xml file, should be included in the config

    Returns:
    - variable_names: list
        List with all variable names in the xml file.
    """
    try:
        # Read the entire file as text
        with open(xml_file, 'r', encoding='utf-8') as f:
            complete_file = f.read()
        # Split into non-empty, stripped lines
        lines = [line.strip() for line in complete_file.splitlines() if line.strip()]
        # Keep only lines containing 'channelShortName'
        lines = [line for line in lines if 'channelShortName' in line]
        variable_names = []
        
        # Extract device and channel short names from each line
        for line in lines:
            channel_match = re.search(r'channelShortName="([^"]+)"', line)
            device_match = re.search(r'deviceShortName="([^"]+)"', line)
        
            if channel_match and device_match:
                channel = channel_match.group(1)
                device = device_match.group(1)
                variable_names.append(f"{device}.{channel}")
        return variable_names
    except:
        print("No xml file found!!!")

def import_data(v, config, flight):
    """
    Import file for single variable using path and prefix for the specific flight from config file

    Parameters:
    - v: str
        Variable name
    - config: dict
        Configuration dictionary
    - flight: int
        Flight number

    Returns:
    - df: pandas.DataFrame
        Dataframe with time and variable data.
    """
    old_name = vars[v]["old"]
    fn = f"{config["flights"][flight]["data_dir"]}/{config["flights"][flight]["prefix"]}{old_name}.dat"
    df = pd.read_csv(fn, header  = 4, sep = r'\s+', names = ["date", "time", f"{v}"])
    return df

def get_timestamps(df):
    """
    Repair timestamps to datetime64 datetimes

    Parameters:
    - df: pandas.DataFrame
        Dataframe with times and data imported from the DMS download.

    Returns:
    - df: pandas.DataFrame
        Dataframe with the same data but datetime64 timestamps.
    """
    dates_fixed = np.array([d.replace('/', '-') for d in df["date"]], dtype=object)
    times = (dates_fixed+"T" +df["time"].values).astype("datetime64[ns]")
    df = pd.DataFrame({"time": times, df.keys()[-1]: df[df.keys()[-1]]})
    return df

def check_sampling(df, v, var_dict, logfile):
    """
    Check if sampling interval in df['time'] matches expected_sampling (seconds).
    95% of samples should match within tolerance. Identified gaps and sampling issues are written in the logfile.

    Parameters:
    - df: pandas.DataFrame
        Dataframe with data
    - v: str
        Variable name
    - var_dict: dict
        Dictionary with varible information
    - logfile: str
        Logfile path
    """
    expected_sampling = var_dict[v]["raw_sampling_time"]
    dt = df["time"].diff().dt.total_seconds().iloc[1:]
    dt_rounded = (dt * 100).round() / 100  # round to 0.01 s

    tolerance = 1e-5
    matches = abs(dt_rounded - expected_sampling) < tolerance
    percentage = matches.mean() * 100

    low = dt_rounded.quantile(0.025)
    high = dt_rounded.quantile(0.975)

    if percentage < 95 or low != expected_sampling or high != expected_sampling:
        print(
            f"{v}: ⚠️ only {percentage:.1f}% match; bounds: {low:.5f}–{high:.5f}s "
        )
        with logfile.open("a") as f:
            f.write(f"{v}: only {percentage:.1f}% match; bounds: {low:.5f}–{high:.5f}s")
    else:
        print(f"{v}: ✅ Sampling OK: {percentage:.1f}% within tolerance")

def find_gaps(df, v, var_dict, logfile, gap_factor=2.0):
    """
    Find data gaps in a time series based on sampling interval.
    
    Parameters
    - df : pandas.DataFrame
        Must contain a 'time' column of datetime64.
    - v: str
        Variable of the dataframe.
    - var_dict : dictionary
        Inlcudes raw_sampling_time for each variable.
    - gap_factor : float, optional
        Threshold factor (default 2.0).
    
    Returns
    - pandas.DataFrame
        Each row is a gap interval: [start_time, end_time, gap_duration].
    """
    expected_sampling = var_dict[v]["raw_sampling_time"]
    # Compute time differences in seconds
    dt = df["time"].diff().dt.total_seconds().iloc[1:]

    # Identify gaps (where sampling interval too long)
    gap_indices = dt[dt > gap_factor * expected_sampling].index

    # Collect gap intervals
    gaps = []
    for idx in gap_indices:
        start_time = df.loc[idx - 1, "time"]
        end_time = df.loc[idx, "time"]
        duration = (end_time - start_time).total_seconds()
        gaps.append((start_time, end_time, duration))

    # Return as DataFrame for easier inspection
    gaps_df = pd.DataFrame(gaps, columns=["start_time", "end_time", "gap_duration_s"])
    if len(gaps_df) > 0:
        add2logfile(logfile, f"{v}: gaps")
        for i in np.arange(len(gaps_df)):
            add2logfile(logfile, f"{i}: period: {gaps_df.iloc[0]["start_time"]} - {gaps_df.iloc[0]["end_time"]}, duration: {gaps_df.iloc[0]["gap_duration_s"]} s")
        
    return gaps_df

def interpolate_time(df, v, var_dict, steps = 0.01):
    """
    Linearly interpolate time on common timestamps with 100 Hz (or choose accordingly), convert to xarray dataset.

    Paramters:
    - df: pandas.DataFrame
        Dataframe containing data.
    - v: str
        Variable name
    - var_dict: dict
        Dictionary containing variable information. It should be specified in the config file.
    - steps: float
        Optional, frequency specification given in s, default is 100 Hz, i.e. 0.01 s

    Returns:
    - ds: xarray.Dataset
        Dataset interpolated to step frequency, default is 100 Hz.
    """
    if type(df[v].values[-1]) == np.datetime64:
        ds = xr.Dataset(coords = {"time": df["time"]}, data_vars = {v: ("time", df[v].values)})
        ds['sec_since1970'] = ds[v].astype('int64')
        ds = ds.resample(time=f"0.01s").interpolate("linear")
        ds[v] = ds['sec_since1970'].astype('datetime64[ns]')
        ds = ds[[v]]
    else:
        ds = df.set_index("time").to_xarray()
        if var_dict[v]["units_old"] == "degree" and v not in ["pitch_inat", "roll_inat", "rolr", "pitr", "pit", "roll", "lat", "lon", "lat_inat", "lon_inat", "lat_gprmc", "lon_gprmc"]: # Circular mean if unit in degrees
            sin = np.sin(np.deg2rad(ds[v])).resample(time=f"{steps}s").interpolate("linear")
            cos = np.cos(np.deg2rad(ds[v])).resample(time=f"{steps}s").interpolate("linear")
            ds = (np.rad2deg(np.arctan2(sin,cos)) % 360).to_dataset(name = v)
        else:
            ds = ds.resample(time=f"{steps}s").interpolate("linear")
    return ds

def add_attrs_var(ds, v, var_dict):
    """
    This function adds all attributes from the variable dictrionary except for the old name which is only necessary to read in the correct file. It adds the original unit because no unit conversion is done so far.

    Parameters:
    - ds: xarray.Dataset
        Input dataset
    - v: str 
        Variable for which attributes are added
    - var_dict: dict
        Dictionary with all the information for each parameter

    Returns:
    - ds: xarray.Dataset
        Dataset with added attributes.
    """
    attrs = var_dict[v]
    for at in list(attrs.keys()):
        if not at in ["old", "units", "units_old", "platform"]:
            att = attrs[at]
            if att == None:
                att = "-"
            ds[v].attrs[at] = att
        unit = attrs["units"]
        if not unit == "UTC":
            ds[v].attrs["units"] = unit
    return ds

def add_global_attrs(ds, config, flight):
    """
    This function adds metadata to the data set as stated in the config file

    Parameters:
    - ds: xarray.Dataset
        Input data set
    - config: dict
        config file containing campaign information
    - flight: int
        research flight which is processed
    """
    attrs = config["metadata"]
    ds.attrs = attrs
    title = f"{config["campaign"]["name"]} RF{flight:02} {config["flights"][flight]["date"]}: Calibrated raw data"
    ds.attrs["title"] = title
    return ds

def convert_unit(ds, var_dict, v):
    """
    This function converts the unit from the original unit to the final (SI) unit specified in the variable information

    Parameters:
    - ds: xarray.Dataset
        Input dataset
    - var_dict: dict
        Dictionary with variable (unit) information
    - v: str
        Variable for which unit should be converted

    Returns:
    - ds: xarray.Dataset
        Dataset with converted units
    """
    if not v in ["t_inat_gpgga", "t_inat_piahs", "t_gpgga", "t_gprmc"]:
        unit_old = str(var_dict[v]["units_old"])
        unit_new = str(var_dict[v]["units"])
        if unit_old[:4] == "9.81":
            ds_unit_old = 9.81 * ds[v].pint.quantify("m/s^2")
        else:
            ds_unit_old = ds[v].pint.quantify(unit_old)
        ds_unit_new = ds_unit_old.pint.to(unit_new)
        ds = ds_unit_new.to_dataset().pint.dequantify()
    return ds

def g_welmec(lat, h):
    '''
    ratio of gravitational acceleration according to Welmec-Formula devided by 9.81

    Parameters:
    - lat: xarray.DataArray
        Latitude in degree
    - h: xarray.DataArray
        height above sea level in m

    Returns:
    - g_ratio: xarray.DataArray
        Dataarray with latitude and height dependent values of g.
    '''
    g_ratio = (9.780318 * (1 + 0.0053024 * np.sin(np.deg2rad(lat))**2 - 0.0000058 * np.sin(2*np.deg2rad(lat))**2) - 0.000003085 * h)
    return g_ratio

def get_global_attributes(ds, config, instrument, flight):
    '''
    assigns attributes to the data set according to the config file.

    Parameters:
    - ds: xarray.Dataset
        data set which should get attributes
    - config: dict
        configuration dictionary
    - instrument: str
        instrument as called in the condig file
    - flight: int
        flight number

    Returns:
    - ds: xarray.Dataset
        data set with attributes.
    '''
    try: 
        attributes = config["instrument_metadata"][instrument]
        ds.attrs = attributes
        ds.attrs["title"] = attributes["title"] + f"RF{flight:02} ({str(config["flights"][flight]["date"])})"
        ds.attrs["campaign"] = config["campaign"]["name"]
        ds.attrs["platform"] = config["campaign"]["platform"]
        return ds
    except:
        print("No instrument metadata in config file")
        return ds

def add_segment_coordinate(ds, config, flight):
    '''
    assigns segment coordinate to a data set

    Parameters:
    - ds: xarray.Dataset
        data set which is supposed to get the segment coordinate
    - config: dict
        configuration dictionary
    - flight: in
        flight number

    Returns:
    - ds: xarray.Dataset
        data set with segment coordinate
    '''
    try:
        seg_fn = config["paths"]["segments"]
        segments = import_dictionary(seg_fn)
        flight_segments = segments["flights"][flight]["segments"]
        
        times  = [np.datetime64(s["time"]) for s in flight_segments]
        labels = [s["segment"] for s in flight_segments]
        
        # For each timestep, find which segment it belongs to
        indices = np.searchsorted(times, ds.time.values, side="right") - 1
        indices = np.clip(indices, 0, len(labels) - 1)
        segment_list = np.array(labels)[indices]
        
        ds = ds.assign_coords(segment=("time", segment_list))
        ds["segment"].attrs = {"long_name": "segment of the research flight"}
        return ds
    except:
        print("No segment file available")
        return ds

def import_device_data(indir, pf, time_offset):
    '''
    import data from different devices

    Parameters:
    - indir: str
        input directory
    - pf: str
        platform
    - time_offset: int
        offset time in ms between device and noseboom to be defined in config file

    Returns:
    - ds: xarray.Dataset
        combined data from files in input directory
    '''
    if pf  == "mcpc":
        files = np.sort([f for f in os.listdir(indir) if f.endswith('.TXT')])
    elif pf == "partector":
        files = np.sort([f for f in os.listdir(indir) if f.endswith('.txt')])
    for fn in files:
        if pf == "mcpc":
            df = pd.read_csv(f"{indir}/{fn}", header = 13, sep = "\t")
            times = pd.to_datetime("20" + df["#YY/MM/DD"]+" "+df["HR:MN:SC"]) - np.timedelta64(time_offset, "ms")
        elif pf == "partector":
            df = pd.read_csv(f"{indir}/{fn}", header = 18, sep = "\t")
            with open(f"{indir}/{fn}") as file:
                date_str = [next(file) for x in range(9)][-1]
            clean = date_str.replace('Start: ', '').strip()
            start_time = np.datetime64(datetime.strptime(clean, '%d.%m.%Y %H:%M:%S'))
            seconds = df["time"].values.astype("timedelta64[s]")
            times = start_time + seconds
        df["time"] = times
        try:
            df_all = pd.concat([df_all, df], ignore_index = True)
        except:
            df_all = df
    df_all = df_all.sort_values(by = "time").set_index("time")
    df_all = df_all[~df_all.index.duplicated(keep='last')]
    ds = df_all.to_xarray()
    return ds