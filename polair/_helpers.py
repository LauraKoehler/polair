import numpy as np
import pandas as pd
import xarray as xr
import yaml
import re
from pathlib import Path
from datetime import datetime

def import_dictionary(yaml_file):
    """
    Import config file for respective campaign and other dictionaries
    
    The config file contains:
    - campaign information
    - Paths to xml-file, data
    """
    try:
        with open(yaml_file, "r") as f:
            config = yaml.safe_load(f)
        return config
    except:
        print("yaml file not found!!!")

def add2logfile(logfile, text):
    with logfile.open("a") as f:
        f.write(text+"\n")

def create_logfile(config):
    """
    Check if processing log file already exists, if not, create it, print date and time
    The config file needs to contain the path information of the log file.
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
    """
    old_name = vars[v]["old"]
    fn = f"{config["flights"][flight]["data_dir"]}/{config["flights"][flight]["prefix"]}{old_name}.dat"
    df = pd.read_csv(fn, header  = 4, sep = r'\s+', names = ["date", "time", f"{v}"])
    return df

def get_timestamps(df):
    """
    Repair timestamps to datetime64 datetimes
    """
    dates_fixed = np.array([d.replace('/', '-') for d in df["date"]], dtype=object)
    times = (dates_fixed+"T" +df["time"].values).astype("datetime64[ns]")
    df = pd.DataFrame({"time": times, df.keys()[-1]: df[df.keys()[-1]]})
    return df

def check_sampling(df, v, var_dict, logfile):
    """
    Check if sampling interval in df['time'] matches expected_sampling (seconds).
    95% of samples should match within tolerance.
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
    ----------
    df : pandas.DataFrame
        Must contain a 'time' column of datetime64.
    v: variable
        Variable of the dataframe.
    var_dict : dictionary
        Inlcudes raw_sampling_time for each variable.
    gap_factor : float, optional
        Threshold factor (default 2.0).
    
    Returns
    -------
    pandas.DataFrame
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
    Linearly interpolate time on common timestamps with 100 Hz (or choose accordingly), convert to xarray dataset

    df: pandas dataframe
    steps: optional
        Frequency specification given in s, default is 100 Hz, i.e. 0.01 s
    """
    if type(df[v].values[-1]) == np.datetime64:
        ds = xr.Dataset(coords = {"time": df["time"]}, data_vars = {v: ("time", df[v].values)})
        ds['sec_since1970'] = ds[v].astype('int64')
        ds = ds.resample(time=f"0.01s").interpolate("linear")
        ds[v] = ds['sec_since1970'].astype('datetime64[ns]')
        ds = ds[[v]]
    else:
        ds = df.set_index("time").to_xarray()
        if var_dict[v]["units_old"] == "degree": # Circular mean if unit in degrees
            sin = np.sin(np.deg2rad(ds[v])).resample(time=f"{steps}s").interpolate("linear")
            cos = np.cos(np.deg2rad(ds[v])).resample(time=f"{steps}s").interpolate("linear")
            ds = (np.rad2deg(np.arctan2(sin,cos)) % 360).to_dataset(name = v)
        else:
            ds = ds.resample(time=f"{steps}s").interpolate("linear")
    return ds