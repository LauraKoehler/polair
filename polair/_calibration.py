import numpy as np
import pandas as pd
import xarray as xr

def cal(v, cal_file, df, fn_prefix, var_dict):
    """
    This function calibrates the DMS data and calculates physical values from analog output when necessary.

    Parameters:
    - v: str
        Variable to be calibrated
    - cal_file: dict
        Dictionary with calibration values for the campaign
    - df: pandas.DataFrame
        Dateframe with data to be calibrated
    - fn_prefix: str
        Filename prefix for loading the data
    - var_dict: dict
        Dictionary with variable information

    Returns:
    - df: pandas.DataFrame
        Dataframe with calibrated values of the variable v.
    """
    if v in ["psT", "psB", "psN"]:
        calibrated=np.interp(df[v],cal_file[v]["TRANSDUCER_OUTPUT"],cal_file[v]["APPLIED_PRESSURE"])
    elif v in ["qaT", "qbT", "qcT", "qcB", "qaN", "qbN", "qcN"]:
        calibrated=np.interp(df[v],cal_file[v]["TRANSDUCER_OUTPUT"],cal_file[v]["APPLIED_PRESSURE"]) * 2.4884
    elif v in ["Te_T", "TejB", "ThuB", "rFHuB", "pssB", "TejN", "Te_N", "ThuN", "rFHuN", "pssN", 
               "Pyau", "Pyao", "Pygu", "Pygo", 
               "axf", "ayf", "azf", "pitr", "rolr", "yawr", "pit", "roll", "thdg", "vew", "vns", "w", "lon", "lat", "h", "azg", "ttrk", "gs"]:
        calibrated=cal_file[v]["a0"] + cal_file[v]["a1"] * df[v] + cal_file[v]["a2"] * df[v]**2
        if v in ["Pygu", "Pygo"]:
            vt = f"T_{v}"
            old_name = var_dict[vt]["old"]
            fn_t = f"{fn_prefix}{old_name}.dat"
            df_t = pd.read_csv(fn_t, header  = 4, sep = r'\s+', names = ["date", "time", f"{v}"])
            temps=(np.sqrt(cal_file[vt]["a"]**2 - 4*cal_file[vt]["b"] * (-df_t[v]/100 + 1)) - cal_file[vt]["a"])/(2*cal_file[vt]["b"]) + 273.15
            calibrated = calibrated + cal_file[vt]["c"] * temps**4
    elif v in ["ttrk_inat", "mtrk_inat", "gs_inat", "t_inat_gpgga", "lat_inat", "lon_inat", "q_inat", "n_inat", "hdop_inat", "h_inat", "geoid_inat", 
               "t_inat_piahs", "gpsi_inat", "status_inat", "thdg_inat", "roll_inat", "pitch_inat",
               "t_gpgga", "n_gpgga", "hdop_gpgga", "h_gpgga", "geoid_gpgga", "t_gprmc", "gs_gprmc", "ttrk_gprmc", "lat_gprmc", "lon_gprmc",
               "gs_bestvel", "ttrk_bestvel", "w_bestvel", "age_bestvel", "latency_bestvel"]:
        message_cols = {"t_inat_gpgga": 1, "lat_inat": 2, "lon_inat": 4, "q_inat": 6, "n_inat": 7, "hdop_inat": 8, "h_inat": 9, "geoid_inat": 11,
                    "t_inat_piahs": 1, "gpsi_inat": 11, "status_inat": 12, "thdg_inat": 2, "roll_inat": 4, "pitch_inat": 5,
                     "t_gpgga": 1, "n_gpgga": 7, "hdop_gpgga": 8, "h_gpgga": 9, "geoid_gpgga": 11, 
                     "ttrk_inat": 1, "mtrk_inat": 3, "gs_inat": 7,
                     "t_gprmc": 1, "gs_gprmc": 7, "ttrk_gprmc": 8,  "lat_gprmc": 3,  "lon_gprmc": 5, 
                     "gs_bestvel": 13, "ttrk_bestvel": 14, "w_bestvel": 15, "age_bestvel": 12, "latency_bestvel": 11}
        if v in ["t_inat_gpgga", "lat_inat", "lon_inat", "q_inat", "n_inat", "hdop_inat", "h_inat", "geoid_inat", 
               "t_gpgga", "n_gpgga", "hdop_gpgga", "h_gpgga", "geoid_gpgga"]:
            df = df[df[v].str.startswith("$GPGGA")].reset_index()
        elif v in ["ttrk_inat", "mtrk_inat", "gs_inat"]:
            df = df[df[v].str.startswith("$GPVTG")].reset_index()
        elif v in ["t_inat_piahs", "gpsi_inat", "status_inat", "thdg_inat", "roll_inat", "pitch_inat"]:
            df = df[df[v].str.startswith("$PIAHS")].reset_index()
        elif v in ["t_gprmc", "gs_gprmc", "ttrk_gprmc", "lat_gprmc", "lon_gprmc"]:
            df = df[df[v].str.startswith("$GPRMC")].reset_index()
        elif v in ["gs_bestvel", "ttrk_bestvel", "w_bestvel", "age_bestvel", "latency_bestvel"]:
            df = df[df[v].str.startswith("#BESTVELA")].reset_index()
        split_cols = df[v].str.split(",", expand=True)
        vals = pd.to_numeric(split_cols[message_cols[v]], errors="coerce").values
        if v in ["t_inat_gpgga", "t_gpgga", "t_gprmc", "t_inat_piahs"]:
            dates = df["time"].values.astype('datetime64[D]').astype("str")
            dates_s = pd.Series(dates)
            times = split_cols[message_cols[v]].values
            times_s = pd.Series(times)
            calibrated = pd.to_datetime(dates_s + ' ' + times_s, format='%Y-%m-%d %H%M%S.%f').values
        elif v in ["lat_inat", "lat_gprmc"]:
            degs = (vals/100).astype("int")
            mins = vals - degs * 100
            vals = degs + mins/60
            orientation = split_cols[message_cols[v]+1].values
            signs = np.where(orientation == 'N', 1, -1)
            calibrated = signs * vals
        elif v in ["lon_inat", "lon_gprmc"]:
            degs = (vals/100).astype("int")
            mins = vals - degs * 100
            vals = degs + mins/60
            orientation = split_cols[message_cols[v]+1].values
            signs = np.where(orientation == 'E', 1, -1)
            calibrated = signs * vals
        else:
            calibrated = vals
    elif v in ["h_rad"]: 
        cond1 = (~(df[v]>10.40077).values).astype("int")
        cond2 = (df[v]>10.40077).values.astype("int")
        calibrated=(cond1 * (cal_file[v]["condition1"]["a0"] + cal_file[v]["condition1"]["a1"] * df[v] + cal_file[v]["condition1"]["a2"] * df[v]**2) + 
                    cond2 * (cal_file[v]["condition2"]["a0"] + cal_file[v]["condition2"]["a1"] * df[v] + cal_file[v]["condition2"]["a2"] * df[v]**2))
    elif v in ["T_Pyau", "T_Pyao", "T_Pygu", "T_Pygo"]:
        calibrated=(np.sqrt(cal_file[v]["a"]**2 - 4*cal_file[v]["b"] * (-df[v]/100 + 1)) - cal_file[v]["a"])/(2*cal_file[v]["b"]) + 273.15
    
    df = pd.DataFrame({"time": df["time"], v: calibrated})
        
    return df