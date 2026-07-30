"""
Definitions needed for aircraft noseboom and tbird processing.

Part of the polair package.

| title: _corr_fcts.py
| author: Laura Köhler
| institution: Alfred-Wegener-Institut, Bremerhaven, Germany
| contact: laura.koehler@awi.de
| date: 2026-04-17
"""

import numpy as np
import pandas as pd
import xarray as xr
import scipy.signal as sig
from scipy.ndimage import binary_dilation
from . import _helpers as h


def sat_correction(ds, ds_corr, t, recovery=1.0):
    """
    Compute static air temperature from TAT using adiabatic correction.

    Recovery is a correction for deiced sensor (recovery=1.00025).

    Args:
        ds: xarray.Dataset
            Dataset with all variables
        ds_corr: xarray.Dataset
            Dataset with corrected variables (not used in docstring but present in code)
        t: str
            Temperature variable name
        recovery: float, optional
            Recovery factor for deiced sensor (default: 1.00025)

    Returns:
        xarray.DataArray: Corrected temperature
    """
    ps_dict = {"Te_T": "ps", "TejB": "psB", "ThuB": "psB", "Te_N": "ps", "ThuN": "ps", "TejN": "ps"} # corresponding static pressures
    qs_dict = {"Te_T": "qc", "TejB": "qcB", "ThuB": "qcB", "Te_N": "qc", "ThuN": "qc", "TejN": "qc"} # corresponding dynamic pressures
    R_over_cp = 0.2858964
    temp = ds[t]
    if t in ["Te_N", "ThuN", "TejN", "Te_T"]:
        ps = ds_corr[ps_dict[t]]
        qs = ds_corr[qs_dict[t]]
    else:
        ps = ds[ps_dict[t]]
        qs = ds[qs_dict[t]]
    da = recovery * temp * (ps / (ps + qs)) ** R_over_cp
    return da


def sat_pressure(temp):
    """
    Calculate saturation pressure using the Magnus formula.

    Args:
        temp: xarray.DataArray
            Temperature in Kelvin

    Returns:
        xarray.DataArray: Saturation pressure in hPa
    """
    es = 6.1094 * np.exp(17.625 * (temp - 273.15)/(temp - 273.15 + 243.04))
    return es


def humidity_correction(rh, T_sensor, T_amb):
    """
    Apply adiabatic correction to relative humidity.

    Cuts values larger than 1.0 (limits of adiabatic correction).

    Args:
        rh: xarray.DataArray
            Relative humidity from humicap
        T_sensor: xarray.DataArray
            Humidity sensor temperature in K
        T_amb: xarray.DataArray
            Ambient temperature in K

    Returns:
        xarray.DataArray: Corrected relative humidity (capped at 1.0)
    """
    es_sensor = sat_pressure(T_sensor)
    es_amb = sat_pressure(T_amb)
    out = rh * (es_sensor/es_amb)
    out = out.where(out <= 1.0, 1.0)
    return out


def reverse_antennas(ds, angle, shift):
    """
    Apply antenna switching correction.

    If shift=True: shifts angle by 180° (possible reason: switched antennas in iNAT).

    Args:
        ds: xarray.Dataset
            Calibrated data
        angle: str
            Angle variable to be switched (e.g., "roll_inat", "pitch_inat", "thdg")
        shift: bool
            Whether to apply 180° shift (True) or keep original (False)

    Returns:
        xarray.DataArray: Shifted data if shift=True, else unchanged data
    """
    if shift:
        delta = 180
        sign = -1
    else:
        delta = 0
        sign = 1
    if angle in ["roll_inat", "pitch_inat"]:
        da = sign * ds[angle]
    else:
        da = (ds[angle] + delta) % 360
    return da


def get_w_ins(data, start, stop, deltat=0.01):
    """
    Calculate vertical velocity from INS vertical acceleration and remove Schuler oscillation.

    Args:
        data: xarray.Dataset
            Calibrated data
        start: numpy.datetime64
            Start time of the flight (from config)
        stop: numpy.datetime64
            End time of the flight (from config)
        deltat: float, optional
            Sampling rate in seconds (default: 0.01 for 100 Hz data)

    Returns:
        xarray.Dataset: Vertical velocity from INS with Schuler oscillation removed
    """
    w_ins = deltat * data["azg"].sel(time=slice(start, stop)).cumsum(dim="time")
    # Remove Schuler oscillation
    fs = 1/deltat  # Hz (sampling rate)
    fc = 1/(84*60)  # Hz (cutoff for Schuler oscillation)
    b, a = sig.butter(N=2, Wn=fc/(fs/2), btype='highpass')
    w_ins_hp = xr.apply_ufunc(sig.filtfilt, b, a, w_ins)
    w_ins_hp = w_ins_hp.to_dataset(name="w_ins")
    return w_ins_hp


def get_h_ins(w, deltat=0.01):
    """
    Calculate aircraft altitude from vertical acceleration and velocity.

    Args:
        w: xarray.Dataset
            Vertical velocity dataset
        deltat: float, optional
            Sampling rate in seconds (default: 0.01 for 100 Hz data)

    Returns:
        xarray.Dataset: Aircraft altitude from INS
    """
    h = (w * deltat).cumsum()
    h = h.to_dataset(name="h_ins")
    return h


def correct_ins_with_gps(data, v):
    """
    Stabilize INS data using GPS measurements.

    Args:
        data: xarray.Dataset
            100 Hz calibrated data
        v: str
            Variable to correct (options: "lon", "lat", "gs", "h_ins", "w_ins", "vew", "vns")

    Returns:
        xarray.Dataset: GPS-corrected data
    """
    gps_var = {
        "lon": "lon_gprmc",
        "lat": "lat_gprmc",
        "gs": "gs_bestvel",
        "h_ins": "h_gpgga",
        "w_ins": "w_bestvel",
        "vew": "vew_gps",
        "vns": "vns_gps"
    }
    data["vew_gps"] = data["gs_bestvel"] * np.sin(np.deg2rad(data["ttrk_bestvel"]))
    data["vns_gps"] = data["gs_bestvel"] * np.cos(np.deg2rad(data["ttrk_bestvel"]))
    gps_v = gps_var[v]

    rolling_ins = data[v].rolling(time=1000, center=True).mean()
    rolling_gps = data[gps_v].rolling(time=1000, center=True).mean()

    difference = rolling_ins - rolling_gps
    corrected = data[v] - difference
    corrected = corrected.to_dataset(name=f"{v}_corr")
    return corrected


def correct_ttrk_ins_with_gps(data, data_corr, v):
    """
    True heading correction from INS by GPS

    Args:
        data: xarray.Dataset
            100 Hz calibrated data
        data_corr: xarray.Dataset
            GPS-corrected data (from correct_ins_with_gps)
        v: str
            Variable to correct (only "ttrk" supported)

    Returns:
        xarray.Dataset: GPS-corrected true track
    """
    gps_var = {"ttrk": "ttrk_bestvel"}
    gps_v = gps_var[v]

    diffsin = np.sin(np.deg2rad(data[v])) - np.sin(np.deg2rad(data[gps_v]))
    diffcos = np.cos(np.deg2rad(data[v])) - np.cos(np.deg2rad(data[gps_v]))

    # For small speeds, we put the difference to zero
    diffsin = diffsin.where(data_corr["gs_corr"] > 30, other=0)
    diffcos = diffcos.where(data_corr["gs_corr"] > 30, other=0)

    rolling_sin = diffsin.rolling(time=1000, center=True).mean()
    rolling_cos = diffcos.rolling(time=1000, center=True).mean()

    sin_corr = np.sin(np.deg2rad(data[v])) - rolling_sin
    cos_corr = np.cos(np.deg2rad(data[v])) - rolling_cos

    corrected = (np.rad2deg(np.arctan2(-sin_corr, -cos_corr)) + 180) % 360
    corrected = corrected.to_dataset(name=f"{v}_corr")
    return corrected


def alignement_correction(data, fhp_params, platform, twist_angle):
    """
    Apply alignment corrections from mounting of the noseboom/t-bird.

    Parameters are determined from calibration segments with manual evaluation.

    Args:
        data: xarray.Dataset
            100 Hz calibrated data
        fhp_params: dict
            Dictionary with parameters for the five-hole probes
        platform: str
            "noseboom" or "tbird"
        twist_angle: float
            Rotation angle of the sonde (in degrees, from config file)

    Returns:
        xarray.DataArray: Dataset with corrected values
    """
    a0 = fhp_params[platform]["a0"]
    a1_qb = fhp_params[platform]["a1_qb"]
    a1_qc = fhp_params[platform]["a1_qc"]
    a1_ps = fhp_params[platform]["a1_ps"]
    a1_qratio = fhp_params[platform]["a1_qratio"]

    if platform == "noseboom":
        if v in ["qc", "ps"]:
            out = a0 + a1_qb * data.qbN + a1_qc * data.qcN + a1_ps * data.psN
        elif v in ["qb"]:
            out = np.cos(twist_angle) * data.qbN - np.sin(twist_angle) * data.qaN
        elif v in ["alpha"]:
            out = a0 + a1_qratio * (np.cos(twist_angle) * data.qaN + np.sin(twist_angle) * data.qbN)/data.qcN
        elif v in ["beta"]:
            b0 = fhp_params[platform]["qb"]["a0"]
            b1_qb = fhp_params[platform]["qb"]["a1_qb"]
            b1_qc = fhp_params[platform]["qb"]["a1_qc"]
            b1_ps = fhp_params[platform]["qb"]["a1_ps"]
            qb = b0 + b1_qb * data.qbN + b1_qc * data.qcN + b1_ps * data.psN
            out = a0 + a1_qratio * qb/data.qcN
    elif platform == "tbird":
        if v in ["qb", "qc", "ps"]:
            out = a0 + a1_qb * data.qbT + a1_qc * data.qcT + a1_ps * data.psT
        elif v in ["alpha"]:
            out = a0 + a1_qratio * data.qaT/data.qcT
        elif v in ["beta"]:
             out = a0 + a1_qratio * data.qbT/data.qcT
#            b0 = fhp_params[platform]["qb"]["a0"]
#            b1_qb = fhp_params[platform]["qb"]["a1_qb"]
#            b1_qc = fhp_params[platform]["qb"]["a1_qc"]
#            b1_ps = fhp_params[platform]["qb"]["a1_ps"]
#            qb = b0 + b1_qb * data.qbT + b1_qc * data.qcT + b1_ps * data.psT
#            out = a0 + a1_qratio * qb/data.qcT
    return out


def get_true_air_speed(data, platform):
    """
    Calculate true airspeed from air density.

    Args:
        data: xarray.Dataset
            Data with corrected variables (adiabatic corrected Te_N_corr and ps)
        platform: str
            "noseboom" or "tbird"

    Returns:
        xarray.DataArray: True airspeed
    """
    if platform == "noseboom":
        temp = "Te_N_corr"
        pres = "ps"
    elif platform == "tbird":
        temp = "Te_T_corr"
        pres = "ps"
    Rs = 287.0528
    rho = (data[pres]) / (Rs * data[temp])
    tas = np.sqrt(2 * data.qc/rho)
    return tas


def true_track_xarray(lat1, lon1, lat2, lon2):
    """
    Calculate the true track (bearing) between two geographic points.

    Args:
        lat1: xarray.DataArray
            Latitude at start point (degrees)
        lon1: xarray.DataArray
            Longitude at start point (degrees)
        lat2: xarray.DataArray
            Latitude at end point (degrees)
        lon2: xarray.DataArray
            Longitude at end point (degrees)

    Returns:
        xarray.DataArray: True track (bearing) in degrees [0, 360)
    """
    phi1 = np.deg2rad(lat1)
    phi2 = np.deg2rad(lat2)
    dl = np.deg2rad(lon2 - lon1)

    x = np.sin(dl) * np.cos(phi2)
    y = np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(dl)

    bearing = np.rad2deg(np.arctan2(x, y))
    bearing = (bearing + 360) % 360
    return bearing


def unwrap_with_nans(da, period=360):
    """
    Unwrap angle data while preserving NaN positions.

    Interpolates over NaNs to enable unwrapping, then restores original NaNs.

    Args:
        da: xarray.DataArray
            Angle data to be unwrapped (degrees)
        period: float, optional
            Period for unwrapping (default: 360 for degrees)

    Returns:
        xarray.DataArray: Unwrapped data array with NaNs at original positions
    """
    values = da.values.copy()
    nans = ~np.isfinite(values)
    # Interpolate over NaNs
    x = np.arange(len(values))
    values[nans] = np.interp(x[nans], x[~nans], values[~nans])
    # Unwrap
    unwrapped = np.unwrap(values, period=period)
    # Restore NaNs
    unwrapped[nans] = np.nan
    out = xr.DataArray(unwrapped, coords=da.coords, dims=da.dims)
    return out


def correct_ttrk_inat_with_gps(data, data_corr):
    """
    Correct INS true track using GPS with unwrapping.

    Args:
        data: xarray.Dataset
            100 Hz calibrated data
        data_corr: xarray.Dataset
            Data including ttrk with switched antenna correction

    Returns:
        xarray.Dataset: GPS-corrected INAT true track
    """
#    v= "ttrk_inat_corr"
#    lat = data.lat_inat
#    lon = data.lon_inat
#    gps_v = true_track_xarray(lat, lon, lat.shift(time=-1), lon.shift(time=-1))

#    rolling_inat = data_corr[v].rolling(time = 1000, center = True).mean()
#    rolling_gps = gps_v.rolling(time = 1000, center = True).mean()
    
#    difference = rolling_inat - rolling_gps
    
#    corrected = (data_corr[v] - difference) % 360
#    corrected = corrected.to_dataset(name = f"{v}")
    v = "ttrk_inat_corr"
    lat = data.lat_inat
    lon = data.lon_inat
    gps_v = true_track_xarray(lat, lon, lat.shift(time=-1), lon.shift(time=-1))
    ttrk_fixed = data_corr.ttrk_inat_corr

    # Unwrap safely and apply rolling mean
    inat_unwrap = unwrap_with_nans(ttrk_fixed)
    gps_unwrap = unwrap_with_nans(gps_v)

    inat_da = xr.DataArray(inat_unwrap, coords=data_corr[v].coords, dims=data_corr[v].dims)
    gps_da = xr.DataArray(gps_unwrap, coords=gps_v.coords, dims=gps_v.dims)

    rolling_inat = inat_da.rolling(time=1000, center=True).mean()
    rolling_gps = gps_da.rolling(time=1000, center=True).mean()

    difference = (rolling_inat - rolling_gps + 180) % 360 - 180

    corrected = (ttrk_fixed - difference) % 360
    corrected = corrected.to_dataset(name=v)
    return corrected


def angle_diff(a, b):
    """
    Calculate the shortest angle difference between two angles to determine peaks.

    Args:
        a: xarray.DataArray
            Angle a (radians)
        b: xarray.DataArray
            Angle b (radians)

    Returns:
        xarray.DataArray: Shortest angle difference in degrees [-180, 180)
    """
    diff = a - b
    angle_diff = np.rad2deg(np.arctan2(np.sin(diff), np.cos(diff)))
    return angle_diff


def mask_ttrk_thdg(ttrk, thdg):
    """
    Mask regions where true track and true heading differ significantly.

    Also masks 2-second windows around unphysical peaks in curves.

    Args:
        ttrk: xarray.DataArray
            True track (radians)
        thdg: xarray.DataArray
            True heading (radians)

    Returns:
        tuple: (masked ttrk, masked thdg) in radians
    """
    ttrk_thdg_diff = angle_diff(ttrk, thdg)
    valid = ttrk.notnull() & thdg.notnull()
    mask = valid & (np.abs(ttrk_thdg_diff) < 30)
    dttrk = angle_diff(ttrk, ttrk.shift(time=100))
    mask_curve = np.abs(dttrk) < 15
    structure = np.ones(201)
    mask_curve = xr.DataArray(
        ~binary_dilation(~mask_curve.values, structure=structure),
        dims=mask_curve.dims,
        coords=mask_curve.coords
    )
    mask_final = mask & mask_curve
    ttrk = ttrk.where(mask_final)
    thdg = thdg.where(mask_final)
    return ttrk, thdg


def get_wind_component(data, data_corr, component, platform):
    """
    Calculate wind components from calibrated raw data and corrected data.

    Args:
        data: xarray.Dataset
            Dataset with raw data
        data_corr: xarray.Dataset
            Dataset with corrected data
        component: str
            Wind component ("u", "v", or "vertwind")
        platform: str
            "noseboom" or "tbird"

    Returns:
        xarray.DataArray: Wind component
    """
    if platform == "noseboom":
        theta = np.deg2rad(data["pit"])
        phi = np.deg2rad(data["roll"])
        alpha = np.deg2rad(data_corr["alpha"])
        beta = np.deg2rad(data_corr["beta"])
        thdg = np.deg2rad(data["thdg"])
        c1, c2, c3 = 1.65, -0.41, 7.34
        vrxf = np.deg2rad(c1 * data["pitr"] - c2 * data["yawr"])
        vryf = np.deg2rad(c3 * data["yawr"] - c1 * data["rolr"])
        vrzf = np.deg2rad(c2 * data["rolr"] - c3 * data["pitr"])
        vns = data_corr["vns_corr"]
        vew = data_corr["vew_corr"]
        vup = data_corr["w_ins_corr"]
    elif platform == "tbird":
        theta = np.deg2rad(data_corr["pitch_inat_corr"])
        theta_rate = theta.diff("time")/0.01
        phi = np.deg2rad(data_corr["roll_inat_corr"])
        phi_rate = phi.diff("time")/0.01
        alpha = np.deg2rad(data_corr["alpha"])
        beta = np.deg2rad(data_corr["beta"])
        ttrk = np.deg2rad(data_corr["ttrk_inat_corr"])
        thdg = np.deg2rad(data_corr["thdg_inat_corr"])
        psi_rate = -(thdg - thdg.shift(time=1))/0.01
        c1, c2, c3 = 0.0, 0.0, 0.0
        vrxf = c1 * theta_rate - c2 * psi_rate
        vryf = c3 * psi_rate - c1 * phi_rate
        vrzf = c2 * phi_rate - c3 * theta_rate
        vns = data["gs_inat"] * np.cos(ttrk)
        vew = data["gs_inat"] * np.sin(ttrk)
        vup = data["h_inat"].rolling(time=100, center=True).mean().diff("time")/0.01

    # Calculate ground-relative wind components
    uKg = (vew
           + vrxf * np.cos(theta) * np.sin(thdg)
           + vryf * (np.sin(phi) * np.sin(theta) * np.sin(thdg) + np.cos(phi) * np.cos(thdg))
           + vrzf * (np.cos(phi) * np.sin(theta) * np.sin(thdg) - np.sin(phi) * np.cos(thdg))
        )
    vKg = (vns +
           vrxf * np.cos(theta) * np.cos(thdg)
           + vryf * (np.sin(phi) * np.sin(theta) * np.cos(thdg) - np.cos(phi) * np.sin(thdg))
           + vrzf * (np.cos(phi) * np.sin(theta) * np.cos(thdg) + np.sin(phi) * np.sin(thdg))
          )
    wKg = (vup
           + vrxf * np.sin(theta)
           - vryf * np.sin(phi) * np.cos(theta)
           - vrzf * np.cos(phi) * np.cos(theta)
          )

    # Calculate air-relative wind components
    ug = (data_corr["tas"] *
              (np.cos(alpha) * np.cos(beta) * np.cos(theta) * np.sin(thdg)
               + np.sin(beta) * (np.sin(phi) * np.sin(theta) * np.sin(thdg) + np.cos(phi) * np.cos(thdg))
               + np.sin(alpha) * np.cos(beta) * (np.cos(phi) * np.sin(theta) * np.sin(thdg) - np.sin(phi) * np.cos(thdg))
             ))
    vg = (data_corr["tas"] *
                (np.cos(alpha) * np.cos(beta) * np.cos(theta) * np.cos(thdg)
              + np.sin(beta) * (np.sin(phi) * np.sin(theta) * np.cos(thdg) - np.cos(phi) * np.sin(thdg))
              + np.sin(alpha) * np.cos(beta) * (np.cos(phi) * np.sin(theta) * np.cos(thdg) + np.sin(phi) * np.sin(thdg))
                ))
    wg = -(data_corr["tas"] *
              (-np.cos(alpha) * np.cos(beta) * np.sin(theta)
               + np.sin(beta) * np.sin(phi) * np.cos(theta)
               + np.sin(alpha) * np.cos(beta) * np.cos(phi) * np.cos(theta)
             ))

    if component == "u":
        out = uKg - ug
    elif component == "v":
        out = vKg - vg
    elif component == "vertwind":
        out = wKg - wg
        out = out - out.mean()
    return out

def mask_out_peaks(ds, refvar="p_amb", threshold=5000, timedelta=1):
    """
    Mask out peaks in data based on pressure changes.

    Args:
        ds: xarray.Dataset
            Data
        refvar: str, optional
            Variable used to check for peaks (default: "p_amb")
        threshold: float, optional
            Threshold to identify peaks (default: 5000 Pa)
        timedelta: int, optional
            Time window in seconds for peak detection (default: 1 s)

    Returns:
        xarray.Dataset: Data with peaks removed
    """
    mask = ds[refvar] - ds[refvar].shift(time=timedelta) < threshold
    structure = np.ones(11)
    mask = xr.DataArray(
        ~binary_dilation(~mask.values, structure=structure),
        dims=mask.dims,
        coords=mask.coords
    )
    ds = ds.where(mask)
    return ds

def check_flow(ds, refvar="flow_rate", variance=0.1):
    """
    Remove flow anomalies where flow rate varies more than specified percentage from average.

    Args:
        ds: xarray.Dataset
            Data
        refvar: str, optional
            Variable used to check flow rate (default: "flow_rate")
        variance: float, optional
            Allowed deviation from mean (default: 0.1 = 10%)

    Returns:
        xarray.Dataset: Data with flow anomalies removed
    """
    mean = ds[refvar].mean()
    mask = np.abs(ds[refvar] - mean) < variance * mean
    structure = np.ones(11)
    mask = xr.DataArray(
        ~binary_dilation(~mask.values, structure=structure),
        dims=mask.dims,
        coords=mask.coords
    )
    ds = ds.where(mask)
    return ds

def ampbox2swr_pyranometer(I):
    """
    Convert pyranometer current to shortwave radiation.

    Based on ampbox manual: 1 mV input → 1 mA output, so 4-20 mA represents 0-16 mV.

    Args:
        I: xarray.DataArray
            Raw current in A

    Returns:
        xarray.DataArray: Radiation in W/m²
    """
    # DMS: "4-20mA correspond to -50 W/m² to + 1950 W/m²"
    I_min, I_max = 4e-3, 20e-3    # A
    R_min, R_max = -50, 1950    # W/m2

    # Linear interpolation current → voltage (V)
    R = (I - I_min) / (I_max - I_min) * (R_max - R_min) + R_min
    return R

def ampbox2lwr_pyrgeometer(I, T):
    """
    Convert pyrgeometer current to longwave radiation.

    Includes temperature correction using Stefan-Boltzmann law.

    Args:
        I: xarray.DataArray
            Raw current in A
        T: xarray.DataArray
            Body temperature in K

    Returns:
        xarray.DataArray: Longwave radiation in W/m²
    """
    sigma = 5.670374 * 10**(-8)

    # DMS: "4-20mA correspond to -400 W/m² to +400 W/m²"
    I_min, I_max = 4e-3, 20e-3    # A
    R_min, R_max = -400, 400    # W/m2

    # Linear interpolation current → power flux (W/m2)
    R = (I - I_min) / (I_max - I_min) * (R_max - R_min) + R_min

    # Temperature correction
    G = R + sigma * T**4
    return G

def resistance2temperature(R):
    """
    Convert PT-100 resistance in Ohm in radiation sensors to temperature in K using Callendar-Van-Dusen equation ((https://de.wikipedia.org/wiki/Callendar-Van-Dusen-Gleichung), a, b from Datasheet.

    Args:
        R: xarray.DataArray
            Resistance in Ohm

    Returns:
        xarray.DataArray: Temperature in K
    """
    a = 3.9080 * 10**(-3)
    b = -5.8019 * 10**(-7)
    t = (-a + np.sqrt(a**2 - 4 * b * (-R/100 + 1)))/(2 * b) + 273.15
    return t

def get_radiation(config, flight, out_vars):
    """
    Calculate body temperatures and radiation from DMS raw data.

    Args:
        config: dict
            Configuration dictionary
        flight: int
            Flight number
        out_vars: dict
            Dictionary with output variables

    Returns:
        xarray.Dataset: Dataset with longwave and shortwave radiation
    """
    fn_prefix = f"{config["flights"][flight]["data_dir"]}/{config["flights"][flight]["prefix"]}"
    vars = list(reversed(np.sort(list(out_vars.keys()))))

    for v in vars:
        old_name = out_vars[v]["old"]
        fn = fn_prefix + old_name + ".dat"
        ds = h.import_radiation_data(fn, old_name)
        if v[0] == "t":
            data = resistance2temperature(ds[old_name]).to_dataset(name=old_name)
        else:
            if old_name[:6] == "Pyrano":
                data = ampbox2swr_pyranometer(ds[old_name]).to_dataset(name=old_name)
            if old_name[:6] == "Pyrgeo":
                t_name = old_name[:-7] + "T_raw"
                data = ampbox2lwr_pyrgeometer(ds[old_name], out[t_name]).to_dataset(name=old_name)
        try:
            out = xr.merge([out, data])
        except:
            out = data

    return out

def stp_conditions(ds, temp="t_amb", pres="p_amb"):
    """
    Add variables reduced to standard temperature and pressure (STP).
    The dataset needs to have the device internal temperature and pressure

    STP conditions: p₀ = 1013 hPa, T₀ = 0°C.

    Args:
        ds: xarray.Dataset
            Dataset with variables to reduce to STP, ambient temperature and pressure
        temp: str, optional
            Name of temperature variable (default: "t_amb")
        pres: str, optional
            Name of pressure variable (default: "p_amb")

    Returns:
        xarray.Dataset: Dataset with additional STP-corrected variables
    """
    stp_vars = {
        "number_conc": {
            'units': '1/m^3',
            'long_name': 'stp number concentration',
            'comment': 'stp conditions: 273.15 K, 1013 hPa'
        }
    }

    for v in ds.keys():
        if v in list(stp_vars.keys()):
            stp_corr = ds[v] * 101300/273.15 * ds[temp]/ds[pres]
            ds[f"{v}_stp"] = stp_corr
            if v == "number_conc":
                ds[f"{v}_stp"].attrs = stp_vars[v]
    return ds