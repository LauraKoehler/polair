"""
titile: _corr_fcts.py
author: Laura Köhler
institution: Alfred-Wegener-Institut, Bremerhaven, Germany
contact: laura.koehler@awi.de
date: 2026-04-17
content: defintions needed for aircraft noseboom and tbird processing
comment: part of polair package
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
    Recovery is a correction for deiced sensor, in this case recovery=1.00025

    Paramters:
    - ds: xarray.Dataset
        dataset with all variables
    - t: str
        temperature name
        
    Returns:
    - da: xarray.DataArray
        DataArray with corrected temperature
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
    Saturation pressure from Magnus formula.

    Parameters:
    - temp: xarray.DataArray
        temperature

    Returns:
    - es: xarray.DataArray
        saturation pressure
    """
    es = 6.1094 * np.exp(17.625 * (temp - 273.15)/(temp - 273.15 + 243.04))
    return es

def humidity_correction(rh, T_sensor, T_amb):
    """
    Adiabatic correction of relative humidity, cut values larger than 1 (limits of adiabatic correction)

    Parameters:
    - rh: xarray.DataArray
        relative humidity from humicap
    - T_sensor: xarray.DataArray
        humidity sensor temperature
    - T_amb: xarray.DataArray
        ambient temperature

    Returns:
    - out: xarray.DataArray
        corrected relative humidity
    """
    es_sensor = sat_pressure(T_sensor)
    es_amb = sat_pressure(T_amb)
    out = rh * (es_sensor/es_amb)
    out = out.where(out <= 1.0, 1.0)
    return out

def reverse_antennas(ds, angle, shift):
    '''
    if shift = True: shifts angle by pi/2, else keep angle 
    (possible reason: switched antennas in iNAT)

    Parameters:
    - ds: xarray.Dataset
        calibrated data
    - angle: str
        angle to be switched
    - shift: bool
        options: True or False

    Returns:
    - da: xarray.DataArray
        shifted data if shift = True, else unchanged data
    '''
    if shift:
        delta = 180
        sign = -1
    else:
        delta = 0
        sign = 1
    if angle in ["roll_inat", "pitch_inat"]:
        da = sign * ds[angle]
    else:
        da = (ds[angle] + delta)
    return da

def get_w_ins(data, start, stop, deltat = 0.01):
    '''
    Calculate w from vertial acceleration from the INS and remove Schuler oscillation.

    Parameters:
    - data: xarrau.Dataset
        calibrated data
    - start: numpy.datetime64
        start of the flight (from config)
    - stop: numpy.datetime64
        end of the flight (from config)
    - deltat: float
        sampling rate in s (should be 0.01 for 100 Hz data)

    Returns:
    - w_ins_hp: xarray.Dataset
        Vertical velocity from INS
    '''
    w_ins = deltat * data["azg"].sel(time = slice(start, stop)).cumsum(dim = "time")
    # Remove Schuler oscillation
    # sampling rate
    fs = 1/deltat # Hz
    
    # cutoff for Schuler oscillation
    fc = 1/(84*60)  # Hz
    
    b, a = sig.butter(N=2, Wn=fc/(fs/2), btype='highpass')
    w_ins_hp = xr.apply_ufunc(sig.filtfilt, b, a, w_ins)
    w_ins_hp = w_ins_hp.to_dataset(name = "w_ins")
    return w_ins_hp

def get_h_ins(w, deltat = 0.01):
    '''
    Get the height from vertical acceleration and velocity

    Parameters:
    - w: xarray.Dataset
        vertical velocity
    - deltat: float
        sampling rate in s (should be 0.01 for 100 Hz data)

    Returns:
    - h: xarray.Dataset
        Aircraft altitude from INS
    '''
    h = (w * deltat).cumsum()
    h = h.to_dataset(name = "h_ins")
    return h

def correct_ins_with_gps(data, v):
    '''
    INS stabilization with GPS

    Parameters:
    - data: xarray.Dataset
        100 Hz calibrated data
    - v: str
        variable, options: lon, lat, gs, h_ins, w_ins, vew, vns

    Returns:
    - corrected: xarray.Dataset
        GPS corrected data
    '''
    gps_var = {"lon": "lon_gprmc",
              "lat": "lat_gprmc",
              "gs": "gs_bestvel",
              "h_ins": "h_gpgga",
              "w_ins": "w_bestvel",
              "vew": "vew_gps",
              "vns": "vns_gps"}
    data["vew_gps"] = data["gs_bestvel"] * np.sin(np.deg2rad(data["ttrk_bestvel"]))
    data["vns_gps"] = data["gs_bestvel"] * np.cos(np.deg2rad(data["ttrk_bestvel"]))
    gps_v = gps_var[v]

    rolling_ins = data[v].rolling(time = 1000, center = True).mean()
    rolling_gps = data[gps_v].rolling(time = 1000, center = True).mean()
    
    difference = rolling_ins - rolling_gps
    
    corrected = data[v] - difference
    corrected = corrected.to_dataset(name = f"{v}_corr")
    return corrected

def correct_ttrk_ins_with_gps(data, data_corr, v):
    '''
    True heading correction from INS by GPS

    Parameters:
    - data: xarray.Dataset
        100 Hz calibrated data
    - data_corr: xarray.Dataset
        GPS corrected data calculated with correct_ins_with_gps
    - v: str
        variable, options: ttrk

    Returns:
    - corrected: xarray.Dataset
        GPS corrected data
    '''
    gps_var = {"ttrk": "ttrk_bestvel"}
    gps_v = gps_var[v]

    diffsin = np.sin(np.deg2rad(data[v])) - np.sin(np.deg2rad(data[gps_v]))
    diffcos = np.cos(np.deg2rad(data[v])) - np.cos(np.deg2rad(data[gps_v]))

    # For small speeds, we put the difference to zero
    diffsin = diffsin.where(data_corr["gs_corr"] > 30, other = 0)
    diffcos = diffcos.where(data_corr["gs_corr"] > 30, other = 0)

    rolling_sin = diffsin.rolling(time = 1000, center = True).mean()
    rolling_cos = diffcos.rolling(time = 1000, center = True).mean()

    sin_corr = np.sin(np.deg2rad(data[v])) - rolling_sin
    cos_corr = np.cos(np.deg2rad(data[v])) - rolling_cos

    corrected = (np.rad2deg(np.arctan2(-sin_corr, -cos_corr)) + 180) % 360

    corrected = corrected.to_dataset(name = f"{v}_corr")
    return corrected

def alignement_correction(data, fhp_params, v, platform, twist_angle):
    '''
    Alignemnet corrections from mounting of the noseboom/t-bird. The used parameters are determined from the calibration segments with a manual evaluation.

    Paramters:
    - data: xarray.Dataset
        100 Hz calibrated data
    - fhb_params: dict
        dictionary with parameters for the five hole probes
    - platform: str
        noseboom or tbird
    - twist_angle: float
        rotation angle of the sonde, to be specified in the condig file

    Returns:
    - out: xarray.DataArray
        Dataset with corrected values.
    '''
    a0 = fhp_params[platform][v]["a0"]
    a1_qb = fhp_params[platform][v]["a1_qb"]
    a1_qc = fhp_params[platform][v]["a1_qc"]
    a1_ps = fhp_params[platform][v]["a1_ps"]
    a1_qratio = fhp_params[platform][v]["a1_qratio"]
    if platform == "noseboom":
        if v in ["qb", "qc", "ps"]:
            out = a0 + a1_qb * data.qbN + a1_qc * data.qcN + a1_ps * data.psN
        elif v in ["alpha"]:
            out = a0 + a1_qratio * (np.cos(twist_angle) * data.qaN + np.sin(twist_angle) * data.qbN)/data.qcN
        elif v in ["beta"]:
            out = a0 + a1_qratio * (np.cos(twist_angle) * data.qbN - np.sin(twist_angle) * data.qaN)/data.qcN
#            b0 = fhp_params[platform]["qb"]["a0"]
#            b1_qb = fhp_params[platform]["qb"]["a1_qb"]
#            b1_qc = fhp_params[platform]["qb"]["a1_qc"]
#            b1_ps = fhp_params[platform]["qb"]["a1_ps"]
#            qb = b0 + b1_qb * data.qbN + b1_qc * data.qcN + b1_ps * data.psN
#            out = a0 + a1_qratio * (np.cos(twist_angle) * qb - np.sin(twist_angle) * data.qaN)/data.qcN
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
    '''
    Calculate true air speed from air density

    Parameters:
    - data: xarray.Dataset
        data with corrected variables (adiabatic corrected Te_N_corr and ps)
    - platform: str
        options: noseboom or tbird

    Returns:
    - tas: xarray.DataArray
        true airspeed
    '''
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
    This calculates the true track from lat and lon.

    Parameters:
    - lat1: xarray.DataArray
        lat at start of time interval
    - lon1: xarray.DataArray
        lon at start of time interval
    - lat2: xarray.DataArray
        lat at end of time interval
    -lon 2: xarray.DataArray
        lon at end of time interval

    Returns:
    - bearing: xarray.DataArray
        The true track between 1 and 2
    """
    phi1 = np.deg2rad(lat1)
    phi2 = np.deg2rad(lat2)
    dl   = np.deg2rad(lon2 - lon1)

    x = np.sin(dl) * np.cos(phi2)
    y = np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(dl)

    bearing = np.rad2deg(np.arctan2(x, y))
    bearing = (bearing + 360) % 360
    return bearing

def correct_ttrk_inat_with_gps(data, data_corr):
    '''
    INS stabilization with GPS

    Parameters:
    - data: xarray.Dataset
        100 Hz calibrated data
    - data_corr: xarray.Dataset
        data including the ttrk with switched antenna correction

    Returns:
    - corrected: xarray.Dataset
        GPS corrected INAT ttrk
    '''
    v= "ttrk_inat_corr"
    lat = data.lat_inat
    lon = data.lon_inat
    gps_v = true_track_xarray(lat, lon, lat.shift(time=-1), lon.shift(time=-1))

    rolling_inat = data_corr[v].rolling(time = 1000, center = True).mean()
    rolling_gps = gps_v.rolling(time = 1000, center = True).mean()
    
    difference = rolling_inat - rolling_gps
    
    corrected = data_corr[v] - difference
    corrected = corrected.to_dataset(name = f"{v}")
    return corrected

def angle_diff(a, b):
    '''
    shorteset angle difference in degree to determine peaks

    Parameters:
    - a: xarray.DataArray
        angle a in rad
    - b: xarray.DataArray
        angle b in rad

    Returns:
    - angle_diff: xarray.DataArray
        shortest angle difference in degree
    '''
    diff = a - b
    angle_diff = np.rad2deg(np.arctan2(np.sin(diff), np.cos(diff)))
    return angle_diff

def mask_ttrk_thdg(ttrk, thdg):
    '''
    The shift between true heading and true track and very big differences in curves lead to unphysical peaks in the wind. Thus, this function masks regions out where ttrk and thdg differ too much and 2 s around the unphysical peaks.

    Parameters:
    - ttrk: xarray.DataArray
        true track in rad
    - thdg: xarray.DataArray
        true heading in rad

    Returns:
    - ttrk: xarray.DataArray
        masked true track in rad
    - thdg: xarray.DataArray
        masked true heading in rad
    '''
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
    '''
    Get wind components from calibrated raw data and corrected data

    Parameters:
    - data: xarray.Dataset
        dataset with raw data
    - data_corr: xarray.Dataset
        dataset with corrected data
    - component: str
        wind component, options "u", "v", "vertwind"
    - platform: str
        options: noseboom or tbird

    Returns:
    - out: xarray.DataArray
        wind component
    '''
    if platform == "noseboom":
        theta = np.deg2rad(data["pit"])
        phi = np.deg2rad(data["roll"])
        alpha = np.deg2rad(data_corr["alpha"])
        beta = np.deg2rad(data_corr["beta"])
        thdg = np.deg2rad(data["thdg"])

        # Difference between five hole probe and INS, five hole probe seems to be c2 m right, c1 m above, and c3 m in front of the INS (looking from the noseboom to the aircraft)
        c1 = 1.65
        c2 = -0.41
        c3 = 7.34
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
        ttrk, thdg = mask_ttrk_thdg(ttrk, thdg)
        psi_rate = -(thdg - thdg.shift(time = 1))/0.01

        c1 = 0.0
        c2 = 0.0
        c3 = 0.0
        vrxf = c1 * theta_rate - c2 * psi_rate
        vryf = c3 * psi_rate - c1 * phi_rate
        vrzf = c2 * phi_rate - c3 * theta_rate
        vns = data["gs_inat"] * np.cos(ttrk)
        vew = data["gs_inat"] * np.sin(ttrk)
        vup = data["h_inat"].diff("time")/0.01

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
    ug = (data_corr["tas"] *
              (np.cos(alpha) * np.cos(beta) * np.cos(theta) * np.sin(thdg)
               + np.sin(beta) * (np.sin(phi) * np.sin(theta) * np.sin(thdg) + np.cos(phi) * np.cos(thdg))
               + np.sin(alpha) * np.cos(beta) * (np.cos(phi) * np.sin(theta) * np.sin(thdg) - np.sin(phi) * np.cos(thdg))
             ))
    vg = (data_corr["tas"] *
                (np.cos(alpha) * np.cos(beta) * np.cos(theta) * np.cos(thdg)
              + np.sin(beta) * (np.sin(phi) * np.sin(theta) * np.cos(thdg) - np.cos(phi) * np.sin(thdg))
              + np.sin(alpha) * np.cos(beta) * (np.cos(phi) * np.sin(theta) * np.cos(thdg) + np.sin(phi) * np.sin (thdg))
                ))
    wg = -(data_corr["tas"] *
              (-np.cos(alpha) * np.cos(beta) * np.sin(theta)
               + np.sin(beta) * np.sin(phi) * np.cos(theta)
               + np.sin(alpha) * np.cos(beta) * np.cos(phi) * np.cos(theta)
             ))

    if component == "u":       
        out = uKg  - ug
    elif component == "v":
        out = vKg - vg
    elif component == "vertwind":
        out = wKg - wg
        out = out - out.mean()
    return out