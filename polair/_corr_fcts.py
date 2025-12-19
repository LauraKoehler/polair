import numpy as np
import pandas as pd
import xarray as xr
import scipy.signal as sig
from . import _helpers as h

def sat_correction(ds, t, recovery=1.0):
    """
    Compute static air temperature from TAT using adiabatic correction.
    t, ps, qc must be xarray.DataArray.
    Recovery is a correction for deiced sensor, in this case recovery=1.00025

    ds: dataset with all variables
    t: temperature name (str)
    """
    ps_dict = {"Te_T": "psT", "TejB": "psB", "ThuB": "psB", "Te_N": "psN", "ThuN": "psN", "TejN": "psN"} # corresponding static pressures
    qs_dict = {"Te_T": "qcT", "TejB": "qcB", "ThuB": "qcB", "Te_N": "qcN", "ThuN": "qcN", "TejN": "qcN"} # corresponding dynamic pressures
    R_over_cp = 0.2858964
    temp = ds[t]
    ps = ds[ps_dict[t]]
    qs = ds[qs_dict[t]]
    da = recovery * temp * (ps / (ps + qs)) ** R_over_cp
    return da

def reverse_antennas(ds, angle, shift = True):
    '''
    if shift = True (default): shifts angle by pi/2, else keep angle 
    (possible reason: switched antennas in iNAT)
    
    ds: data set
    angle: 
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
        da = (ds[angle] + delta) % 360
    return da

def bearing_from_latlon(lat, lon):
    """
    Compute bearing between consecutive lat/lon points.
    lat, lon must be 1-D xarray DataArrays in degrees.
    Returns a DataArray of same length with first element NaN.
    """
    
    # Convert to radians
    lons = np.deg2rad(lon)
    lats = np.deg2rad(lat)

    # Differences
    dlon = lons.diff("time").values
    lat1 = lats.isel(time=slice(None, -1)).values
    lat2 = lats.isel(time=slice(1, None)).values

    # Bearing formula
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    
    bearing = np.arctan2(x, y)

    # Convert to degrees and normalize to [0, 360)
    bearing = (np.rad2deg(bearing) + 360) % 360
    bearing = np.append(bearing, np.array([np.nan]))

    da = xr.DataArray(data = bearing, dims=["time"], coords = {"time": (["time"], lat.time.values)})
    
    return da

def get_w_ins(data, start, stop, deltat = 0.01):
    '''
    Calculate w from vertial acceleration from the INS

    data: calibrated data
    start: start of the flight (from config)
    stop: end of the flight (from config)
    deltat: sampling rate in s (should be 0.01 for 100 Hz data)
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
    
    w: vertical velocity
    deltat: sampling rate in s (should be 0.01 for 100 Hz data)
    '''
    h = (w * deltat).cumsum()
    h = h.to_dataset(name = "h_ins")
    return h

def correct_ins_with_gps(data, v):
    '''
    INS stabilization with GPS

    data: 100 Hz data
    v: variable, options: lon, lat, gs, h_ins, w_ins, vew, vns
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

    data: 100 Hz calibrated data
    data_corr: GPS corrected data calculated with correct_ins_with_gps
    v: variable, options: ttrk
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

def alignement_correction(data, fhp_params, v, platform):
    '''
    Alignemnet corrections from mounting of the noseboom/t-bird

    data: 100 Hz calibrated data
    fhb_params: dictionary with parameters for the five hole probes
    platform: noseboom or tbird
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
            out = a0 + a1_qratio * data.qaN/data.qcN
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
            b0 = fhp_params[platform]["qb"]["a0"]
            b1_qb = fhp_params[platform]["qb"]["a1_qb"]
            b1_qc = fhp_params[platform]["qb"]["a1_qc"]
            b1_ps = fhp_params[platform]["qb"]["a1_ps"]
            qb = b0 + b1_qb * data.qbT + b1_qc * data.qcT + b1_ps * data.psT
            out = a0 + a1_qratio * qb/data.qcT
    return out

def get_true_air_speed(data, platform):
    '''
    Calculate true air speed from air density

    data: data with corrected variables (adiabatic corrected Te_N_corr and ps)
    '''
    if platform == "noseboom":
        temp = "Te_N_corr"
    elif platform == "tbird":
        temp = "Te_T_corr"
    Rs = 287.058
    rho = data.ps / (Rs * data[temp])
    tas = np.sqrt(2 * data.qc/rho)
    return tas

def hbaro_icao(data):
    '''
    Barometric height from ICAO standard atmosphere

    data: data with ps
    '''
    T0 = 288.15
    L = 0.0065
    p0 = 101315
    R = 287.05287
    g = 9.80665
    h = T0/L * ((data.ps/p0)**(-R * L / g) -1)
    return h

def hbaro(data, start, stop, platform):
    '''
    Get barometric pressure from pressure and temperature

    data: dataset with ps and Te_N_corr
    config: config file defining start and stop of the flight
    '''
    if platform == "noseboom":
        t_air = "Te_N_corr"
    elif platform == "tbird":
        t_air = "Te_T_corr"
    Rs = 287.058
    hicao = hbaro_icao(data)
    g = h.g_welmec(data.lat_corr, hicao)
    deltap = data.ps.diff("time")
    deltah = (- Rs * data[t_air] * deltap)/ (data.ps * g)
    hbaro = deltah.sel(time = slice(start, stop)).cumsum()
    return hbaro

def deltah(data):
    '''
    height differences from pressure. It includes a rolling mean over 10 s to reduce the noice from the ps sensors.

    data: dataset including ps and Te_N_corr
    '''
    Rs = 287.058
    hicao = hbaro_icao(data)
    g = h.g_welmec(data.lat_corr, hicao)
    p = data.ps.rolling(time = 1000).mean()
    deltah = (- Rs * data.Te_N_corr * p.diff("time"))/ (p * g)
    return deltah

def get_wind_component(data, data_corr, component, platform):
    '''
    Get wind components from calibrated raw data and corrected data

    data: dataset with raw data
    data_corr: dataset with corrected data
    component: wind component, options "u", "v", "vertwind"
    platform: noseboom or tbird
    '''
    if platform == "noseboom":
        theta = np.deg2rad(data["pit"])
        phi = np.deg2rad(data["roll"])
        alpha = np.deg2rad(data_corr["alpha"])
        beta = np.deg2rad(data_corr["beta"])
        thdg = np.deg2rad(data["thdg"])
    
        c1 = 1.65
        c2 = -0.41
        c3 = 7.34
        vrxf = np.deg2rad(c1 * data["pitr"] - c2 * data["yawr"])
        vryf = np.deg2rad(c3 * data["yawr"] - c1 * data["rolr"])
        vrzf = np.deg2rad(c2 * data["rolr"] - c3 * data["pitr"])

        vKg = (data_corr["vns_corr"] + 
               vrxf * np.cos(theta) * np.cos(thdg)
               + vryf * (np.sin(phi) * np.sin(theta) * np.cos(thdg) - np.cos(phi) * np.sin(thdg))
               + vrzf * (np.cos(phi) * np.sin(theta) * np.cos(thdg) + np.sin(phi) * np.sin(thdg))
              )
        uKg = (data_corr["vew_corr"]
               + vrxf * np.cos(theta) * np.sin(thdg)
               + vryf * (np.sin(phi) * np.sin(theta) * np.sin(thdg) + np.cos(phi) * np.cos(thdg))
               + vrzf * (np.cos(phi) * np.sin(theta) * np.sin(thdg) - np.sin(phi) * np.cos(thdg))
            )
        wKg = (-data_corr["w_ins_corr"]
               - vrzf * np.sin(theta)
               + vryf * np.sin(phi) * np.cos(theta)
               + vrzf * np.cos(phi) * np.cos(theta)
              )
        
    elif platform == "tbird":
        theta = np.deg2rad(data["pitch_inat"])
        phi = np.deg2rad(data["roll_inat"])
        alpha = np.deg2rad(data_corr["alpha"])
        beta = np.deg2rad(data_corr["beta"])
        thdg = np.deg2rad(data["thdg_inat"])
        ttrk = np.deg2rad(data_corr["tang_track_inat"])

        uKg = 0 #data["gs_inat"] * np.cos(ttrk)
        vKg = 0 #data["gs_inat"] * np.sin(ttrk)
        wKg = 0 #- data["h_inat"].diff("time")/0.01

    vg = (data_corr["tas"] *
                (np.cos(alpha) * np.cos(beta) * np.cos(theta) * np.cos(thdg)
              + np.sin(beta) * (np.sin(phi) * np.sin(theta) * np.cos(thdg) - np.cos(phi) * np.sin(thdg))
              + np.sin(alpha) * np.cos(beta) * (np.cos(phi) * np.sin(theta) * np.cos(thdg) + np.sin(phi) * np.sin (thdg))
                ))
    ug = (data_corr["tas"] *
              (np.cos(alpha) * np.cos(beta) * np.cos(theta) * np.sin(thdg)
               + np.sin(beta) * (np.sin(phi) * np.sin(theta) * np.sin(thdg) + np.cos(phi) * np.cos(thdg))
               + np.sin(alpha) * np.sin(beta) * (np.cos(phi) * np.sin(theta) * np.sin(thdg) - np.sin(phi) * np.cos(thdg))
             ))
    wg = (data_corr["tas"] *
              (-np.cos(alpha) * np.cos(beta) * np.sin(theta)
               + np.sin(beta) * np.sin(phi) * np.cos(theta)
               + np.sin(alpha) * np.cos(beta) * np.cos(phi) * np.cos(theta)
             ))

    if component == "v":       
        out = vKg - vg
    elif component == "u":
        out = uKg - ug
    elif component == "vertwind":
        out = -wKg + wg
    return out