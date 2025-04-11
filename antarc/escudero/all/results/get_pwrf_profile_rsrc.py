#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:35:22 2024

@author: prowe
"""


from netCDF4 import Dataset
import datetime as dt
import numpy as np
from os.path import exists


# My modules
from antarc.pwrf_rsrc import get_inds_to_latlon

# Parameters
from antarc import params


def get_temp_from_pottemppert(pot_temp_perts, pressures):
    """
    Get the temperature from the potential temperature perturbation
    as a function of time and pressure
    """
    temp = np.nan + np.ones(np.shape(pot_temp_perts))
    for i, press in enumerate(pressures):
        temp[:, i] = get_temp_from_pottemppert_at_level(
            pot_temp_perts[:, i], press
        )
    return temp


def get_temp_from_pottemppert_at_level(pot_temp_pert, pressure):
    """
    Total pot. temp. in K = T + 300 (T is the perturbation pot. temp.)
    temp = pot. temp. * (p/p_0)^kappa

      p_0 = 1000 mbar
      p is the pressure in mbar

      kappa is the Poisson constant (kappa = R/c_p), the ratio of the gas
      constant R to the specific heat at constant pressure c_p. For dry air
      kappa = 0.2854.
    """

    # Get and plot the near-surface temperature
    press0 = 1000.0  # mbar
    kappa = 0.2854
    pot_temp = pot_temp_pert + 273.15 + 300
    temp = pot_temp * (pressure / press0) ** kappa
    return temp


def get_pwrf_profiles_at_dtime(
    pwrf_temp_file, pwrf_rh_file, lat, lon, dtime, press
):
    """
    # Directories and files (dirs with "base_dir" require /year)
    pwrf_dir = f"{measdir}Escudero/pwrf/T/"
    rh_dir = f"{measdir}Escudero/pwrf/rh/"
    pwrf_dir2 = f"{measdir}Escudero/pwrf/T2/"

    # Formats and files
    fmt_pwrf = "wrfout_T_levels_d03_%Y%m%d.nc"
    rh_fmt = "wrfout_RH_levels_d03_%Y%m%d.nc"

    files_pwrf = ["wrfout_T_levels_d03_20220510.nc"]
    rh_files = ["wrfout_RH_levels_d03_20220510.nc"]
    pwrf_dates = [dt.datetime.strptime(x, fmt_pwrf) for x in files_pwrf]
    """

    if not exists(pwrf_temp_file) or not exists(pwrf_rh_file):
        return {
            "press": np.nan * np.ones(len(press)),
            "temp": np.nan * np.ones(len(press)),
            "rhw": np.nan * np.ones(len(press)),
            "hasdata": False,
        }

    # Inputs from parameter files
    measdir = params.MEAS_DIR
    latlon_dir = f"{measdir}Escudero/pwrf/T/eta_levels/"
    latlon_file = latlon_dir + "wrfout_T_d03_20220604.nc"

    # Get the lat and lon
    with Dataset(latlon_file) as nc:
        xlons = nc.variables["XLONG"][:]
        xlats = nc.variables["XLAT"][:]

    # Get polar wrf temperature profiles
    # T:description = "perturbation potential temperature theta-t0" ;
    # T2:description = "TEMP at 2 M" ;

    ilat, ilon = get_inds_to_latlon(xlats, xlons, lat, lon)

    # Check out relative humidities
    with Dataset(pwrf_rh_file) as nc:
        pwrf_press = nc["level"][:].data
        if np.any(pwrf_press != press):
            raise ValueError("Pressures differ!")

        # nc.variables.keys()
        # dict_keys(['time', 'level', 'RH_levels'])

        # Get all temps at closest grid point to Escudero
        rhw = nc["RH_levels"][:, :, ilat, ilon].data
        rhw[rhw == nc["RH_levels"]._FillValue] = np.nan

        if np.max(rhw) > 100:
            print("pause here")

    # Check out temperatures
    with Dataset(pwrf_temp_file) as nc:
        pwrf_press = nc["level"][:].data
        if np.any(pwrf_press != press):
            raise ValueError("Pressures differ!")

        # nc.variables.keys()
        # dict_keys(['time', 'level', 'T_levels'])

        # Get all temps at closest grid point to Escudero
        pot_temp_pert = nc["T_levels"][:, :, ilat, ilon]
        pot_temp_pert = pot_temp_pert.data
        pot_temp_pert[pot_temp_pert == nc["T_levels"]._FillValue] = np.nan
        temp = get_temp_from_pottemppert(pot_temp_pert, press)

    # Interpolate to dtime
    pwrf_time = [
        dt.datetime(dtime.year, dtime.month, dtime.day, x) for x in range(24)
    ]

    ddt = [np.abs(x - dtime) for x in pwrf_time]
    imin = np.argmin(ddt)
    if ddt[imin] > dt.timedelta(hours=1):
        raise ValueError("Time difference is too large!")

    pwrf = {
        "time": pwrf_time[imin],
        "temp": temp[imin, :],
        "rhw": rhw[imin, :],
        "hasdata": True,
    }

    # Missing: 'alt', 'ozone', 'uswind', 'vwind'
    return pwrf
