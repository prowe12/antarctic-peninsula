#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 13:36:00 2022

@author: prowe
"""

import cdsapi
import numpy as np
import datetime as dt

from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import esc202202

LATITUDE = esc_params.LATITUDE
LONGITUDE = esc_params.LONGITUDE
ALTITUDE = esc_params.ALTITUDE

DATE1 = esc202202.DATE1
DATE2 = esc202202.DATE2

# Back up to the 4th
DATE1 = dt.datetime(2022, 2, 2, 0, 0, tzinfo=dt.timezone.utc)


outdir = "/Users/prowe/Sync/measurements/Escudero/era5/"

north = np.round(LATITUDE / 0.25) * 0.25 + 0.25
south = north - 0.25
east = np.round(LONGITUDE / 0.25) * 0.25 + 0.25
west = east - 0.5

# Pressure levels
press_levs = np.hstack(
    [
        np.arange(1000, 925 - 5, -25),
        np.arange(900, 200, -50),
        np.arange(200, 100 - 5, -25),
    ]
)
press_levs_ints = list(press_levs) + [70, 50, 30, 20, 10, 7, 5, 3, 2]
press_levs = [str(p) for p in press_levs_ints]

# Times
times = [f"0{i}:00" for i in range(10)] + [f"{i}:00" for i in range(10, 24)]
download_flag = True

yearstr = str(DATE1.year)
monthstr = str(DATE1.month).zfill(2)
days = list(range(DATE1.day - 1, DATE2.day + 2))  # pad out a day on each side

# Get downwelling shortwave and longwave fluxes
dataset = "reanalysis-era5-single-levels"
c = cdsapi.Client()
for day in days[:3]:
    daystr = str(day).zfill(2)
    outfile = "era5_esc_broadband_202202" + daystr + ".nc"

    # api parameters:
    # area: [north, west, south, east]; west is negative
    params = {
        "format": "netcdf",
        "product_type": "reanalysis",
        "variable": [
            "surface_net_solar_radiation",
            "surface_net_solar_radiation_clear_sky",
            "surface_net_thermal_radiation",
            "surface_net_thermal_radiation_clear_sky",
            "surface_solar_radiation_downward_clear_sky",
            "surface_solar_radiation_downwards",
            "surface_thermal_radiation_downward_clear_sky",
            "surface_thermal_radiation_downwards",
        ],
        "year": [yearstr],
        "month": [monthstr],
        "day": [daystr],
        "time": times,
        "grid": [0.25, 0.25],
        "area": [north, west, south, east],
    }
    # retrieves the path to the file
    fl = c.retrieve(dataset, params)
    # download the file
    if download_flag:
        fl.download(outdir + outfile)

# TODO: delete following line
dataset = "reanalysis-era5-pressure-levels"
download_flag = False
# Get met variables for pressure levels
c = cdsapi.Client()
for day in days[1:]:
    daystr = str(day).zfill(2)
    outfile = "era5_esc_202202" + daystr + ".nc"

    # api parameters:
    # area: [north, west, south, east]; west is negative
    params = {
        "format": "netcdf",
        "product_type": "reanalysis",
        "variable": [
            "pressure",
            "Geopotential",
            "Temperature",
            "Relative humidity",
            "Ozone mass mixing ratio",
        ],
        "pressure_level": press_levs,
        "year": [yearstr],
        "month": [monthstr],
        "day": [daystr],
        "time": times,
        "grid": [0.25, 0.25],
        "area": [north, west, south, east],
    }
    # retrieves the path to the file
    fl = c.retrieve(dataset, params)
    # download the file
    if download_flag:
        fl.download(outdir + outfile)
