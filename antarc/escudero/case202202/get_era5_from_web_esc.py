#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:14:02 2023

@author: prowe
"""

import numpy as np

from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import esc202202
from antarc.get_era5_from_web import get_paramsets, get_era5

# Parameter modules, specific to user
from antarc import params


LATITUDE = esc_params.LATITUDE
LONGITUDE = esc_params.LONGITUDE
ALTITUDE = esc_params.ALTITUDE

DATE1 = esc202202.DATE1
DATE2 = esc202202.DATE2


MEAS_DIR = params.MEAS_DIR
outdir = MEAS_DIR + "Escudero/era5/"

north = np.round(LATITUDE / 0.25) * 0.25 + 0.25
east = np.round(LONGITUDE / 0.25) * 0.25 + 0.25
latlon = {
    "north": north,
    "south": north - 0.25,
    "east": east,
    "west": east - 0.5,
}

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
dataset = "reanalysis-era5-pressure-levels"
download_flag = True

yearstr = str(DATE1.year)
monthstr = str(DATE1.month).zfill(2)
days = list(range(DATE1.day - 1, DATE2.day + 2))  # pad out a day on each side

prefix = "era5_esc_202202"

paramsets = get_paramsets(
    prefix,
    yearstr,
    monthstr,
    days,
    times,
    latlon,
    press_levs,
)
get_era5(outdir, paramsets, dataset)
