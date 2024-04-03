#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:14:02 2023

@author: prowe

Get the profile data (P, T, H2O, ozone) from ERA5 for the
entire Escudero period.
"""

import numpy as np

from antarc.escudero.parameters import esc_params
from antarc.get_era5_from_web import get_paramsets, get_era5
import calendar

# Parameter modules, specific to user
from antarc import params


LATITUDE = esc_params.LATITUDE
LONGITUDE = esc_params.LONGITUDE
ALTITUDE = esc_params.ALTITUDE
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

# 2017: 01: Done
# 2017: 11, 12: Done and transferred;
# 2018: 01, 02, 03: Done and transferred
# 2018: 11, 12: Done and transferred
# 2019: 01, 02: 03: Done and transferred
# 2019: 09, 10, 11, 12: Done and transferred
# 2020: 01, 02, 03: Done and transferred
# 2020: 12: Done and transferred
# 2021: 01, 02: Done and transferred
# 2021: 12: Done and transferred
# 2022: 1, 2, 3, 4, 5, 6, 7, 8: Done and transferred
# 2022: 11, 12: Done and transferred
# 2023: 1, 2, 3, 4, 5, 6, 7: Done and transferred, in progress: 8
years = [2023]
for year in years:
    for month in range(8, 9):

        yearstr = str(year)
        monthstr = str(month).zfill(2)
        num_days = calendar.monthrange(year, month)[1]
        days = list(range(1, num_days + 1))

        prefix = "era5_esc_" + yearstr + monthstr

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
