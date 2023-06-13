#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:14:02 2023

@author: prowe
"""

import numpy as np
from os.path import exists
import math

from antarc.get_era5_from_web import get_paramsets, get_era5


def get_era5_press_levs(lat, lon, date1, date2, outdir, prefix, redo):

    # MEAS_DIR = params.MEAS_DIR
    # outdir = MEAS_DIR + "Escudero/era5/"
    # prefix = "era5_esc_" (+ date1.strftime("%Y%m"))

    north = math.ceil(lat * 4) / 4
    east = math.ceil(lon * 4) / 4
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
    times = [f"0{i}:00" for i in range(10)] + [
        f"{i}:00" for i in range(10, 24)
    ]
    dataset = "reanalysis-era5-pressure-levels"
    # download_flag = True

    # Get date; pad out four days on each side
    yearstr = str(date1.year)
    monthstr = str(date1.month).zfill(2)
    days = list(
        range(date1.day, date2.day + 1)
    )  # list(range(date1.day - 4, date2.day + 4))
    prefix += date1.strftime("%Y%m")

    paramsets = get_paramsets(
        prefix,
        yearstr,
        monthstr,
        days,
        times,
        latlon,
        press_levs,
    )

    # If redo is false, purge the days that were already downloaded
    if not redo:
        for i in range(len(paramsets) - 1, -1, -1):
            if exists(outdir + paramsets[i]["outfile"]):
                paramsets.pop(i)

    if len(paramsets) > 0:
        get_era5(outdir, paramsets, dataset)
