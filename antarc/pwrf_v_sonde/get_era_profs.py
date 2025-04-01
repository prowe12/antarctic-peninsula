#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:25 2022

@author: prowe
"""

import datetime as dt


# Modules
from antarc.pwrf_v_sonde.get_atm_profs_rsrc import FileInfo, Era
from antarc.get_era5_press_levs_from_web import get_era5_press_levs


# # # # # #   ERA5   # # # # #
def get_era(prefix, params, atm_prof_params, date_str, date_fmt):
    lat = params.LATITUDE
    lon = params.LONGITUDE
    dtime = dt.datetime.strptime(date_str, date_fmt)

    direc = f"{atm_prof_params.ERA_DIR}{dtime.year}/"
    era_fileformat = atm_prof_params.ERA_FILEFORMAT

    # Get the ERA5 files form the web if needed
    date1 = dt.datetime(dtime.year, dtime.month, dtime.day)
    get_era5_press_levs(lat, lon, date1, date1, direc, prefix, False)

    Files = FileInfo(direc, era_fileformat)

    # Get the ERA5 data
    era = Era(Files, era_fileformat, dtime, lat, lon)
    return era
