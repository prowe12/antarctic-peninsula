#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 14:57:55 2023

@author: prowe

Purpose: Get cloud forcing at Escudero for 2017/01/31

"""

# Dependencies
import datetime as dt
import numpy as np
import pytz
import matplotlib.pyplot as plt
from os.path import exists
import pandas as pd
import calendar


# My modules
from antarc.run_radtran import get_clear_fluxes
from antarc.load_rad_flux import load_rad_flux
from antarc.escudero.get_lwd_flux import get_lwd_clear_date
from antarc.escudero.read_era5_rad import read_era5_broadband_down_by_file

from antarc.escudero.get_swd_flux import get_obs_swd, do_run_libradtran
from antarc.escudero.get_era5_from_web import get_era5_from_web
from antarc.escudero.all.get_met_from_web import get_met_15_minutes
from antarc.escudero.all.get_era5_rad_from_web import get_era5_rad_from_web

from antarc.escudero.case201701.get_atm_profs import get_atm_prof


# Parameter modules specific to the user: the user must create params.py
# in the directory antarc with the variables MEAS_DIR and PROJ_DIR set
from antarc import params

# Parameter modules
from antarc.escudero.parameters import esc202205 as esc_case
from antarc.escudero.parameters import esc_params, atm_prof_params
from antarc.escudero.parameters import swd_params, lwd_params, radtran_params


def to_timestamp(d):
    return calendar.timegm(d.timetuple())


MEAS_DIR = params.MEAS_DIR
dir_lwdclear = MEAS_DIR + "Escudero/lwd_clear/"
out_dir = params.PROJECT_DIR + "case_studies/model_performance_KGI_2017_01/"
prof_dir = out_dir + "ground_measurements/Escudero/profiles/"

sonde_file = "esc_sonde_dd_v2_2017013118.txt"

# Run flags
SAVEFIGS = False
REDO = True  # redo calculations even if they already exist


# Params needed from param files
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
# alt_surf = esc_params.ALTITUDE
# location = esc_params.LOCATION


# A5) Get ERA5 data
# get_era5_from_web(esc_params, esc_case, REDO)

date1 = dt.datetime(2017, 1, 31, 12, 0, tzinfo=pytz.utc)
date2 = dt.datetime(2017, 1, 31, 23, 59, tzinfo=pytz.utc)


# # Load in clear sky calculation from profile
# # (No pyrgeometer for this time period)
# lwd_clear_date, lwd_clear = get_lwd_clear_date(
#     lwd_params.LWD_CLEAR_DIR,
#     lwd_params.LWD_CLEAR_FILEFORMAT_LOWNU,
#     lwd_params.LWD_CLEAR_FILEFORMAT_HINU,
#     lwd_params.LWD_CLEAR_FILEFORMAT_OUTNU,
#     date1,
#     date2,
# )

# C) Longwave
# C2) Create the atmospheric profiles
get_atm_prof(
    atm_prof_params, lat, lon, date1, date2, sonde_file, prof_dir, REDO
)


# lwd_tstamp = [x.timestamp() for x in lwd_date]
# lwd_clear_tstamp = [x.timestamp() for x in lwd_clear_date]
# lwd_clear = np.array(lwd_clear)
# # lwd_clear_interp = np.interp(lwd_tstamp, lwd_clear_tstamp, lwd_clear)
# lwd_interp = np.interp(lwd_clear_tstamp, lwd_tstamp, lwd)

# # C5) Subtract cloudy and clear-sky fluxes to get LW forcing
# # lwd_force = lwd - lwd_clear_interp
# lwd_force = lwd_interp - lwd_clear

# # D) SWD
# # Setup
# year_str = date1.strftime("%Y")
# if year_str != date2.strftime("%Y"):
#     raise ValueError("Different years must be processed separately")

# D1) Load pyranometer data
swd_dir = f"{swd_params.SWD_DIR}{date1.strftime('%Y')}/"
swd_fmt = swd_params.SWD_FILEFORMAT
swd_date0, swd = load_rad_flux(swd_dir, swd_fmt, date1, date2)
swd_date = np.array(swd_date0)
# swd_date, swd = get_obs_swd(swd_params, esc_case)
swd[swd < 0] = 0
# # swd_date_long, swd_long = load_rad_flux(swd_dir, swd_fmt, date1, date2)


# # D2) If needed, calculate the libradtran results and save to file
# # Option 1: Use libradtran to calculate
# # # Spacing: every 15 minutes
# # timestep_s = range(0, 86400, 900)

# # REDO = False
# # for iday in range(366):
# #     this_date = date1 + dt.timedelta(days=iday)
# #     if this_date > date2:
# #         print(this_date)
# #         break

# #     swd_clear_date = [this_date + dt.timedelta(seconds=x) for x in timestep_s]

# #     datestr = f"{this_date.strftime('%Y%m%d')}"
# #     clrfile = f"{swd_params.SWD_CLEAR_DIR}{year_str}/libradtran_{datestr}.csv"

# #     if not exists(clrfile) or REDO:
# #         libclear = do_run_libradtran(
# #             esc_params.LATITUDE, esc_params.LONGITUDE, swd_clear_date, clrfile
# #         )
# #         # do_run_libradtran(lat, lon, dtimes, outfile)
# #         swd_clear = libclear[:, 0]
# #     else:
# #         # Load the libradtran results, which have columns:
# #         libclear = pd.read_csv(clrfile, delimiter=",")
# #         swd_clear_date = pd.to_datetime(libclear["Date"])
# #         swd_clear = libclear["global"]

# #     # libdate, swd_clear = get_libradtran(esc_params, swd_date, clrfile, REDO)
# #


# # # uncomment following
# # # lib_tstamp = [x.timestamp() for x in libdate]
# # swd_clear_tstamp = [x.timestamp() for x in swd_clear_date]
# # swd_tstamp = [x.timestamp() for x in swd_date]
# # # swd_clear = np.interp(swd_tstamp, lib_tstamp, swd_clear)
# # swd_interp = np.interp(swd_clear_tstamp, swd_tstamp, swd)

# # # D3) Subtract cloudy and clear-sky fluxes to get SW forcing
# # # swd_force = swd - swd_clear
# # swd_force = swd_interp - swd_clear
