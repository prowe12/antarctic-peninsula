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


# Load in LWD and SWD clear and cloud from ERA5
meas_dir = f"{params.MEAS_DIR}Escudero"
era5rad_dir = f"{meas_dir}/era5/broadband/{date1.strftime('%Y')}/"
fmt = "era5_esc_broadband_*[0-9]*.nc"
era5 = read_era5_broadband_down_by_file(era5rad_dir, fmt, date1, date2)

# Closest grid point
i = np.where(era5["lat"] == -62.25)[0][0]
j = np.where(era5["lon"] == -59.0)[0][0]
era_swd = era5["swd"][:, i, j]
era_swdc = era5["swd_clr"][:, i, j]
era_lwd = era5["lwd"][:, i, j]
era_lwdc = era5["lwd_clr"][:, i, j]

# Closest grid point


launch_time = dt.datetime(2017, 1, 31, 18, 20, 0, tzinfo=pytz.utc)
lwd_clear_date[1] = launch_time
ilaunch = np.where(swd_date >= launch_time)[0][0]

lim1 = [dt.datetime(2017, 1, 30, 23, 45, 0), dt.datetime(2017, 2, 1), -10, 910]
lim2 = [dt.datetime(2017, 1, 30, 23, 45, 0), dt.datetime(2017, 2, 1), 200, 340]

# Plot the results
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 6), sharex=True, num=1)
wid = 0.85
hght = 0.42
lft = 0.1
ax1.set_position([lft, 0.52, wid, hght])
ax2.set_position([lft, 0.06, wid, hght])

ax1.plot([launch_time, launch_time], [-20, 1000], "k:")
ax1.plot([lwd_clear_date[0], lwd_clear_date[0]], [-20, 1000], "k:")
ax1.plot(era5["dates"], era_swdc, "b--", label="ERA5 clear")
ax1.plot(swd_date, swd, color="orange", label="measured")
ax1.plot(era5["dates"], era_swd, "o", label="ERA5")
ax1.set_ylabel("SWD (W/m$^{2}$)")
ax1.axis(lim1)
ax1.legend()

ax2.plot(era5["dates"], era_lwd, "o", label="ERA5")
ax2.plot(era5["dates"], era_lwdc, "b--", label="ERA5 clear")
ax2.plot([launch_time, launch_time], [195, 345], "k:")
ax2.plot([lwd_clear_date[0], lwd_clear_date[0]], [195, 345], "k:")
ax2.plot(lwd_clear_date, lwd_clear, "s", label="Clear (sonde)")
ax2.axis(lim2)
ax2.set_ylabel("LWD (W/m$^{2}$)")
ax2.legend(
    bbox_to_anchor=(0.65, 1.0),
    loc="upper center",
    borderaxespad=0,
)
# frameon=False,

outfile = "downwelling_radiation.png"
plt.savefig(outfile, format="png")


era_date = [
    dt.datetime.strptime(str(x), "%Y-%m-%d %H:%M:%S") for x in era5["dates"]
]
era_ts = np.array([to_timestamp(x) for x in era_date])
swd_ts = np.array([to_timestamp(x) for x in swd_date])
era_swdc_interp = np.interp(swd_ts, era_ts, era_swdc)

lwd_clear_ts = np.array([to_timestamp(x) for x in lwd_clear_date])
era_lwd_interp = np.interp(lwd_clear_ts, era_ts, era_lwd)

# ax1.plot(swd_date, era_swdc_interp, "+")


# Plot the cloud forcing

# def convert_to_dt(x):
#     return dt.strptime(str(x), '%Y-%m-%d %H:%M:%S')

lim1[2:4] = [-500, 200]
lim2[2:4] = [0, 90]

fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 6), sharex=True, num=2)
ax1.plot([lwd_clear_date[0], lwd_clear_date[0]], [-500, 200], "k:")
ax1.plot([lwd_clear_date[1], lwd_clear_date[1]], [-500, 200], "k:")
ax1.plot(era5["dates"], era_swd - era_swdc, "o", label="ERA5 - ERA5 clear")
ax1.plot(swd_date, swd - era_swdc_interp, label="Measured - ERA5 clear")
ax1.axis(lim1)
ax1.legend()

ax1.legend(
    loc="lower left",
)

ax1.set_ylabel("SWD Cloud Forcing (W/m$^{2}$)")
ax2.plot([lwd_clear_date[0], lwd_clear_date[0]], [0, 100], "k:")
ax2.plot([lwd_clear_date[1], lwd_clear_date[1]], [0, 100], "k:")
ax2.plot(era5["dates"], era_lwd - era_lwdc, "o", label="ERA5 - ERA5 clear")
ax2.plot(
    lwd_clear_date,
    era_lwd_interp - lwd_clear,
    "s",
    label="ERA5 - Clear (sonde)",
)
ax2.set_ylabel("LWD Cloud Forcing (W/m$^{2}$)")
ax2.legend()
ax2.axis(lim2)

outfile = "cloud_forcing.png"
plt.savefig(outfile, format="png")

# # plt.figure()
# # plt.subplot(212)
# # plt.plot(lwd_date, lwd, label="measured")
# # plt.plot(lwd_clear_date, lwd_clear, "s", label="clear")
# # plt.plot(lwd_clear_date, lwd_force, "+", label="forcing")
# # plt.ylabel("Downward longwave (W/m$^{2}$)")
# # plt.legend()

# # plt.subplot(211)
# # plt.plot(swd_date, swd, label="measured")
# # plt.plot(swd_clear_date, swd_clear, "s", label="clear")
# # plt.plot(swd_clear_date, swd_force, "+", label="forcing")
# # plt.ylabel("Downward shortwave (W/m$^{2}$)")
# # plt.legend()


# # # Cumulative forcing for SWD
# # swd_tstamp = [x.timestamp() for x in swd_date]
# # swd_clear_lib_interp = np.interp(swd_tstamp, lib_tstamp, swd_clear_lib)
# # i1 = np.where(np.array(swd_tstamp) >= start_date.timestamp())[0][0]
# # i2 = np.where(np.array(swd_tstamp) <= final_date.timestamp())[0][-1]
# # swd_force = swd - swd_clear_lib_interp
# # ntimes = i2 - i1
# # swf_tot = np.zeros([ntimes])
# # swf_tot_mwh = np.zeros([ntimes])
# # count = 0
# # for i in range(i1, i2):
# #     swf = swd_force[i] * (swd_tstamp[i] - swd_tstamp[i - 1])
# #     swf_mwh = swd_force[i] * (swd_tstamp[i] - swd_tstamp[i - 1]) / 3600 / 1000
# #     swf_tot[count] = swf_tot[count - 1] + swf
# #     swf_tot_mwh[count] = swf_tot_mwh[count - 1] + swf_mwh
# #     count += 1
# # is1 = i1
# # is2 = i2


# # # lwd_tot = np.interp(lwd_tstamp, lwd_clear_interp)
# # # swd_tot = np.interp(lwd_tstamp, lwd_clear_interp)
# # lwf_tot_sw_date = np.interp(swd_tstamp[is1:is2], lwd_tstamp[il1:il2], lwf_tot)
# # lwd_force_sw_date = np.interp(swd_tstamp, lwd_tstamp, lwd_force)

# # # ax1 = plt.subplot(221)
# # # ax2 = plt.subplot(222, sharex=ax1)
# # # ax3 = plt.subplot(223, sharex=ax1)
# # # ax4 = plt.subplot(224, sharex=ax1)

# # # ax2.plot(lwd_date[:ilw], lwd_force[:ilw], color="blue", label="LWD Forcing")
# # # ax2.plot(swd_date[:isw], swd_force[:isw], color="red", label="SWD Forcing")
# # # ax2.legend()
# # # ax2.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# # # ax2.tick_params(axis="x", rotation=rot)

# # # ax4.plot(lwd_date[:ilw], lwf_tot[:ilw] / 1e6, color="blue", label="LWD")
# # # ax4.plot(swd_date[:isw], swf_tot[:isw] / 1e6, color="red", label="SWD")
# # # ax4.plot(
# # #     swd_date[:isw],
# # #     (swf_tot[:isw] + lwf_tot_sw_date[:isw]) / 1e6,
# # #     "k",
# # #     label="Total",
# # # )
# # # ax4.legend()
# # # ax4.set_ylabel("Cumulative Forcing (MJ/m$^{2}$)")
# # # ax4.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# # # ax4.set_xlabel("Day of 2022/02")
# # # ax4.tick_params(axis="x", rotation=rot)

# # totf = swf_tot + lwf_tot_sw_date
# # totfstr = str(round(totf[-1] / 1e6))

# # plot_measured_and_clear()
# # # plot_meas_clear_forcings_cumulative_shorttime()
# # plot_meas_clear_forcings()

# # if plotfigs:
# #     # .. Plot it
# #     plt.figure(1)
# #     plt.clf()
# #     plt.plot(dtimes, rad, label="pyranometer")
# #     plt.plot(dtimes, diffuse, label="pysolar")
# #     plt.plot(dtimes, diffuse * diffuse_fac, label="clear")

# #     # .. Cloud forcing
# #     sw_forcing = diffuse * diffuse_fac - rad
# #     plt.plot(dtimes, sw_forcing, label="forcing")
# #     plt.legend()
