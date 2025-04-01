#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 08:56:46 2022

@author: prowe

Purpose: Plot broadband radiation at Escudero on 2022/02


Steps for Analysis:

The user must create a file called params.py within the git repo
antarctic_peninsula/antarc/. This file is ignored by git, since it is
unique to each compuater. In params.py, set the measurement directory 
for the computer, e.g. MEAS_DIR = "/my-user-path/measurements"
In the description that follows, params.MEAS_DIR is referred to as measdir

Steps indicated with asterisks are performed in this code. Other steps
must be performed manually or using other code before running this code, 
as indicated.

A) Collect the data and put it in the directories within 'measdir/Escudero/'
  1)  LWD (Longwave downward radiation): pyranometer/orig_v1. (Do manually)
  2)  SWD (Shortwave downward radiation): pyrgeometer/orig_v2 (Do manually)
  3)  Sonde: GRAW_radiosondes/simulation/yopp2022_escudero_profiles (Manually)
  4*) ERA5: era5/; Downloaded using /antarctic_peninsula/antarc/escudero/
      get_era5_from_web.py (done here; see code below)
  5)  CO2: co2/; See readme_carbontracker.rtf, readme_co2_surface_flask.rtf,
      and readme_co2_noaa_esrl.rtf (Manually)

B) Convert raw data files into standardized data files, using code in 
   antarctic_peninsula/antarc/escudero/get_stand_files:
  1)  LWD: pyr_save_stand_files_lwd.py, uses antarc.pyr_save_stand_files and
      antarc/escudero/parameters/lwd_orig_params.py. Results are saved in 
      measdir/Escudero/pyranometer/stand/yyyy/. (Do before running this code)
  2)  SWD: pyr_save_stand_files_swd.py, uses antarc.pyr_save_stand_files and
      antarc/escudero/parameters/swd_orig_params. Results are saved in in
      measdir/Escudero/pyrgeometer/stand/yyyy/ (Do before running this code)
  3*) Sonde: sonde_raw_to_data_denial_format.py. Results are saved in
      measdir/measurements/Escudero/radiosondes/datadenial. (done here)
  4)  ERA5: These files are left in their original format. (Nothing to do)
  5)  CO2: These files are left in their original format. (Nothing to do)

C) Longwave:
  1)  Requires steps A1, A3, A4, A5, B1, and B3 above.
  2*) get_atm_profs.py: Create atmospheric profiles from radiosondes and 
      ancillary data (done here; see below)
  3*) Do LBLRTM runs to get clear-sky fluxes. These go in to lwd_clear directory.
      (Done here using get_clear_fluxes.)
  4*) Load in the lwd measurements and the clear sky files (done here)
  5*) Subtract cloudy and clear-sky fluxes to get LW forcing (done here)

D) Shortwave:
  1*) Load pyranometer data using get_obs_swd (done here)
  2*) Compute clear-sky SWD flux using get_libradtran (done here)
  3*) Subtract cloudy and clear-sky fluxes to get SW forcing (done here)

E) miniMPL: Please see the repos mpl2nc and mpl-processing.

"""

# Dependencies
import datetime as dt
import numpy as np
import pytz
import matplotlib.pyplot as plt
from os.path import exists
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from netCDF4 import Dataset

# from scipy.interpolate import CubicSpline

# My modules
from antarc.run_radtran import get_clear_fluxes
from antarc.escudero.get_lwd_flux import get_lwd, get_lwd_clear
from antarc.escudero.get_swd_flux import get_obs_swd, get_libradtran
from antarc.escudero.get_era5_from_web import get_era5_from_web
from antarc.escudero.case202205.get_atm_profs import get_atm_profs
from antarc.escudero.case202205.make_plots import plot_measured_and_clear
from antarc.escudero.case202205.make_plots import plot_forcings
from antarc.escudero.case202205.make_plots import plot_meas_clear_forcings
from antarc.escudero.case202205.make_plots import plot_forcings_and_lidar
from antarc.escudero.case202205.make_plots import get_cmap


# Parameter modules specific to the user: the user must create params.py
# in the directory antarc with the variables MEAS_DIR and PROJ_DIR set
from antarc import params

# Parameter modules
from antarc.escudero.parameters import esc202205 as esc_case
from antarc.escudero.parameters import esc_params, atm_prof_params
from antarc.escudero.parameters import swd_params, lwd_params, radtran_params
from antarc.escudero.all.read_era5_bbnd import read_era5_broadband

utc = dt.timezone.utc


def get_lidar_data(lidar_file):
    with Dataset(lidar_file) as ncid:
        # dict_keys(['DataTime', 'Range', 'CloudMask', 'PhaseMask',
        #   'CloudBaseAltitude', 'CloudBaseTemperature', 'CloudTopAltitude',
        #   'CloudTopTemperature', 'ColumnType', 'Temperature', 'CountRate',
        #   'Depolarization', 'BackscatterRatio'])

        #% This colortable makes the following values in the data mask a single color
        #% Take the log10 of the mask, you will get data that is colored as:
        #% Good data = green, clipped data = red, missing = blue, SNR filtered = white
        #% and density filtered = yellow. If multiple filters are tripped, the
        #% highest number is plotted
        colortable = [
            [0, 1, 1],  # Clear
            [0.75, 0.75, 0.75],  # Cloud
            [0, 0, 1],  # Liquid
            [1, 0, 0],  # Ice
        ]

        alt = ncid["Range"][:] / 1000
        timestamps = ncid["DataTime"][:]

        # Create a combined mask
        combo = -999 * np.ones(ncid["PhaseMask"].shape)
        for i in range(ncid["PhaseMask"].shape[1]):
            combined = ncid["PhaseMask"][:, i] + 2
            # TODO: what is the following for?
            # if not np.all(ncid["CloudMask"][np.isnan(combined), i]):
            #     raise ValueError("Found cloud where NaN expected")
            inan = np.where(np.isnan(combined))[0]
            combined[inan] = ncid["CloudMask"][inan, i] + 0.6
            combo[:, i] = combined

        val = [combo]

    # cbarbnds = np.array([[0, 10]])
    cmaps = [[np.arange(0.5, 5, 1), colortable]]
    cticks = cmaps[0][0][:-1] + 0.5
    cticklabels = ["Clear", "Cloud", "Liquid", "Ice"]

    # Discrete pixels
    i = 0
    bounds = cmaps[i][0]
    mymap = get_cmap(cmaps[i][1])
    norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=mymap.N)

    return timestamps, alt, val[i].T, cticks, cticklabels, norm, mymap


# TODO: get this info from the parameter file
MEAS_DIR = f"{params.MEAS_DIR}Escudero/"
dir_lwdclear = MEAS_DIR + "lwd_clear/"
out_dir = params.PROJECT_DIR + "case_studies/May_2022/figures/"
era_broadband_dir = f"{MEAS_DIR}era5/broadband/"
# lidar_dir = f"{MEAS_DIR}mpl/reprocessing/ResultNetCDF/with_era5_profiles/"
lidar_dir = f"{MEAS_DIR}mpl/reprocessing_stillwell/v3/"


# Run flags
SAVEFIGS = False
savefig = False
REDO = False  # redo calculations even if they already exist
# redo_profs = False

esc_case.DATE1 = esc_case.DATE1.replace(tzinfo=None)
esc_case.DATE2 = dt.datetime(2022, 5, 18)  # , tzinfo=pytz.UTC)
date1 = esc_case.DATE1
date2 = esc_case.DATE2


# Already done
#
# # A5) Get ERA5 data if it is not found or if REDO=TRUE
# get_era5_from_web(esc_params, esc_case, REDO)
#
# # C) Longwave
# # C2) Create the atmospheric profiles
# get_atm_profs(atm_prof_params, esc_params, esc_case, REDO)
#

# C3)
# Check for the clear-sky lwd fluxes and if the files have not been made,
# create them. We expect one per profile
get_clear_fluxes(dir_lwdclear, radtran_params, (date1, date2), 3, REDO)

# C4) Load in measured LWD and clear sky files created in previous step
lwd_date, lwd = get_lwd(lwd_params, date1, date2)
lwd_clear_date, lwd_clear = get_lwd_clear(esc_case, lwd_params)
lwd_tstamp = [x.timestamp() for x in lwd_date]
lwd_clear_tstamp = [x.timestamp() for x in lwd_clear_date]
lwd_clear = np.array(lwd_clear)
lwd_clear_interp = np.interp(lwd_tstamp, lwd_clear_tstamp, lwd_clear)


# D) SWD
# Setup
year1 = date1.strftime("%Y")
datestr = f"{date1.strftime('%Y%m%d')}_{date2.strftime('%Y%m%d')}"
if year1 != date2.strftime("%Y"):
    raise ValueError("Different years must be processed separately")
swd_dir = f"{swd_params.SWD_DIR}{year1}/"
swd_fmt = swd_params.SWD_FILEFORMAT
clrfile = f"{swd_params.SWD_CLEAR_DIR}libradtran_{datestr}.csv"

# D1) Load pyranometer data
swd_date, swd = get_obs_swd(swd_params, esc_case)
swd[swd < 0] = 0
swd_tstamp = [x.timestamp() for x in swd_date]
# swd_date_long, swd_long = load_rad_flux(swd_dir, swd_fmt, date1, date2)

# D2)
# Option a: use libradtran results
# swd_date, swd, swd_clear_lib, swd_clear, diffuse, diffuse_fac = get_swd()

# Spacing: every 15 minutes
timestep_s = range(0, 86400, 900)
year_str = date1.strftime("%Y")

REDO = False
swd_clear = []
libdate = []
this_date = date1
iday = 0
while this_date <= date2:
    this_date = date1 + dt.timedelta(days=iday)
    iday += 1
    #     if this_date > date2:
    #         print(this_date)
    #         break

    datestr = f"{this_date.strftime('%Y%m%d')}"
    clrfile = f"{swd_params.SWD_CLEAR_DIR}{year_str}/libradtran_{datestr}.csv"

    if not exists(clrfile) or REDO:
        raise ValueError("modify code to create this file")
        # run_date = [this_date + dt.timedelta(seconds=x) for x in timestep_s]
        # libclear = do_run_libradtran(
        #     esc_params.LATITUDE,
        #     esc_params.LONGITUDE,
        #     run_date,
        #     clrfile,
        # )
        # # do_run_libradtran(lat, lon, dtimes, outfile)
        # swd_clear = libclear[:, 0]
    else:
        # Load the libradtran results, which have columns:
        libclear = pd.read_csv(clrfile, delimiter=",")
        swd_clear_date0 = pd.to_datetime(libclear["Date"]).to_list()
        swd_clear0 = libclear["global"].to_list()
    libdate += swd_clear_date0
    swd_clear += swd_clear0


# libdate, swd_clear = get_libradtran(esc_params, swd_date, clrfile, REDO)
lib_tstamp = [x.timestamp() for x in libdate]
swd_tstamp = [x.timestamp() for x in swd_date]
swd_clear_interp = np.interp(swd_tstamp, lib_tstamp, swd_clear)

# Option b: Get from ERA5
era5rad_outdir = f"{era_broadband_dir}/{date1.strftime('%Y')}/"
# get_era5_rad_from_web(era5rad_outdir, date1, date2, lat, lon)
fmt = "era5_esc_broadband_*[0-9]*.nc"
era5bb = read_era5_broadband(era5rad_outdir, fmt, date1, date2)

era5_date_list = [x for x in era5bb["date"]]
# [x.replace(tzinfo=pytz.utc) for x in era5bb["date"]]
era5_tstamp = [x.timestamp() for x in era5_date_list]
era5_date = np.array(era5_date_list)

# Get closest grid point
i = np.where(era5bb["lat"] == -62.25)[0][0]
j = np.where(era5bb["lon"] == -59.0)[0][0]

era5_swd = era5bb["swd"][:, i, j]
era5_swd_clr = era5bb["swd_clr"][:, i, j].data
era5_swd_clr_intp = np.interp(swd_tstamp, era5_tstamp, era5_swd_clr)

# QC the interpolation
plt.figure()
plt.plot(era5_tstamp, era5_swd_clr)
plt.plot(swd_tstamp, era5_swd_clr_intp)

era5_lwd = era5bb["lwd"][:, i, j].data
era5_lwd_clr = era5bb["lwd_clr"][:, i, j].data

era5_swnet = era5bb["swnet"][:, i, j]
era5_swnet_clr = era5bb["swnet_clr"][:, i, j].data
era5_lwnet = era5bb["lwnet"][:, i, j].data
era5_lwnet_clr = era5bb["lwnet_clr"][:, i, j].data

# Demonstrate that ERA5 needs to be shifted back by 1/2 hour
halfhr = dt.timedelta(hours=0.5)
plt.figure()
plt.subplot(211)
plt.plot(era5_date, era5_swd_clr, ".", label="ERA5 SWD clear")
plt.plot(swd_date, era5_swd_clr_intp, label="ERA5 SWD clear interp")
plt.plot(libdate, swd_clear, label="SWD clear libradtran")
plt.plot(swd_date, swd, label="SWD")
plt.legend()
plt.subplot(212)
plt.plot(swd_date - halfhr, era5_swd_clr_intp, label="ERA5 SWD clear interp")
plt.plot(era5_date - halfhr, era5_swd_clr, ".", label="ERA5 SWD clear offset")
plt.plot(libdate, swd_clear, label="SWD clear libradtran")
plt.plot(swd_date, swd, label="SWD")
plt.legend()


# Shift era5 dates back by 1/2 hour
era5_date = np.array(era5_date_list) - halfhr
era5_tstamp = [x.timestamp() for x in era5_date]

# For the smoothed SWD clear, interpolate the libradtran results.
# (Do not use spline because it is poorly behaved at the edges
# and sources of uncertainty are high anyway.) Then scale to
# match ERA5 for this day
# swd_clear_interp = CubicSpline(lib_tstamp, swd_clear)(swd_tstamp)
swd_clear_interp = np.interp(swd_tstamp, lib_tstamp, swd_clear)

# splinefun = CubicSpline(era5_tstamp, era5_swd_clr)
# era5_swd_clr_intp = splinefun(swd_tstamp)
era5_swd_clr_intp = np.interp(swd_tstamp, era5_tstamp, era5_swd_clr)
# era5_swd_clr_intp = swd_clear_interp * 0.98


# D3) Subtract cloudy and clear-sky fluxes to get SW forcing
swd_force = swd - era5_swd_clr_intp
era5_swd_force = era5_swd - era5_swd_clr

# C5) Subtract cloudy and clear-sky fluxes to get LW forcing
era5_lwd_clr_intp = np.interp(lwd_tstamp, era5_tstamp, era5_lwd_clr)
lwd_force = lwd - era5_lwd_clr_intp
era5_lwd_force = era5_lwd - era5_lwd_clr


# Plot the fluxes
start_date = dt.datetime(2022, 5, 14)
final_date = dt.datetime(2022, 5, 18)

plt.figure()
plt.subplot(211)
plt.plot(lwd_date, lwd, label="LWD")
plt.plot(lwd_date, era5_lwd_clr_intp, label="LWD clear ERA5 interp")
plt.plot(era5_date, era5_lwd_clr, "+", label="LWD clear ERA5")
plt.plot(lwd_date, lwd_clear_interp, label="LWD clear interp")
plt.plot(lwd_clear_date, lwd_clear, "o", label="LWD clear")
plt.legend()
plt.xlim([start_date, final_date])

plt.subplot(212)
plt.plot(era5_date, era5_swd_clr, ".", label="SWD clear ERA5")
plt.plot(swd_date, swd_clear_interp, label="SWD clear interp")
plt.plot(swd_date, era5_swd_clr_intp, label="ERA5 SWD clear interp")
plt.plot(libdate, swd_clear, label="SWD clear libradtran ")
plt.plot(swd_date, swd, label="SWD measured")
plt.legend()
plt.xlim([start_date, final_date])


# Figure of downward fluxes
fontsize = 12
start_date = dt.datetime(2022, 5, 15, 20)
final_date = dt.datetime(2022, 5, 17, 12)


fig, axs = plt.subplots(2, 1, constrained_layout=True, num=8, clear=True)

ax = axs[0]
ax.set_position((0.11, 0.6, 0.83, 0.35))
ax.plot(lwd_date, lwd, label="Sky, Meas.")
ax.plot(era5_date, era5_lwd, ".", label="Sky, ERA5")
ax.plot(era5_date, era5_lwd_clr, "g.", label="Clear, ERA5")
ax.plot(lwd_clear_date, lwd_clear, "cs", label="Clear, sonde")
ax.set_ylabel("LWD flux (W/m$^2$)")
ax.legend(loc="upper right", bbox_to_anchor=[1.07, 1.1])
# " left", bbox_to_anchor=[.15,.3])
ax.set_xlim([start_date, final_date])
ax.set_ylim([200, 350])
ax.tick_params(axis="x", rotation=20)

ax = axs[1]
ax.set_position((0.11, 0.1, 0.83, 0.35))
ax.plot(libdate, swd_clear, "c", label="Clear, Simulated")
ax.plot(era5_date, era5_swd_clr, "g.", label="Clear, ERA5")
ax.plot(swd_date, swd, label="Sky, Measured")
ax.plot(era5_date, era5_swd, ".", label="Sky, ERA5")
ax.set_ylabel("SWD flux (W/m$^2$)")
ax.legend()
ax.set_xlim([start_date, final_date])
ax.set_ylim([-5, 120])
ax.tick_params(axis="x", rotation=20)
ax.set_xlabel("Hour on 2022/05/16", fontsize=fontsize)

if savefig:
    plt.savefig("flux.png", format="png")
    plt.savefig("flux.eps", format="eps")


# Save the data for Dave
import pandas as pd

data_dir = "antarc/escudero/case202205/data/"

# DLW
fname = data_dir + "meas_dlw.csv"
data = pd.DataFrame(
    np.hstack((lwd_date[:, None], lwd[:, None])),
    columns=["Date", "DLW"],
)
data.to_csv(fname, index=False)

# DSW
fname = data_dir + "meas_dsw.csv"
data = pd.DataFrame(
    np.hstack((swd_date[:, None], swd[:, None])),
    columns=["Date", "DLW"],
)
data.to_csv(fname, index=False)

# DSW clear simulated with libradtran
fname = data_dir + "libradtran_dsw_clear.csv"
data = pd.DataFrame(
    np.hstack((np.array(libdate)[:, None], np.array(swd_clear)[:, None])),
    columns=["Date", "DSW_clear"],
)
data.to_csv(fname, index=False)


# ERA5
fname = data_dir + "era5_dflux.csv"
data = pd.DataFrame(
    np.hstack(
        (
            era5_date[:, None],
            era5_swd[:, None],
            era5_swd_clr[:, None],
            era5_lwd[:, None],
            era5_lwd_clr[:, None],
        )
    ),
    columns=["Date", "DSW", "DSW_clear", "DLW", "DLW_clear"],
)
data.to_csv(fname, index=False)


# start_date = swd_date[1900]
# final_date = swd_date[-1]
# plot_measured_and_clear(
#     swd_date,
#     swd,
#     swd_clear,
#     lwd_date,
#     lwd,
#     lwd_clear_date,
#     lwd_clear,
#     SAVEFIGS,
#     out_dir + "lwd_swd_202205_10_17.eps",
# )

# plot_meas_clear_forcings(
#     swd_date,
#     swd,
#     swd_clear,
#     swd_force,
#     lwd_date,
#     lwd,
#     lwd_clear_date,
#     lwd_clear,
#     lwd_force,
#     SAVEFIGS,
#     out_dir + "lwd_swd_forcing_202205_10_17.eps",
#     start_date,
#     final_date,
# )

# plot_forcings(
#     swd_date,
#     swd,
#     swd_clear,
#     swd_force,
#     lwd_date,
#     lwd,
#     lwd_clear_date,
#     lwd_clear,
#     lwd_force,
#     start_date,
#     final_date,
#     "%d",
#     "2022/05/16",
#     SAVEFIGS,
#     out_dir + "forcing_202205_10_17.eps",
# )

# Plot the forcing for just 2022/05/16
# start_date = dt.datetime(2022, 5, 16)
# final_date = swd_date[-1]
# plot_forcings(
#     swd_date,
#     swd,
#     swd_clear,
#     swd_force,
#     lwd_date,
#     lwd,
#     lwd_clear_date,
#     lwd_clear,
#     lwd_force,
#     start_date,
#     final_date,
#     "%H",
#     "Hour on 2022/05/16",
#     SAVEFIGS,
#     out_dir + "forcing_20220516.eps",
# )


# Plot the forcing and the lidar data for 2022/05/16
#                    vmin=cbar[i,0], vmax=cbar[i,1],


tot_force = swd_force + lwd_force

start_date = dt.datetime(2022, 5, 16)
final_date = dt.datetime(2022, 5, 17)
datestring = "20220516"
# fname = f"KGI_MPLDataV2{datestring}.nc"
fname = f"KGI_MPLData_{datestring}.cdf"
titlestr = "King George Island MPL Data Processed " + datestring
rot = 0
fontsize = 12
datefmt = "%H"


timestamps, alt, pixelval, cticks, cticklabels, norm, mymap = get_lidar_data(
    lidar_dir + fname
)


# plt.figure(num=3)  # , figsize=[6.5, 5.6])
# plt.clf()
# figsize=(plotwidth, plotheight)


# 2-panel plot of lidar and cloud forcing for Bromwich et al BAMS article
fig, axs = plt.subplots(1, 1, constrained_layout=True, clear=True, num=9)

ax = axs
pcm = ax.pcolor(timestamps, alt, pixelval, norm=norm, cmap=mymap)
ax.set_xlim([0, 24])
ax.set_xticks([0, 3, 6, 9, 12, 15, 18, 21, 24])
ax.set_position((0.11, 0.6, 0.91, 0.35))
ylim = ax.set_ylim([0, 2.5])
ax.set_ylabel("Altitude (km)", fontsize=fontsize)
ax.text(0.5, 2.08, "e)", fontsize=fontsize)
cbar = fig.colorbar(pcm, ax=ax, pad=0.02, location="right", shrink=0.8)
cbar.set_ticks(cticks)
cbar.set_ticklabels(cticklabels, fontsize=fontsize)

# ax3 = axs[1]
# ax3.plot(lwd_date, lwd_force, color="orange", linewidth=3, label="LWD")
# ax3.plot(swd_date, swd_force, color="red", label="SWD")
# ax3.plot(swd_date, tot_force, color="blue", label="Total")
# # ax3.plot(era5_date, era5_lwd_force, "k.", label="ERA5 LWD")
# # ax3.plot(era5_date, era5_swd_force, ".", color="brown", label="ERA5 SWD")
# # ax3.plot(era5_date, era5_lwd_force+era5_swd_force, "ms", label="ERA5 Total")
# ax3.plot(lwd_date, np.zeros(len(lwd_date)), "k:")
# ax3.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# ax3.tick_params(axis="x", rotation=rot)
# ax3.set_xlabel("Hour on 2022/05/16", fontsize=fontsize)

# ax3.set_xlim(start_date, final_date)
# ax3.set_ylim([-105, 120])
# ax3.tick_params(axis="x", rotation=rot)
# ax3.set_ylabel("Forcing (W/m$^{2}$)", fontsize=fontsize)
# ax3.yaxis.set_label_coords(-0.08, 0.5)
# ax3.legend(loc="upper center", bbox_to_anchor=(1.093, 0.8))
# ax3.set_position((0.11, 0.1, 0.75, 0.43))
# ax3.text(dt.datetime(2022, 5, 16, 0, 40), 100, "f)", fontsize=fontsize)
# ax3.legend(bbox_to_anchor=(1., 1), loc="upper left", borderaxespad=0.0)
# ax3.legend(loc=[1.3, 0.62])

# if savefig:
#     plt.savefig("forcing.png", format="png")
#     plt.savefig("forcing.eps", format="eps")


# Figure of lidar for Bromwich et al 2024 BAMS paper
fig, ax = plt.subplots(num=9, figsize=(7, 3), clear=True, layout="constrained")

pcm = ax.pcolor(timestamps, alt, pixelval, norm=norm, cmap=mymap)
ax.set_xlim([0, 24])
ax.set_xticks((0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24))
ax.set_position((0.09, 0.15, 0.94, 0.8))
ylim = ax.set_ylim([0, 2.6])
ax.set_ylabel("Altitude (km)", fontsize=fontsize)
ax.set_xlabel("Time (UTC)")
cbar = fig.colorbar(pcm, ax=ax, pad=0.02, location="right", shrink=0.8)
cbar.set_ticks(cticks)
cbar.set_ticklabels(cticklabels, fontsize=fontsize)

if savefig:
    plt.savefig("lidar_20220516.png", format="png")
    plt.savefig("lidar_20220516.eps", format="eps")


ax.text(0.5, 2.35, "e)", fontsize=fontsize)
if savefig:
    plt.savefig("lidar_20220516_e.png", format="png")
    plt.savefig("lidar_20220516_e.eps", format="eps")


base_date = dt.datetime(2022, 5, 16)  # , tzinfo=pytz.utc)
lwd_hour = [(x - base_date).total_seconds() / 60 / 60 for x in lwd_date]
swd_hour = [(x - base_date).total_seconds() / 60 / 60 for x in swd_date]
era_hour = [(x - base_date).total_seconds() / 60 / 60 for x in era5_date]

# Difference between net and downward forcings
# plt.figure()
# plt.plot(era_hour, era5_lwnet - era5_lwnet_clr - era5_lwd_force)
# plt.plot(era_hour, era5_swnet - era5_swnet_clr - era5_swd_force)
# plt.plot(era_hour, -.16*era5_swd_force)

# plt.figure()
# plt.plot(era_hour, era5_swnet - era5_swnet_clr, label="net SW forcing")
# plt.plot(era_hour, era5_swd_force, label="SWD forcing")
# plt.legend()

# plt.figure()
# plt.plot(era_hour, era5_swd, label = "SW down")
# plt.plot(era_hour, era5_swnet - era5_swd, label = "SW up")
# plt.plot(era_hour, -era5_swd*.1, label = "-0.1*SW down")
# plt.legend()

# Figure of cloud forcing for Bromwich et al 2024 BAMS paper

# fmt: off
fig, ax = plt.subplots(num=10, figsize=(7, 3), clear = True, layout="constrained")
ax3 = ax
ax3.plot(lwd_hour, lwd_force, color="orange", linewidth=3, label="Measured, Longwave")
ax3.plot(swd_hour, swd_force, color="red", label="Measured, Shortwave")
ax3.plot(swd_hour, tot_force, color="blue", label="Measured, Total")
ax3.plot(era_hour, era5_lwd_force, ".-", color='brown', label="ERA5, Longwave")
ax3.plot(era_hour, era5_swd_force, ".-", color='pink', label="ERA5, Shortwave")
ax3.plot(era_hour, era5_swd_force + era5_lwd_force, "c.-", label="ERA5, Total")

# fmt: on

ax3.plot((0, 24), (0, 0), "k:")
ax3.set_xlabel("Time (UTC)", fontsize=fontsize)
ax3.set_xlim([0, 24])
ax3.set_ylim([-105, 120])
ax3.set_xticks((0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24))
ax3.set_ylabel("Cloud Forcing (W/m$^{2}$)", fontsize=fontsize)
ax3.yaxis.set_label_coords(-0.08, 0.5)
ax3.legend(loc="center left", bbox_to_anchor=(0, 0.37))
ax3.set_position((0.115, 0.155, 0.86, 0.8))

if savefig:
    plt.savefig("forcing.png", format="png")
    plt.savefig("forcing.eps", format="eps")

ax3.text(dt.datetime(2022, 5, 16, 0, 40), 100, "f)", fontsize=fontsize)
if savefig:
    plt.savefig("forcing_e.png", format="png")
    plt.savefig("forcing_e.eps", format="eps")


# What fraction of light was blocked in the SW?
ind1 = np.argmin(np.abs(np.array(swd_hour)))
ind2 = np.argmin(np.abs(np.array(swd_hour) - 24)) + 1
swd_tot = np.trapz(swd[ind1:ind2], swd_hour[ind1:ind2])
swd_clr_tot = np.trapz(era5_swd_clr_intp[ind1:ind2], swd_hour[ind1:ind2])
print(f"Sunlight making it through: {100*swd_tot/swd_clr_tot}%")

ind1 = np.argmin(np.abs(np.array(lwd_hour)))
ind2 = np.argmin(np.abs(np.array(lwd_hour) - 24)) + 1
lwd_tot = np.trapz(lwd[ind1:ind2], lwd_hour[ind1:ind2])
print(f"Shortwave: {swd_tot} W/m2")
print(f"Longwave: {lwd_tot} W/m2")
print(f"Fraction Longwave: {100*lwd_tot/(lwd_tot+swd_tot)}")

# ERA5 net
era_net = era5_lwnet - era5_lwnet_clr + era5_swnet - era5_swnet_clr
plt.figure()
plt.plot(era_hour, era5_lwnet - era5_lwnet_clr, label="LW net forcing")
plt.plot(era_hour, era_net, label="Net forcing")
plt.plot(era_hour, era5_swnet - era5_swnet_clr, label="SW net forcing")
plt.legend()

# # Plot the lidar data with the measured data
# fig, axs = plt.subplots(3, 1, constrained_layout=True, num=10)

# # Other formatting
# ax = axs[0]
# pcm = ax.pcolor(timestamps, alt, pixelval, norm=norm, cmap=mymap)
# ax.set_xlim([0, 24])
# ax.set_xticks([0, 3, 6, 9, 12, 15, 18, 21, 24])
# ax.set_position((0.11, 0.6, 0.914, 0.35))
# ylim = ax.set_ylim([0, 2.5])
# ax.set_ylabel("Altitude (km)", fontsize=fontsize)
# ax.text(0.5, 2.08, "e)", fontsize=fontsize)
# cbar = fig.colorbar(pcm, ax=ax, pad=0.02, location="right", shrink=0.8)
# cbar.set_ticks(cticks)
# cbar.set_ticklabels(cticklabels, fontsize=fontsize)

# start_date = dt.datetime(2022, 5, 15, 22, 0, 0)
# final_date = dt.datetime(2022, 5, 17, 13)
# ax2 = axs[1]
# ax2.plot(lwd_date, lwd, label="lwd")
# ax2.plot(lwd_clear_date, lwd_clear, "+", label="lwd clear")
# ax2.plot(lwd_date, lwd_clear_interp, label="lwd clear interp")
# ax2.set_xlim([start_date, final_date])
# ax2.legend()

# ax3 = axs[2]
# ax3.plot(era5_date, era5_swd_clr, ".", label="SWD clear ERA5")
# ax3.plot(libdate, swd_clear, label="SWD clear libradtran ")
# ax3.plot(swd_date, swd, label="SWD measured")
# ax3.legend()
# ax3.set_xlim([start_date, final_date])
# ax3.legend(bbox_to_anchor=(1., 1), loc="upper left", borderaxespad=0.0)
# ax3.legend(loc=[1.3, 0.62])


# # Cumulative forcing for SWD
# swd_tstamp = [x.timestamp() for x in swd_date]
# swd_clear_lib_interp = np.interp(swd_tstamp, lib_tstamp, swd_clear_lib)
# i1 = np.where(np.array(swd_tstamp) >= start_date.timestamp())[0][0]
# i2 = np.where(np.array(swd_tstamp) <= final_date.timestamp())[0][-1]
# swd_force = swd - swd_clear_lib_interp
# ntimes = i2 - i1
# swf_tot = np.zeros([ntimes])
# swf_tot_mwh = np.zeros([ntimes])
# count = 0
# for i in range(i1, i2):
#     swf = swd_force[i] * (swd_tstamp[i] - swd_tstamp[i - 1])
#     swf_mwh = swd_force[i] * (swd_tstamp[i] - swd_tstamp[i - 1]) / 3600 / 1000
#     swf_tot[count] = swf_tot[count - 1] + swf
#     swf_tot_mwh[count] = swf_tot_mwh[count - 1] + swf_mwh
#     count += 1
# is1 = i1
# is2 = i2


# # lwd_tot = np.interp(lwd_tstamp, lwd_clear_interp)
# # swd_tot = np.interp(lwd_tstamp, lwd_clear_interp)
# lwf_tot_sw_date = np.interp(swd_tstamp[is1:is2], lwd_tstamp[il1:il2], lwf_tot)
# lwd_force_sw_date = np.interp(swd_tstamp, lwd_tstamp, lwd_force)

# # ax1 = plt.subplot(221)
# # ax2 = plt.subplot(222, sharex=ax1)
# # ax3 = plt.subplot(223, sharex=ax1)
# # ax4 = plt.subplot(224, sharex=ax1)

# # ax2.plot(lwd_date[:ilw], lwd_force[:ilw], color="blue", label="LWD Forcing")
# # ax2.plot(swd_date[:isw], swd_force[:isw], color="red", label="SWD Forcing")
# # ax2.legend()
# # ax2.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# # ax2.tick_params(axis="x", rotation=rot)

# # ax4.plot(lwd_date[:ilw], lwf_tot[:ilw] / 1e6, color="blue", label="LWD")
# # ax4.plot(swd_date[:isw], swf_tot[:isw] / 1e6, color="red", label="SWD")
# # ax4.plot(
# #     swd_date[:isw],
# #     (swf_tot[:isw] + lwf_tot_sw_date[:isw]) / 1e6,
# #     "k",
# #     label="Total",
# # )
# # ax4.legend()
# # ax4.set_ylabel("Cumulative Forcing (MJ/m$^{2}$)")
# # ax4.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# # ax4.set_xlabel("Day of 2022/02")
# # ax4.tick_params(axis="x", rotation=rot)

# totf = swf_tot + lwf_tot_sw_date
# totfstr = str(round(totf[-1] / 1e6))

# plot_measured_and_clear()
# # plot_meas_clear_forcings_cumulative_shorttime()
# plot_meas_clear_forcings()

# if plotfigs:
#     # .. Plot it
#     plt.figure(1)
#     plt.clf()
#     plt.plot(dtimes, rad, label="pyranometer")
#     plt.plot(dtimes, diffuse, label="pysolar")
#     plt.plot(dtimes, diffuse * diffuse_fac, label="clear")

#     # .. Cloud forcing
#     sw_forcing = diffuse * diffuse_fac - rad
#     plt.plot(dtimes, sw_forcing, label="forcing")
#     plt.legend()
