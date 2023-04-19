#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 08:56:46 2022

@author: prowe

Purpose: Plot broadband radiation at Escudero on 2022/02
"""

# Dependencies
import numpy as np
import datetime as dt
import pandas as pd
from os.path import exists
import os

# My modules
from antarc.load_rad_flux import load_rad_flux
from antarc.run_radtran import get_clear_fluxes
from antarc.escudero.get_lwd_flux import get_lwd, get_lwd_clear
from antarc.escudero.get_atm_profs import get_atm_profs
from antarc.escudero.get_swd_flux import get_obs_swd, get_pysolar_swd
from antarc.escudero.get_swd_flux import get_libradtran
from antarc.escudero.get_era5_from_web import get_era5_from_web
from antarc.escudero.case202205.make_plots import plot_measured_and_clear

# Parameter modules
from antarc.escudero.parameters import esc202205 as esc_case
from antarc.escudero.parameters import esc_params, atm_prof_params
from antarc.escudero.parameters import swd_params, lwd_params, radtran_params

# from antarc.escudero.parameters import radtran_params
# Parameter modules
# from antarc.escudero.parameters import esc_params, swd_params
# from antarc.escudero.parameters import esc202205 as esc_case


# TODO: get this info from the parameter file
dir_lwdclear = "/Users/prowe/sync/measurements/Escudero/lwd_clear/"
out_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/May_2022/figures/"
fig1name = out_dir + "lwd_swd.png"
fig2name = out_dir + "lwd_swd_zoom.png"

utc = dt.timezone.utc

# Run flags
savefigs = False

date1 = esc_case.DATE1
date2 = esc_case.DATE2

# B4) Get ERA5 data
get_era5_from_web(esc_params, esc_case, False)

# C) Longwave
# C2) Create the atmospheric profiles
get_atm_profs(atm_prof_params, esc_params, esc_case)

# C3)
# Check for the clear-sky lwd fluxes and if the files have not been made,
# create them. We expect one per profile
# dir_prof = "/Users/prowe/sync/measurements/Escudero/profiles/"
# samplefile = "prof20220514_1205.nc"
# fmt = "prof%Y%m%d_%H%M.nc"
# lwdfmt = "flux%Y%m%d_%H%M_1500_2222.txt"
# allfiles = os.listdir(dir_prof)
# files = []
# for file in allfiles:
#     if len(file) == len(samplefile) and file[:4] == samplefile[:4]:
#         pfile = dt.strptime(file, fmt).stftime(lwdfmt)
#         lwdfile = f"{dir_lwdclear}{pfile}"
get_clear_fluxes(dir_lwdclear, radtran_params, (date1, date2), 3, False)

# C4) Load in measd LWD and  clear sky files created in previous step
lwd_date, lwd = get_lwd(lwd_params, date1, date2)
lwd_clear_date, lwd_clear = get_lwd_clear(esc_case, lwd_params)
lwd_tstamp = [x.timestamp() for x in lwd_date]
lwd_clear_tstamp = [x.timestamp() for x in lwd_clear_date]
lwd_clear = np.array(lwd_clear)
lwd_clear_interp = np.interp(lwd_tstamp, lwd_clear_tstamp, lwd_clear)

# C5) Subtract cloudy and clear-sky fluxes to get LW forcing
lwd_force = lwd - lwd_clear_interp

# # Cumulative forcing for LWD
# # dt.datetime(fdate.year, fdate.month, fdate.day, fdate.hour + 1, 0)
# i1 = np.where(np.array(lwd_tstamp) >= start_date.timestamp())[0][0]
# i2 = np.where(np.array(lwd_tstamp) <= final_date.timestamp())[0][-1] + 1
# ntimes = i2 - i1
# lwf_tot = np.zeros([ntimes])
# lwf_tot_mwh = np.zeros([ntimes])
# count = 0
# for i in range(i1, i2):
#     lwf = lwd_force[i] * (lwd_tstamp[i] - lwd_tstamp[i - 1])
#     lwf_mwh = lwd_force[i] * (lwd_tstamp[i] - lwd_tstamp[i - 1]) / 3600 / 1000
#     lwf_tot[count] = lwf_tot[count - 1] + lwf
#     lwf_tot_mwh[count] = lwf_tot_mwh[count - 1] + lwf_mwh
#     count += 1
# il1 = i1
# il2 = i2


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
# swd_date_long, swd_long = load_rad_flux(swd_dir, swd_fmt, date1, date2)

# D2) If needed, calculate the libradtran results and save to file
# # _, diffuse, diffuse_fac = get_pysolar_swd(swd_params, esc_params, swd_date)
# # swd_date, swd, swd_clear_lib, swd_clear, diffuse, diffuse_fac = get_swd()
redo = False
libdate, swd_clear = get_libradtran(esc_params, swd_date, clrfile, redo)
lib_tstamp = [x.timestamp() for x in libdate]
swd_tstamp = [x.timestamp() for x in swd_date]
swd_clear = np.interp(swd_tstamp, lib_tstamp, swd_clear)

# D3) Subtract cloudy and clear-sky fluxes to get SW forcing
swd_force = swd - swd_clear

start_date = swd_date[0]
final_date = swd_date[-1]
savef = False
fname = ""

plot_measured_and_clear(
    swd_date,
    swd,
    swd_clear,
    lwd_date,
    lwd,
    lwd_clear_date,
    lwd_clear,
    savef,
    fname,
)


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
