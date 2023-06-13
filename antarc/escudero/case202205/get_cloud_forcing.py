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
  2*) Dompute clear-sky SWD flux using get_libradtran (done here)
  3*) Subtract cloudy and clear-sky fluxes to get SW forcing (done here)

E) miniMPL: Please see the repos mpl2nc and mpl-processing.

"""

# Dependencies
import datetime as dt
import numpy as np

# My modules
from antarc.run_radtran import get_clear_fluxes
from antarc.escudero.get_lwd_flux import get_lwd, get_lwd_clear
from antarc.escudero.get_atm_profs import get_atm_profs
from antarc.escudero.get_swd_flux import get_obs_swd, get_libradtran
from antarc.escudero.get_era5_from_web import get_era5_from_web
from antarc.escudero.case202205.make_plots import plot_measured_and_clear
from antarc.escudero.case202205.make_plots import plot_forcings
from antarc.escudero.case202205.make_plots import plot_meas_clear_forcings
from antarc.escudero.case202205.make_plots import plot_forcings_and_lidar


# Parameter modules specific to the user: the user must create params.py
# in the directory antarc with the variables MEAS_DIR and PROJ_DIR set
from antarc import params

# Parameter modules
from antarc.escudero.parameters import esc202205 as esc_case
from antarc.escudero.parameters import esc_params, atm_prof_params
from antarc.escudero.parameters import swd_params, lwd_params, radtran_params

utc = dt.timezone.utc

# TODO: get this info from the parameter file
MEAS_DIR = params.MEAS_DIR
dir_lwdclear = MEAS_DIR + "Escudero/lwd_clear/"
out_dir = params.PROJECT_DIR + "case_studies/May_2022/figures/"

# Run flags
SAVEFIGS = True
REDO = False  # redo calculations even if they already exist
# redo_profs = False

date1 = esc_case.DATE1
date2 = esc_case.DATE2


# A5) Get ERA5 data if it is not found or if REDO=TRUE
get_era5_from_web(esc_params, esc_case, REDO)

# C) Longwave
# C2) Create the atmospheric profiles
get_atm_profs(atm_prof_params, esc_params, esc_case, REDO)


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

# C5) Subtract cloudy and clear-sky fluxes to get LW forcing
lwd_force = lwd - lwd_clear_interp


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
libdate, swd_clear = get_libradtran(esc_params, swd_date, clrfile, REDO)
lib_tstamp = [x.timestamp() for x in libdate]
swd_tstamp = [x.timestamp() for x in swd_date]
swd_clear = np.interp(swd_tstamp, lib_tstamp, swd_clear)

# D3) Subtract cloudy and clear-sky fluxes to get SW forcing
swd_force = swd - swd_clear

start_date = swd_date[1900]
final_date = swd_date[-1]


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
start_date = dt.datetime(2022, 5, 16)
final_date = swd_date[-1]
plot_forcings_and_lidar(
    swd_date,
    swd,
    swd_clear,
    swd_force,
    lwd_date,
    lwd,
    lwd_clear_date,
    lwd_clear,
    lwd_force,
    start_date,
    final_date,
    "%H",
    "Hour on 2022/05/16",
    SAVEFIGS,
    out_dir + "forcing_20220516.eps",
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
