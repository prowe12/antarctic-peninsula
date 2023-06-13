#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:25 2022

@author: prowe
"""

import os

import numpy as np
import datetime as dt
from matplotlib import pyplot as plt
from netCDF4 import Dataset

# Parameter modules
from antarc import params
from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import atm_prof_params as esc_atm_prof_params
from antarc.rothera.parameters import roth_params, roth_atm_prof_params
from antarc.marambio.parameters import mbio_atm_prof_params, mbio_params

# Modules
from antarc.pwrf_v_sonde.get_era_profs import get_era
from antarc.pwrf_v_sonde.load_pwrf import LoadPwrf, Pwrf
from antarc.pwrf_v_sonde.load_sonde import Sonde, IgraSonde
from antarc.pwrf_v_sonde.load_sonde import DenialSonde, DenialSondePR

# from antarc.get_atm_profs_rsrc import FileInfo, Era
# from antarc.get_era5_press_levs_from_web import get_era5_press_levs


# # # # # #   ERA5   # # # # #
# def get_era(prefix, params, atm_prof_params, date_str, date_fmt):
#     lat = params.LATITUDE
#     lon = params.LONGITUDE
#     direc = atm_prof_params.ERA_DIR
#     era_fileformat = atm_prof_params.ERA_FILEFORMAT

#     # Get the ERA5 files form the web if needed
#     dtime = dt.datetime.strptime(date_str, date_fmt)
#     date1 = dt.datetime(dtime.year, dtime.month, dtime.day)
#     get_era5_press_levs(lat, lon, date1, date1, direc, prefix, False)

#     Files = FileInfo(direc, era_fileformat)
#     # Get the ERA5 data
#     return Era(Files, era_fileformat, dtime, lat, lon)


# Base directories
meas_dir = params.MEAS_DIR
main_dir = params.PROJECT_DIR + "case_studies/"
era5_dir = meas_dir + "Escudero/era5/"
save_dir = main_dir

# Flag for whether to save the resulting figure
# See the if savefig block to see where the figure will be saved to
# (usually savedir)
savefig = True


pwrfdirs = [
    "Dec_2018/pwrf_stations/",
    "Dec_2018/pwrf_stations/",
    "Dec_2018/pwrf_stations/",
    "Feb_2022/pwrf_stations/",
    "Feb_2022/pwrf_stations/",
]
pwrf_dates = [
    "20181205_08",
    "20181205_08",
    "20181205_08",
    "20220206_09",
    "20220206_09",
]
pwrf_date_fmts = [
    "%Y%m%d_08",
    "%Y%m%d_08",
    "%Y%m%d_08",
    "%Y%m%d_09",
    "%Y%m%d_09",
]

# Modify idate from 0 to 4 to select the plot to create
idate = 4
pwrfdir = main_dir + pwrfdirs[idate]
pwrf_date_fmt = pwrf_date_fmts[idate]
pwrf_date = pwrf_dates[idate]

# Upwind: Either Escudero or Rothera
upwnd_dirs = [
    meas_dir + "Rothera/radiosonde/",
    meas_dir + "Rothera/radiosonde/",
    meas_dir + "Rothera/igra/derived/",
    meas_dir + "Rothera/igra/derived/",
    meas_dir + "Escudero/radiosondes/datadenial/v0/",
]
upwnd_snd_files = [
    "Rothera_2018120511.dat",  # idate = 0
    "Rothera_2018120611.dat",  # idate = 1; Figure S3b
    "Rothera_igra_20181207_12.txt",  # idate = 2
    "Rothera_igra_20220209_12.txt",  # idate = 3
    "esc_sonde_dd_v0_2022020823.txt",  # idate = 4; Figure S3a
]
upwnd_sndfile_fmts = [
    "Rothera_%Y%m%d%H.dat",
    "Rothera_%Y%m%d%H.dat",
    "Rothera_igra_%Y%m%d_%H.txt",
    "Rothera_igra_%Y%m%d_%H.txt",
    "esc_sonde_dd_v0_%Y%m%d%H.txt",
]
upwnd_date_strs = [
    "2018/12/05 12 UT",
    "2018/12/06 12 UT",
    "2018/12/07 12 UT",
    "2022/02/09 12 UT",
    "2022/02/08 23 UT",
]
upwnd_snd_classes = [
    DenialSonde,
    DenialSonde,
    IgraSonde,
    IgraSonde,
    DenialSondePR,
]
upwnd_locs = ["Rothera", "Rothera", "Rothera", "Rothera", "Escudero"]


# Leeside: Marambio
mbio_dirs = [
    "Marambio/radiosonde/",
    "Marambio/radiosonde/",
    "Marambio/radiosonde/",
    "Marambio/radiosonde/",
    "Marambio/radiosonde/",
]
mbio_snd_files = [
    "Marambio_2018120511.dat",
    "Marambio_2018120511.dat",
    "Marambio_2018120511.dat",
    "Marambio_2022020712.nc",
    "Marambio_2022020712.nc",
]
mbio_date_strs = [
    "2018/12/05 12 UT",
    "2018/12/05 12 UT",
    "2018/12/05 12 UT",
    "2022/02/07 12 UT",
    "2022/02/07 12 UT",
]
mbio_sndfile_fmts = [
    "Marambio_%Y%m%d%H.dat",
    "Marambio_%Y%m%d%H.dat",
    "Marambio_%Y%m%d%H.dat",
    "Marambio_%Y%m%d%H.nc",
    "Marambio_%Y%m%d%H.nc",
]
mbio_snd_classes = [DenialSonde, DenialSonde, DenialSonde, Sonde, Sonde]


dhour = 1.0


mbio_dir = meas_dir + mbio_dirs[idate]
mbio_sndfile_fmt = mbio_sndfile_fmts[idate]
mbio_snd_file = mbio_snd_files[idate]
mbio_sndfile_fmt = mbio_sndfile_fmts[idate]
mbio_date_str = mbio_date_strs[idate]
mbioSndClass = mbio_snd_classes[idate]
upwnd_loc = upwnd_locs[idate]

# The radiosonde data
mbioSnd = mbioSndClass(mbio_dir, mbio_snd_file, mbio_sndfile_fmt)


# Files and dir for upwind (can be either Escudero or Rothera)
upwnd_dir = upwnd_dirs[idate]
upwnd_snd_file = upwnd_snd_files[idate]
upwnd_sndfile_fmt = upwnd_sndfile_fmts[idate]
upwnd_datestr = upwnd_date_strs[idate]
upwndSndClass = upwnd_snd_classes[idate]


# ERA5 interpolated to location: upwind (Escudero or Rothera)
fmt = "%Y/%m/%d %H UT"
# Location
era_pfx = {"Escudero": "era5_esc_", "Rothera": "era5_roth_"}
station_params = {
    "Escudero": [esc_params, esc_atm_prof_params],
    "Rothera": [roth_params, roth_atm_prof_params],
}
pfx = era_pfx[upwnd_locs[idate]]
sta_params, sta_atm_prof_params = station_params[upwnd_locs[idate]]

upwndEra = get_era(pfx, sta_params, sta_atm_prof_params, upwnd_datestr, fmt)
# lat = esc_params.LATITUDE
# lon = esc_params.LONGITUDE
# direc = esc_atm_prof_params.ERA_DIR
# era_fileformat = esc_atm_prof_params.ERA_FILEFORMAT
# Files = FileInfo(direc, era_fileformat)
# era_fileformat = esc_atm_prof_params.ERA_FILEFORMAT
# dtime = dt.datetime.strptime(upwnd_datestr, "%Y/%m/%d %H UT")
# Get the ERA5 files form the web if needed
# date1 = dt.datetime(dtime.year, dtime.month, dtime.day)
# get_era5_press_levs(lat, lon, date1, date1, direc, pfx, False)
# Get the ERA5 data
# eraUpwind = Era(Files, era_fileformat, dtime, lat, lon)

# ERA5 interpolated to location: leeside (Marambio)
pfx = "era5_mbio_"
mbioEra = get_era(pfx, mbio_params, mbio_atm_prof_params, mbio_date_str, fmt)
# lat = mbio_params.LATITUDE
# lon = mbio_params.LONGITUDE
# direc = mbio_atm_prof_params.ERA_DIR
# era_fileformat = mbio_atm_prof_params.ERA_FILEFORMAT
# Files = FileInfo(mbio_dir, era_fileformat)
# era_fileformat = esc_atm_prof_params.ERA_FILEFORMAT
# dtime = dt.datetime.strptime(mbio_date_str, "%Y/%m/%d %H UT")
# # Get the ERA5 files form the web if needed
# date1 = dt.datetime(dtime.year, dtime.month, dtime.day)
# get_era5_press_levs(lat, lon, date1, date1, direc, pfx, False)
# # Get the ERA5 data
# eraMbio = Era(Files, era_fileformat, dtime, lat, lon)

# # # # # # # # # # # # # # #


# Create the PWRF paths
loc = "Marambio"
tc_file = pwrfdir + "PWRF_D03_" + loc + "_tc_" + pwrf_date + "_hly.nc"
rh_file = pwrfdir + "PWRF_D03_" + loc + "_rh_" + pwrf_date + "_hly.nc"
wnd_file = pwrfdir + "PWRF_D03_" + loc + "_uvmet_" + pwrf_date + "_hly.nc"

# The PWRF data for Marambio, trimmed to remove heights outside the
# radiosounding, and including the index to the PWRF time of interest
starttime = dt.datetime.strptime(pwrf_date, pwrf_date_fmt)
allPwrf = LoadPwrf(tc_file, rh_file, wnd_file, starttime, dhour)
mbioPwrf = Pwrf(allPwrf, mbioSnd.dtime)
mbioPwrf.trim(mbioSnd.lev)

# Get difference statistics: max, min, mean, rms difference
mbio_dtemp, mbio_drh, mbio_dwind = mbioSnd.get_diffs(
    mbioPwrf.level, mbioPwrf.temp, mbioPwrf.rh, mbioPwrf.wind
)


# Create the PWRF paths for the upwind location
loc = upwnd_loc
tc_file = pwrfdir + "PWRF_D03_" + loc + "_tc_" + pwrf_date + "_hly.nc"
rh_file = pwrfdir + "PWRF_D03_" + loc + "_rh_" + pwrf_date + "_hly.nc"
wnd_file = pwrfdir + "PWRF_D03_" + loc + "_uvmet_" + pwrf_date + "_hly.nc"

# The radiosonde data
upwndSnd = upwndSndClass(upwnd_dir, upwnd_snd_file, upwnd_sndfile_fmt)

# The PWRF data for Rothera, including the index to the PWRF time of interest
allPwrf = LoadPwrf(tc_file, rh_file, wnd_file, starttime, dhour)
upwndPwrf = Pwrf(allPwrf, upwndSnd.dtime)
upwndPwrf.trim(upwndSnd.lev)

# Get difference statistics: max, min, mean, rms difference
upwnd_dtemp, upwnd_drh, upwnd_dwind = upwndSnd.get_diffs(
    upwndPwrf.level, upwndPwrf.temp, upwndPwrf.rh, upwndPwrf.wind
)


# What we will plot
# ax1.semilogy(upwndPwrf.temp, upwndPwrf.level, label="PWRF")
# ax1.semilogy(upwndSnd.temp, upwndSnd.lev, label="Radiosonde")

# ax2.semilogy(mbioPwrf.temp, mbioPwrf.level, label="PWRF")
# ax2.semilogy(mbioSnd.temp, mbioSnd.lev, label="Radiosonde")

# ax3.semilogy(upwndPwrf.rh, upwndPwrf.level, label="PWRF")
# ax3.semilogy(upwndSnd.rh, upwndSnd.lev, label="Radiosonde")
# ax4.semilogy(mbioPwrf.rh, mbioPwrf.level, label="PWRF")
# ax4.semilogy(mbioSnd.rh, mbioSnd.lev, label="Radiosonde")
# ax5.semilogy(upwndPwrf.wind, upwndPwrf.level, label="PWRF")
# ax5.semilogy(upwndSnd.wspeed[iwnd], upwndSnd.lev[iwnd], label="Radiosonde")
# ax6.semilogy(mbioPwrf.wind, mbioPwrf.level, label="PWRF")
# ax6.semilogy(mbioSnd.wspeed, mbioSnd.lev, label="Radiosonde")


# Make the figure with winds
plt.close(3)
plt.figure(num=3, figsize=[6, 6])
ax1 = plt.subplot(321)
ax2 = plt.subplot(322, sharey=ax1)
ax3 = plt.subplot(323, sharey=ax1)
ax4 = plt.subplot(324, sharey=ax1)
ax5 = plt.subplot(325, sharey=ax1)
ax6 = plt.subplot(326, sharey=ax1)

wid = 0.37
hit = 0.24
ax1.set_position([0.12, 0.72, wid, hit])
ax2.set_position([0.6, 0.72, wid, hit])
ax3.set_position([0.12, 0.4, wid, hit])
ax4.set_position([0.6, 0.4, wid, hit])
ax5.set_position([0.12, 0.08, wid, hit])
ax6.set_position([0.6, 0.08, wid, hit])

ax1.semilogy(upwndPwrf.temp, upwndPwrf.level, label="PWRF")
ax1.semilogy(upwndSnd.temp, upwndSnd.lev, label="Radiosonde")
ax1.semilogy(upwndEra.T - 273.15, upwndEra.P, label="ERA5")
ax1.set_ylim([1000, 100])
ax1.set_xlabel("Temperature ($^o$C)")
ax1.set_xlim([-65, 20])
# ax1.legend()
ax1.text(-60, 850, "a")
ax1.set_title(upwnd_loc + " " + upwnd_datestr)

ax2.semilogy(mbioPwrf.temp, mbioPwrf.level, label="PWRF")
ax2.semilogy(mbioSnd.temp, mbioSnd.lev, label="Radiosonde")
ax2.semilogy(mbioEra.T - 273.15, mbioEra.P, label="ERA5")
ax2.set_ylim([1000, 100])
ax2.set_xlim([-65, 20])
ax2.set_xlabel("Temperature ($^o$C)")
ax2.legend()
ax2.text(-60, 850, "d")
ax2.set_title("Marambio " + mbio_date_str)

ax3.semilogy(upwndPwrf.rh, upwndPwrf.level, label="PWRF")
ax3.semilogy(upwndSnd.rh, upwndSnd.lev, label="Radiosonde")
ax3.semilogy(upwndEra.rh, upwndEra.P, label="ERA5")
ax3.set_ylim([1000, 100])
ax3.set_xlabel("RH (%)")
ax3.set_xlim([0, 120])
ax3.text(5, 850, "b")

upwnd_era_windspeed = np.sqrt(upwndEra.uwind**2 + upwndEra.vwind**2)
iwnd = np.where(~np.isnan(upwndSnd.wspeed))[0]
ax5.semilogy(upwndPwrf.wind, upwndPwrf.level, label="PWRF")
ax5.semilogy(upwndSnd.wspeed[iwnd], upwndSnd.lev[iwnd], label="Radiosonde")
ax5.semilogy(upwnd_era_windspeed, upwndEra.P, label="ERA5")
# TODO: did not get ERA5 wind speed
# ax5.semilogy(upwndEra.rh, upwndEra.P, label="ERA5")
ax5.set_ylim([1000, 100])
# ax5.set_xlim([0, 50])
ax5.set_xlabel("Wind speed (m/s)")
ax5.text(1.5, 800, "c")

ax4.semilogy(mbioPwrf.rh, mbioPwrf.level, label="PWRF")
ax4.semilogy(mbioSnd.rh, mbioSnd.lev, label="Radiosonde")
ax4.semilogy(mbioEra.rh, mbioEra.P, label="ERA5")
ax4.set_ylim([1000, 100])
ax4.set_xlabel("RH (%)")
ax4.text(5, 850, "e")
ax4.set_xlim([0, 120])

mbio_era_windspeed = np.sqrt(mbioEra.uwind**2 + mbioEra.vwind**2)
ax6.semilogy(mbioPwrf.wind, mbioPwrf.level, label="PWRF")
ax6.semilogy(mbioSnd.wspeed, mbioSnd.lev, label="Radiosonde")
ax6.semilogy(mbio_era_windspeed, mbioEra.P, label="ERA5")
ax6.set_ylim([1000, 100])
ax6.set_xlabel("Wind speed (m/s)")
ax6.text(1.5, 800, "f")
# ax6.set_xlim([0, 25])

ax1.set_ylabel("Pressure (hPa)")
ax3.set_ylabel("Pressure (hPa)")
ax5.set_ylabel("Pressure (hPa)")


ax1.set_yticks(
    [1000, 600, 400, 300, 200, 100],
    labels=["1000", "600", "400", "300", "200", "100"],
)

if savefig:
    figname = "era_pwrf_v_sonde_" + pwrf_date + ".png"
    plt.savefig(save_dir + figname, dpi=300, format="png")


print()
print(upwnd_snd_file)
print(f'Max T: {max(abs(upwnd_dtemp["max"]),abs(upwnd_dtemp["min"]))}')
print(f"Mean T: {upwnd_dtemp['mean']}")
print(f"RMS T: {upwnd_dtemp['rms']}")

print(upwnd_dwind["max"])
print(upwnd_dwind["min"])
print(upwnd_dwind["rms"])

print(f'Max RH: {max(abs(upwnd_drh["max"]),abs(upwnd_drh["min"]))}')
print(f"Mean RH: {upwnd_drh['mean']}")
print(f'RMS RH: {upwnd_drh["rms"]}')


print()
print(mbio_snd_file)
print(mbio_dtemp["max"])
print(mbio_dtemp["min"])
print(mbio_dtemp["rms"])

print(mbio_dwind["max"])
print(mbio_dwind["min"])
print(mbio_dwind["rms"])

print(mbio_drh["max"])
print(mbio_drh["min"])
print(mbio_drh["rms"])
