#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:25 2022

@author: prowe
"""

import numpy as np
import datetime as dt
from matplotlib import pyplot as plt

from load_pwrf import LoadPwrf, Pwrf
from load_sonde import Sonde, IgraSonde, DenialSonde, DenialSondePR


# Base directories
meas_dir = "/Users/prowe/Sync/measurements/"
main_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/"

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
    "Rothera_2018120511.dat",
    "Rothera_2018120611.dat",
    "Rothera_igra_20181207_12.txt",
    "Rothera_igra_20220209_12.txt",
    "esc_sonde_dd_v0_2022020823.txt",
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


# Files and dir
upwnd_dir = upwnd_dirs[idate]
upwnd_snd_file = upwnd_snd_files[idate]
upwnd_sndfile_fmt = upwnd_sndfile_fmts[idate]
upwnd_date_str = upwnd_date_strs[idate]
upwndSndClass = upwnd_snd_classes[idate]

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
mbio_dtemp, mbio_drh, mbio_dwind = mbioSnd.get_diffs(mbioPwrf)


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
upwnd_dtemp, upwnd_drh, upwnd_dwind = upwndSnd.get_diffs(upwndPwrf)

# Make the figure with winds
plt.close(3)
plt.figure(num=3)
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

ax = ax1
ax.semilogy(upwndPwrf.temp, upwndPwrf.level, label="PWRF")
ax.semilogy(upwndSnd.temp, upwndSnd.lev, label="Radiosonde")
ax.set_ylim([1000, 100])
ax.set_xlabel("Temperature ($^o$C)")
ax.set_xlim([-65, 20])
# ax.legend()
ax.text(-60, 850, "a")
ax.set_title(upwnd_loc + " " + upwnd_date_str)

ax = ax3
ax.semilogy(upwndPwrf.rh, upwndPwrf.level, label="PWRF")
ax.semilogy(upwndSnd.rh, upwndSnd.lev, label="Radiosonde")
ax.set_ylim([1000, 100])
ax.set_xlabel("RH (%)")
ax.set_xlim([0, 100])
ax.text(5, 850, "b")

iwnd = np.where(~np.isnan(upwndSnd.wspeed))[0]
ax = ax5
ax.semilogy(upwndPwrf.wind, upwndPwrf.level, label="PWRF")
ax.semilogy(upwndSnd.wspeed[iwnd], upwndSnd.lev[iwnd], label="Radiosonde")
ax.set_ylim([1000, 100])
# ax6.set_xlim([0, 50])
ax.set_xlabel("Wind speed (m/s)")
ax.text(1.5, 800, "c")

ax2.semilogy(mbioPwrf.temp, mbioPwrf.level, label="PWRF")
ax2.semilogy(mbioSnd.temp, mbioSnd.lev, label="Radiosonde")
ax2.set_ylim([1000, 100])
ax2.set_xlim([-65, 20])
ax2.set_xlabel("Temperature ($^o$C)")
ax2.legend()
ax2.text(-60, 850, "d")
ax2.set_title("Marambio " + mbio_date_str)

ax = ax4
ax.semilogy(mbioPwrf.rh, mbioPwrf.level, label="PWRF")
ax.semilogy(mbioSnd.rh, mbioSnd.lev, label="Radiosonde")
ax.set_ylim([1000, 100])
ax.set_xlabel("RH (%)")
ax.text(5, 850, "e")
ax.set_xlim([0, 100])

ax = ax6
ax.semilogy(mbioPwrf.wind, mbioPwrf.level, label="PWRF")
ax.semilogy(mbioSnd.wspeed, mbioSnd.lev, label="Radiosonde")
ax.set_ylim([1000, 100])
ax.set_xlabel("Wind speed (m/s)")
ax.text(1.5, 800, "f")
# ax.set_xlim([0, 25])

ax1.set_ylabel("Pressure (hPa)")
ax3.set_ylabel("Pressure (hPa)")
ax5.set_ylabel("Pressure (hPa)")


ax1.set_yticks(
    [1000, 600, 400, 300, 200, 100],
    labels=["1000", "600", "400", "300", "200", "100"],
)


print()
print(upwnd_snd_file)
print(f'Max: {max(abs(upwnd_dtemp["max"]),abs(upwnd_dtemp["min"]))}')
print(f"Mean: {upwnd_dtemp['mean']}")
print(f"RMS: {upwnd_dtemp['rms']}")

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
