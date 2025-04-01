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
from load_sonde import DenialSondePR


# Base directories
meas_dir = "/Users/prowe/Sync/measurements/"
main_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/"

idate = 0  # 0 2 4 6
year = 2022
dhour = 1.0

if year == 2018:
    pwrfdir = "Dec_2018/pwrf_stations/"
    pwrf_date = "20181205_08"
    pwrf_date_fmt = "%Y%m%d_08"
    esc_dir = meas_dir + "Escudero/radiosondes/datadenial/v0/"
    version = 0
    esc_snd_files = [
        f"esc_sonde_dd_v{version}_2018120523.txt",
        f"esc_sonde_dd_v{version}_2018120612.txt",
        f"esc_sonde_dd_v{version}_2018120614.txt",
        f"esc_sonde_dd_v{version}_2018120617.txt",
        f"esc_sonde_dd_v{version}_2018120623.txt",
        f"esc_sonde_dd_v{version}_2018120723.txt",
        f"esc_sonde_dd_v{version}_2018120823.txt",
        f"esc_sonde_dd_v{version}_2018120923.txt",
    ]
    esc_date_strs = [
        "2018/12/05 23 UT",
        "2018/12/06 12 UT",
        "2018/12/06 14 UT",
        "2018/12/06 17 UT",
        "2018/12/06 23 UT",
        "2018/12/07 23 UT",
        "2018/12/08 23 UT",
        "2018/12/09 23 UT",
    ]
    esc_sndfile_fmt = "esc_sonde_dd_v0_%Y%m%d%H.txt"


else:
    pwrfdir = "Feb_2022/pwrf_stations/"
    pwrf_date = "20220206_09"
    pwrf_date_fmt = "%Y%m%d_09"
    esc_dir = f"{meas_dir}Escudero/radiosondes/datadenial/{year}/"
    version = 2

    esc_snd_files = [
        f"esc_sonde_dd_v{version}_2022020700.txt",
        f"esc_sonde_dd_v{version}_2022020712.txt",
        f"esc_sonde_dd_v{version}_2022020800.txt",
        f"esc_sonde_dd_v{version}_2022020823.txt",
        f"esc_sonde_dd_v{version}_2022020912.txt",
    ]

    esc_date_strs = [
        "2022/02/07 00 UT",
        "2022/02/07 12 UT",
        "2022/02/08 00 UT",
        "2022/02/08 23 UT",
        "2022/02/09 12 UT",
    ]
    esc_sndfile_fmt = f"esc_sonde_dd_v{version}_%Y%m%d%H.txt"


pwrfdir = main_dir + pwrfdir
esc_snd_file1 = esc_snd_files[idate]
esc_date_str1 = esc_date_strs[idate]
esc_snd_file2 = esc_snd_files[idate + 1]
esc_date_str2 = esc_date_strs[idate + 1]

# The radiosonde data
escSnd1 = DenialSondePR(esc_dir, esc_snd_file1, esc_sndfile_fmt)
escSnd2 = DenialSondePR(esc_dir, esc_snd_file2, esc_sndfile_fmt)


# Create the PWRF paths
loc = "Escudero"
tc_file = pwrfdir + "PWRF_D03_" + loc + "_tc_" + pwrf_date + "_hly.nc"
rh_file = pwrfdir + "PWRF_D03_" + loc + "_rh_" + pwrf_date + "_hly.nc"
wnd_file = pwrfdir + "PWRF_D03_" + loc + "_uvmet_" + pwrf_date + "_hly.nc"

# The PWRF data for Marambio, trimmed to remove heights outside the
# radiosounding, and including the index to the PWRF time of interest
starttime = dt.datetime.strptime(pwrf_date, pwrf_date_fmt)
escPwrfAll = LoadPwrf(tc_file, rh_file, wnd_file, starttime, dhour)
escPwrf1 = Pwrf(escPwrfAll, escSnd1.dtime)
escPwrf1.trim(escSnd1.lev)

escPwrf2 = Pwrf(escPwrfAll, escSnd2.dtime)
escPwrf2.trim(escSnd2.lev)

# Make the figure with winds
plt.figure(num=3, clear=True)
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

ax1.semilogy(escPwrf1.temp, escPwrf1.level, label="PWRF")
ax1.semilogy(escSnd1.temp, escSnd1.lev, label="Radiosonde")
ax1.set_ylim([1000, 100])
ax1.set_xlim([-65, 20])
ax1.set_ylabel("Pressure (hPa)")
ax1.set_xlabel("Temperature ($^o$C)")
ax1.legend()
ax1.text(-60, 850, "a")
ax1.set_title(loc + " " + esc_date_str1)

ax3.semilogy(escPwrf1.rh, escPwrf1.level, label="PWRF")
ax3.semilogy(escSnd1.rh, escSnd1.lev, label="Radiosonde")
ax3.set_ylim([1000, 100])
ax3.set_ylabel("Pressure (hPa)")
ax3.set_xlabel("RH (%)")
ax3.text(5, 850, "b")
ax3.set_xlim([0, 110])

ax5.semilogy(escPwrf1.wind, escPwrf1.level, label="PWRF")
ax5.semilogy(escSnd1.wspeed, escSnd1.lev, label="Radiosonde")
ax5.set_ylim([1000, 100])
ax5.set_ylabel("Pressure (hPa)")
ax5.set_xlabel("Wind speed (m/s)")
ax5.text(1.5, 800, "c")
# ax5.set_xlim([0, 25])

ax2.semilogy(escPwrf2.temp, escPwrf2.level, label="PWRF")
ax2.semilogy(escSnd2.temp, escSnd2.lev, label="Radiosonde")
ax2.set_ylim([1000, 100])
ax2.set_xlabel("Temperature ($^o$C)")
ax2.set_xlim([-65, 20])
# ax2.legend()
ax2.text(-60, 850, "d")
ax2.set_title(loc + " " + esc_date_str2)

ax4.semilogy(escPwrf2.rh, escPwrf2.level, label="PWRF")
ax4.semilogy(escSnd2.rh, escSnd2.lev, label="Radiosonde")
ax4.set_ylim([1000, 100])
ax4.set_xlabel("RH (%)")
ax4.set_xlim([0, 110])
ax4.text(5, 850, "e")

iwnd = np.where(~np.isnan(escSnd2.wspeed))[0]
ax6.semilogy(escPwrf2.wind, escPwrf2.level, label="PWRF")
ax6.semilogy(escSnd2.wspeed[iwnd], escSnd2.lev[iwnd], label="Radiosonde")
ax6.set_ylim([1000, 100])
# ax6.set_xlim([0, 50])
ax6.set_xlabel("Wind speed (m/s)")
ax6.text(1.5, 800, "f")

ax1.set_yticks(
    [1000, 600, 400, 300, 200, 100],
    labels=["1000", "600", "400", "300", "200", "100"],
)


# iesc1 = escPwrf.get_timeind(escSnd1.dtime)
# iesc2 = escPwrf.get_timeind(escSnd2.dtime)

# Get difference statistics: max, min, mean, rms difference
esc1_dtemp, esc1_drh, esc1_dwind = escSnd1.get_diffs(escPwrf1)
esc2_dtemp, esc2_drh, esc2_dwind = escSnd2.get_diffs(escPwrf2)


print()
print(esc_snd_file1)
print()
print(f'Max: {max(abs(esc1_drh["max"]),abs(esc1_drh["min"]))}')
print(f'Mean: {esc1_drh["mean"]}')
print(f'RMS: {esc1_drh["rms"]}')


print()
print(esc_snd_file2)
print()
print(f'Max: {max(abs(esc2_drh["max"]),abs(esc2_drh["min"]))}')
print(f'Mean: {esc2_drh["mean"]}')
print(f'RMS: {esc2_drh["rms"]}')

# print()
# print(esc_snd_file1)
# print()
# for x in esc1_dtemp:
#     print(x)
# for x in esc1_drh:
#     print(x)
# for x in esc1_dwind:
#     print(x)

# print()
# print(esc_snd_file2)
# print()
# for x in esc2_dtemp:
#     print(x)
# for x in esc2_drh:
#     print(x)
# for x in esc2_dwind:
#     print(x)
