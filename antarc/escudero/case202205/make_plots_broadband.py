#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 08:56:46 2022

@author: prowe

Purpose: Plot broadband radiation at Escudero on 2022/02
"""

# Dependencies
from matplotlib import pyplot as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import pandas as pd
from os.path import exists

# My modules
from antarc.escudero.parameters import esc202205 as esc_case
from antarc.load_rad_flux import load_rad_flux
from antarc.escudero.case202205.get_swd_flux import get_obs_swd
from antarc.escudero.case202205.get_swd_flux import get_pysolar_swd
from antarc.escudero.case202205.get_swd_flux import get_libradtran
from antarc.escudero.case202205.get_lwd_flux import get_lwd, get_lwd_clear
from antarc.escudero.parameters import swd_params


def plot_measured_and_clear():
    rot = 0
    datefmt = "%d"
    plt.figure(num=1)
    ax1 = plt.subplot(211)
    ax3 = plt.subplot(212, sharex=ax1)
    wid = 0.85
    hit = 0.4
    ax1.set_position([0.12, 0.57, wid, hit])
    ax3.set_position([0.12, 0.10, wid, hit])

    ax1.plot(swd_date_long, swd_long, "r", label="Measured")
    ax1.plot(swd_date_long, swd_clear_lib_long, "k--", label="Clear")
    # global down
    # ax1.plot(swd_date[ind], swd_clear_lib[:, 0], "ks")  # direct
    # ax1.plot(swd_date[ind], swd_clear_lib[:, 2], "c.")  # diffuse down
    # ax1.plot(swd_date[ind], swd_clear_lib[:, 3], "m+")  # diffuse up
    ax1.legend()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
    ax1.tick_params(axis="x", rotation=rot)
    ax1.set_ylabel("SWD (W/m$^{2}$)")

    ax3.plot(lwd_date, lwd, color="blue", label="Measured")
    ax3.plot(lwd_clear_date, lwd_clear, "o:", color="blue", label="Clear")
    ax3.legend()
    ax3.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_xlabel("Day of 2022/05")
    ax3.set_ylabel("LWD (W/m$^{2}$)")
    ax3.legend(loc="upper left")

    if savefigs:
        plt.savefig(fig1name)


def plot_meas_clear_forcings():
    rot = 40
    # datefmt = "%d %H"
    # fdate = swd_date[isw - 1]
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    plt.figure(num=3, figsize=[6.5, 5.6])
    plt.clf()
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312, sharex=ax1)
    ax3 = plt.subplot(313, sharex=ax1)
    wid = 0.84
    ax1.set_position([0.13, 0.67, wid, 0.30])
    ax2.set_position([0.13, 0.45, wid, 0.18])
    ax3.set_position([0.13, 0.09, wid, 0.33])

    ax1.plot(swd_date, swd, color="red", label="Meas")
    ax1.plot(swd_date, swd_clear_lib_interp, "k--", label="Clear")
    ax1.legend()
    ax1.set_ylim([-2, 1000])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax1.tick_params(axis="x", rotation=rot, labelcolor="white")
    ax1.set_xlim(start_date, final_date)
    ax1.set_ylabel("SWD (W/m$^{2}$)")
    ax1.legend()

    ax2.plot(lwd_date, lwd, color="blue", label="Measured")
    ax2.plot(lwd_clear_date, lwd_clear, "o:", color="blue", label="Clear")
    ax2.set_ylim([0, 400])
    ax2.legend()
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax2.tick_params(axis="x", rotation=rot, labelcolor="white")
    ax2.set_xlabel("Day and Hour in 2022")
    ax2.set_ylabel("LWD (W/m$^{2}$)")

    tot_force = swd_force + lwd_force_sw_date
    ax3.plot(lwd_date, lwd_force, color="blue", label="LWD")
    ax3.plot(swd_date, swd_force, color="red", label="SWD")
    ax3.plot(swd_date, tot_force, "k", label="Total")
    ax3.plot(lwd_date, np.zeros(len(lwd_date)), "k:")
    ax3.legend()
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax3.set_ylim([-600, 600])
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_ylabel("Forcing (W/m$^{2}$)")
    ax3.legend(loc=[0.01, 0.62])

    if savefigs:
        plt.savefig(out_dir + "lwd_swd_forcing_short_timespan.png")


pyr_dir = swd_params.STAND_DIR
pyr_fmt = swd_params.PYR_FILEFORMAT

out_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/May_2022/figures/"
fig1name = out_dir + "lwd_swd.png"
fig2name = out_dir + "lwd_swd_zoom.png"

utc = dt.timezone.utc

# Run flags
get_libradtran_and_savetofile = False  # Slow step
savefigs = False

date1 = esc_case.DATE1  # dt.datetime(2022, 2, 5, tzinfo=utc)
date2 = esc_case.DATE2  # dt.datetime(2022, 2, 11, tzinfo=utc)

# SWD
swd_date, swd = get_obs_swd()
swd_pysolar, diffuse, diffuse_fac = get_pysolar_swd(swd_date)


# swd_date, swd, swd_clear_lib, swd_clear, diffuse, diffuse_fac = get_swd()
pyr_dir += date1.strftime("%Y") + "/"
swd_date_long, swd_long = load_rad_flux(pyr_dir, pyr_fmt, date1, date2)

# If needed, calculate the libradtran results and save to file
outfile = out_dir + "libradtran.csv"
if get_libradtran_and_savetofile:
    libdate, clear_lib = get_libradtran(swd_date_long, outfile)
    swd_clear_lib = clear_lib[:, 0]
elif exists(outfile):
    # Load the libradtran results, which have columns: Date,swd,o2,o3,o4,o5,o6
    clear_lib = pd.read_csv(outfile, delimiter=",\s")
    libdate = pd.to_datetime(clear_lib["Date"])
    swd_clear_lib = clear_lib["swd"]
else:
    msg1 = "Change get_libradtran_and_savetofile to True to calculate and "
    msg2 = "save clear-sky sw down from libradtran"
    raise ValueError(msg1 + msg2)

lib_tstamp = [x.timestamp() for x in libdate]
swd_long_tstamp = [x.timestamp() for x in swd_date_long]
swd_clear_lib_long = np.interp(swd_long_tstamp, lib_tstamp, swd_clear_lib)

start_date = swd_date[0]
final_date = swd_date[-1]

# LWD
lwd_clear_date, lwd_clear = get_lwd_clear()
lwd_date, lwd = get_lwd(date1, date2)
lwd_tstamp = [x.timestamp() for x in lwd_date]
lwd_clear_tstamp = [x.timestamp() for x in lwd_clear_date]
lwd_clear = np.array(lwd_clear)

# Get LWD forcing
lwd_clear_interp = np.interp(lwd_tstamp, lwd_clear_tstamp, lwd_clear)

# Cumulative forcing for LWD
# dt.datetime(fdate.year, fdate.month, fdate.day, fdate.hour + 1, 0)
i1 = np.where(np.array(lwd_tstamp) >= start_date.timestamp())[0][0]
i2 = np.where(np.array(lwd_tstamp) <= final_date.timestamp())[0][-1] + 1
lwd_force = lwd - lwd_clear_interp
ntimes = i2 - i1
lwf_tot = np.zeros([ntimes])
lwf_tot_mwh = np.zeros([ntimes])
count = 0
for i in range(i1, i2):
    lwf = lwd_force[i] * (lwd_tstamp[i] - lwd_tstamp[i - 1])
    lwf_mwh = lwd_force[i] * (lwd_tstamp[i] - lwd_tstamp[i - 1]) / 3600 / 1000
    lwf_tot[count] = lwf_tot[count - 1] + lwf
    lwf_tot_mwh[count] = lwf_tot_mwh[count - 1] + lwf_mwh
    count += 1
il1 = i1
il2 = i2


# Cumulative forcing for SWD
swd_tstamp = [x.timestamp() for x in swd_date]
swd_clear_lib_interp = np.interp(swd_tstamp, lib_tstamp, swd_clear_lib)
i1 = np.where(np.array(swd_tstamp) >= start_date.timestamp())[0][0]
i2 = np.where(np.array(swd_tstamp) <= final_date.timestamp())[0][-1]
swd_force = swd - swd_clear_lib_interp
ntimes = i2 - i1
swf_tot = np.zeros([ntimes])
swf_tot_mwh = np.zeros([ntimes])
count = 0
for i in range(i1, i2):
    swf = swd_force[i] * (swd_tstamp[i] - swd_tstamp[i - 1])
    swf_mwh = swd_force[i] * (swd_tstamp[i] - swd_tstamp[i - 1]) / 3600 / 1000
    swf_tot[count] = swf_tot[count - 1] + swf
    swf_tot_mwh[count] = swf_tot_mwh[count - 1] + swf_mwh
    count += 1
is1 = i1
is2 = i2


# lwd_tot = np.interp(lwd_tstamp, lwd_clear_interp)
# swd_tot = np.interp(lwd_tstamp, lwd_clear_interp)
lwf_tot_sw_date = np.interp(swd_tstamp[is1:is2], lwd_tstamp[il1:il2], lwf_tot)
lwd_force_sw_date = np.interp(swd_tstamp, lwd_tstamp, lwd_force)

# ax1 = plt.subplot(221)
# ax2 = plt.subplot(222, sharex=ax1)
# ax3 = plt.subplot(223, sharex=ax1)
# ax4 = plt.subplot(224, sharex=ax1)

# ax2.plot(lwd_date[:ilw], lwd_force[:ilw], color="blue", label="LWD Forcing")
# ax2.plot(swd_date[:isw], swd_force[:isw], color="red", label="SWD Forcing")
# ax2.legend()
# ax2.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# ax2.tick_params(axis="x", rotation=rot)

# ax4.plot(lwd_date[:ilw], lwf_tot[:ilw] / 1e6, color="blue", label="LWD")
# ax4.plot(swd_date[:isw], swf_tot[:isw] / 1e6, color="red", label="SWD")
# ax4.plot(
#     swd_date[:isw],
#     (swf_tot[:isw] + lwf_tot_sw_date[:isw]) / 1e6,
#     "k",
#     label="Total",
# )
# ax4.legend()
# ax4.set_ylabel("Cumulative Forcing (MJ/m$^{2}$)")
# ax4.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# ax4.set_xlabel("Day of 2022/02")
# ax4.tick_params(axis="x", rotation=rot)

totf = swf_tot + lwf_tot_sw_date
totfstr = str(round(totf[-1] / 1e6))

plot_measured_and_clear()
plot_meas_clear_forcings_cumulative_shorttime()
plot_meas_clear_forcings()

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
