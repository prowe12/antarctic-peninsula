#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 11:32:46 2023

@author: prowe
"""

from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import numpy as np


def plot_measured_and_clear(
    swd_date, swd, swd_clr, lwd_date, lwd, lwd_clr_date, lwd_clr, savef, fname
):
    rot = 0
    datefmt = "%d"
    plt.figure(num=1)
    ax1 = plt.subplot(211)
    ax3 = plt.subplot(212, sharex=ax1)
    wid = 0.85
    hit = 0.4
    ax1.set_position([0.12, 0.57, wid, hit])
    ax3.set_position([0.12, 0.10, wid, hit])

    ax1.plot(swd_date, swd, "r", label="Measured")
    ax1.plot(swd_date, swd_clr, "k--", label="Clear")
    # global down
    # ax1.plot(swd_date[ind], swd_clear_lib[:, 0], "ks")  # direct
    # ax1.plot(swd_date[ind], swd_clear_lib[:, 2], "c.")  # diffuse down
    # ax1.plot(swd_date[ind], swd_clear_lib[:, 3], "m+")  # diffuse up
    ax1.legend()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
    ax1.tick_params(axis="x", rotation=rot)
    ax1.set_ylabel("SWD (W/m$^{2}$)")

    ax3.plot(lwd_date, lwd, color="blue", label="Measured")
    ax3.plot(lwd_clr_date, lwd_clr, "o:", color="blue", label="Clear")
    ax3.legend()
    ax3.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_xlabel("Day of 2022/05")
    ax3.set_ylabel("LWD (W/m$^{2}$)")
    ax3.legend(loc="upper left")

    if savef:
        plt.savefig(fname)


def plot_meas_clear_forcings(
    swd_date,
    swd,
    swd_clr,
    swd_force,
    lwd_date,
    lwd,
    lwd_clr_date,
    lwd_clr,
    lwd_force,
    savef,
    fname,
    start_date,
    final_date,
):

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
    ax1.plot(swd_date, swd_clr, "k--", label="Clear")
    ax1.legend()
    ax1.set_ylim([-2, 1000])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax1.tick_params(axis="x", rotation=rot, labelcolor="white")
    ax1.set_xlim(start_date, final_date)
    ax1.set_ylabel("SWD (W/m$^{2}$)")
    ax1.legend()

    ax2.plot(lwd_date, lwd, color="blue", label="Measured")
    ax2.plot(lwd_clr_date, lwd_clr, "o:", color="blue", label="Clear")
    ax2.set_ylim([0, 400])
    ax2.legend()
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax2.tick_params(axis="x", rotation=rot, labelcolor="white")
    ax2.set_xlabel("Day and Hour in 2022")
    ax2.set_ylabel("LWD (W/m$^{2}$)")

    tot_force = swd_force + lwd_force
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

    # fname = out_dir + "lwd_swd_forcing_short_timespan.png"
    if savef:
        plt.savefig(fname)
