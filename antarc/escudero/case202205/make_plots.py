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
    plt.figure(num=2, figsize=[6.5, 5.6])
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
    ax1.set_ylim([-2, 120])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax1.tick_params(axis="x", rotation=rot, labelcolor="white")
    ax1.set_xlim(start_date, final_date)
    ax1.set_ylabel("SWD (W/m$^{2}$)")
    ax1.legend()

    ax2.plot(lwd_date, lwd, color="blue", label="Measured")
    ax2.plot(lwd_clr_date, lwd_clr, "o:", color="blue", label="Clear")
    ax2.set_ylim([150, 350])
    ax2.legend()
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax2.tick_params(axis="x", rotation=rot, labelcolor="white")
    ax2.set_xlabel("Day and Hour in 2022")
    ax2.set_ylabel("LWD (W/m$^{2}$)")

    tot_force = swd_force + lwd_force
    ax3.plot(lwd_date, lwd_force, color="blue", linewidth=2, label="LWD")
    ax3.plot(swd_date, swd_force, color="red", label="SWD")
    ax3.plot(swd_date, tot_force, color="orange", label="Total")
    ax3.plot(lwd_date, np.zeros(len(lwd_date)), "k:")
    ax3.legend()
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax3.set_ylim([-100, 150])
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_ylabel("Forcing (W/m$^{2}$)")
    ax3.legend(loc=[0.01, 0.62])

    if savef:
        plt.savefig(fname)


def plot_forcings(
    swd_date,
    swd,
    swd_clr,
    swd_force,
    lwd_date,
    lwd,
    lwd_clr_date,
    lwd_clr,
    lwd_force,
    start_date,
    final_date,
    datefmt: str = "%d",
    xlab: str = "Day of 2022/05",
    savef: bool = False,
    fname: str = "",
):

    rot = 0
    plt.figure(num=3)  # , figsize=[6.5, 5.6])
    plt.clf()
    ax3 = plt.subplot(111)
    ax3.set_position([0.13, 0.11, 0.84, 0.84])

    tot_force = swd_force + lwd_force
    ax3.plot(lwd_date, lwd_force, color="orange", linewidth=3, label="LWD")
    ax3.plot(swd_date, swd_force, color="red", label="SWD")
    ax3.plot(swd_date, tot_force, color="blue", label="Total")
    ax3.plot(lwd_date, np.zeros(len(lwd_date)), "k:")
    ax3.legend()
    ax3.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_xlabel(xlab)

    ax3.set_xlim(start_date, final_date)
    ax3.set_ylim([-75, 150])
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_ylabel("Forcing (W/m$^{2}$)")

    if savef:
        plt.savefig(fname)
