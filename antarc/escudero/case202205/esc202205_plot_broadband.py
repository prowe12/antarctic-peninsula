#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:21:08 2022

@author: prowe
"""

# Built-in modules
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

# My modules
from antarc.plot_broadband import plot_broadband
from antarc.escudero.parameters import esc202205 as esc_case
from antarc.getfilenames import getfilenames_daterange as getfnames


def plot_broadband_esc(savefigs: bool, outdir: str):
    """Plot broadband data"""
    # TODO: remove code that is redundant here and above
    # TODO: Shade in foehn region here (start with commented-out code)

    figname1 = outdir + "esc_broadband.png"

    date1 = esc_case.DATE1 - dt.timedelta(days=2)
    date2 = esc_case.DATE2

    # Set up the figure
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(7, 5), sharex=True, num=1)
    wid = 0.84
    height = 0.43
    lft = 0.11
    bot2 = 0.08
    bot1 = bot2 + height + 0.03
    ax1.set_position([lft, bot1, wid, height])
    ax2.set_position([lft, bot2, wid, height])

    # # # # # # #   BROADBAND RADIATION   # # # # # # # # # # #
    # Add the broadband radiation
    fnamelen = len("esc_swd20220206_0301.csv")

    fnames_swd = getfnames(dir_swd, fnamelen, "esc_swd", date1, date2)
    fnames_lwd = getfnames(dir_lwd, fnamelen, "esc_lwd", date1, date2)

    plot_broadband(fnames_swd, fnames_lwd, ax1, ax2)
    ax1.autoscale(enable=True, axis="x", tight=True)
    ax1.autoscale(enable=True, axis="x", tight=True)

    # Format the dates
    date_form = DateFormatter("%m/%d")
    ax2.xaxis.set_major_formatter(date_form)

    if savefigs:
        plt.savefig(figname1)

    # Same plot, but with ERA5 interpolated to location
    # figname2 = outdir + "esc_broadband_with_era5.png"
    # direc = PROJ_DIR + "/NSF_AP/case_studies/Feb_2022/ERA5_at_stations/"
    # filename = "Escudero_SWD_LWD_Feb_06_10_era5.nc"
    # era_date, era_swd, era_lwd = read_era5_broadband_down(direc + filename)
    # ax1.plot(era_date, era_swd, ".", color="orange", label="ERA5")
    # ax1.legend()
    # ax2.plot(era_date, era_lwd, ".", color="orange", label="ERA5")
    # ax2.legend()
    # if savefigs:
    #     plt.savefig(figname2)


def one_plot_broadband_esc(savefigs: bool, outdir: str):
    """Plot broadband data on a single panel"""
    # TODO: do we want to keep this?
    figname1 = outdir + "esc_broadband.png"

    date1 = esc_case.DATE1 - dt.timedelta(days=2)
    date2 = esc_case.DATE2

    # Set up the figure
    fig, (ax1) = plt.subplots(nrows=1, figsize=(7, 5), num=1)
    wid = 0.84
    height = 0.43
    lft = 0.11
    bot2 = 0.08
    bot1 = bot2 + height + 0.03
    ax1.set_position([lft, bot1, wid, height])
    # ax2.set_position([lft, bot2, wid, height])

    # # # # # # #   BROADBAND RADIATION   # # # # # # # # # # #
    # Add the broadband radiation
    fnamelen = len("esc_swd20220206_0301.csv")

    fnames_swd = getfnames(dir_swd, fnamelen, "esc_swd", date1, date2)
    fnames_lwd = getfnames(dir_lwd, fnamelen, "esc_lwd", date1, date2)

    plot_broadband(fnames_swd, fnames_lwd, ax1, ax1)
    ax1.autoscale(enable=True, axis="x", tight=True)
    ax1.autoscale(enable=True, axis="x", tight=True)

    # Format the dates
    date_form = DateFormatter("%m/%d")
    ax1.xaxis.set_major_formatter(date_form)

    if savefigs:
        plt.savefig(figname1)


dir_swd = "/Users/prowe/sync/measurements/Escudero/pyranometer/stand/2022/"
dir_lwd = "/Users/prowe/sync/measurements/Escudero/pyrgeometer/stand/2022/"

outdir = "antarc/escudero/case202205/"
# plot_broadband_with_solar_elevation_esc(False, outdir)
# plot_broadband_esc(True, outdir)
one_plot_broadband_esc(False, outdir)
