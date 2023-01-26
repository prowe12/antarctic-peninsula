#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:21:08 2022

@author: prowe
"""

from antarc.plot_broadband import plot_broadband

import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter


from antarc.plot_sun_elev import PlotSunElevation
from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import esc202202
from antarc.getfilenames import getfilenames_daterange
from antarc.era5_read_station_broadband import read_era5_broadband_down


def plot_broadband_with_solar_elevation_esc(savefigs: bool, outdir: str):
    """Plot that includes solar elevation"""

    figname1 = outdir + "esc_broadband_with_sun_elev.png"
    figname2 = outdir + "esc_broadband_with_sun_elev_era5.png"

    # Location
    lat = esc_params.LATITUDE
    lon = esc_params.LONGITUDE
    alt = esc_params.ALTITUDE

    dtime = esc202202.DTIME
    date1 = esc202202.DATE1
    date2 = esc202202.DATE2
    dtfoehn1 = esc202202.DT_FOEHN1
    dtfoehn2 = esc202202.DT_FOEHN2
    dt_foehnmax = esc202202.DT_FOEHNMAX

    # Set up the figure
    fig, (ax1, ax2, ax3) = plt.subplots(
        nrows=3, figsize=(8, 8), sharex=True, num=1
    )
    wid = 0.88
    hght = 0.28
    lft = 0.09
    ax1.set_position([lft, 0.68, wid, hght])
    ax2.set_position([lft, 0.37, wid, hght])
    ax3.set_position([lft, 0.06, wid, hght])

    # # # # # # #   SUN   # # # # # # # # # # #
    # Get the class instance
    plotSunElevation = PlotSunElevation(lat, lon, alt, dtime)
    ax = ax1
    dtimes, elev = plotSunElevation.daterange(date1, date2, ax)

    # Put a star at the peak warming
    tstamps = [x.timestamp() for x in dtimes]
    elev_foehn_max = np.interp(dt_foehnmax.timestamp(), tstamps, elev)
    ax.plot(dt_foehnmax, elev_foehn_max, "s", color="red", label="Foehn max")

    # Shade the foehn warming period
    i = 0
    dtime = dtimes[i]
    while dtime < dtfoehn1:
        i += 1
        dtime = dtimes[i]
    i1 = i

    while dtime < dtfoehn2:
        i += 1
        dtime = dtimes[i]
    i2 = i
    ax.axvspan(
        dtimes[i1],
        dtimes[i2],
        color="gray",
        alpha=0.5,
        label="strong foehn warming",
    )
    ax.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(ax.xaxis.get_major_locator())
    )
    ax.legend(loc="upper right")
    ax.set_title("Escudero Station")
    ax.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(ax.xaxis.get_major_locator())
    )

    # # # # # # #   BROADBAND RADIATION   # # # # # # # # # # #
    # Add the broadband radiation
    dir_swd = "/Users/prowe/sync/measurements/Escudero/pyranometer/stand/2022/"
    dir_lwd = "/Users/prowe/sync/measurements/Escudero/pyrgeometer/stand/2022/"
    fnamelen = len("esc_swd20220206_0301.csv")

    fnames_swd = getfilenames_daterange(
        dir_swd, fnamelen, "esc_swd", date1, date2
    )
    fnames_lwd = getfilenames_daterange(
        dir_lwd, fnamelen, "esc_lwd", date1, date2
    )

    plot_broadband(fnames_swd, fnames_lwd, ax2, ax3)

    if savefigs:
        plt.savefig(figname1)

    # Same plot, but with ERA5 interpolated to location
    direc = "/Users/prowe/Sync/projects/NSF_AP/case_studies/Feb_2022/ERA5_at_stations/"
    filename = "Escudero_SWD_LWD_Feb_06_10_era5.nc"

    era_date, era_swd, era_lwd = read_era5_broadband_down(direc + filename)

    # ax1.plot(era_date, era_swd * 0.06, ".-", color="green", label="ERA5")
    ax2.plot(era_date, era_swd, ".", color="orange", label="ERA5")
    ax2.legend()

    ax3.plot(era_date, era_lwd, ".", color="orange", label="ERA5")
    ax3.legend()

    if savefigs:
        plt.savefig(figname2)


def plot_broadband_esc(savefigs: bool, outdir: str):
    """Plot broadband data"""
    # TODO: remove code that is redundant here and above
    # TODO: Shade in foehn region here (start with commented-out code)

    figname1 = outdir + "esc_broadband.png"
    figname2 = outdir + "esc_broadband_with_era5.png"

    # dtime = esc202202.DTIME
    date1 = esc202202.DATE1 - dt.timedelta(days=2)
    date2 = esc202202.DATE2
    # dtfoehn1 = esc202202.DT_FOEHN1
    # dtfoehn2 = esc202202.DT_FOEHN2
    # dt_foehnmax = esc202202.DT_FOEHNMAX

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
    dir_swd = "/Users/prowe/sync/measurements/Escudero/pyranometer/stand/2022/"
    dir_lwd = "/Users/prowe/sync/measurements/Escudero/pyrgeometer/stand/2022/"
    fnamelen = len("esc_swd20220206_0301.csv")

    fnames_swd = getfilenames_daterange(
        dir_swd, fnamelen, "esc_swd", date1, date2
    )
    fnames_lwd = getfilenames_daterange(
        dir_lwd, fnamelen, "esc_lwd", date1, date2
    )

    plot_broadband(fnames_swd, fnames_lwd, ax1, ax2)
    ax1.autoscale(enable=True, axis="x", tight=True)
    ax1.autoscale(enable=True, axis="x", tight=True)

    # Format the dates
    date_form = DateFormatter("%m/%d")
    # ax1.xaxis.set_major_formatter(date_form)
    ax2.xaxis.set_major_formatter(date_form)

    if savefigs:
        plt.savefig(figname1)

    # Put a star at the peak warming
    # tstamps = [x.timestamp() for x in dtimes]
    # elev_foehn_max = np.interp(dt_foehnmax.timestamp(), tstamps, elev)
    # ax.plot(dt_foehnmax, elev_foehn_max, "s", color="red", label="Foehn max")

    # # Shade the foehn warming period
    # i = 0
    # dtime = dtimes[i]
    # while dtime < dtfoehn1:
    #     i += 1
    #     dtime = dtimes[i]
    # i1 = i

    # while dtime < dtfoehn2:
    #     i += 1
    #     dtime = dtimes[i]
    # i2 = i
    # ax.axvspan(
    #     dtimes[i1],
    #     dtimes[i2],
    #     color="gray",
    #     alpha=0.5,
    #     label="strong foehn warming",
    # )
    # ax.xaxis.set_major_formatter(
    #     mdates.ConciseDateFormatter(ax.xaxis.get_major_locator())
    # )
    # ax.legend(loc="upper right")
    # ax.set_title("Escudero Station")
    # ax.xaxis.set_major_formatter(
    #     mdates.ConciseDateFormatter(ax.xaxis.get_major_locator())
    # )

    # Same plot, but with ERA5 interpolated to location

    direc = "/Users/prowe/Sync/projects/NSF_AP/case_studies/Feb_2022/ERA5_at_stations/"
    filename = "Escudero_SWD_LWD_Feb_06_10_era5.nc"

    era_date, era_swd, era_lwd = read_era5_broadband_down(direc + filename)

    # ax1.plot(era_date, era_swd * 0.06, ".-", color="green", label="ERA5")
    ax1.plot(era_date, era_swd, ".", color="orange", label="ERA5")
    ax1.legend()

    ax2.plot(era_date, era_lwd, ".", color="orange", label="ERA5")
    ax2.legend()

    if savefigs:
        plt.savefig(figname2)


outdir = "antarc/escudero/case202202/"
# plot_broadband_with_solar_elevation_esc(False, outdir)
plot_broadband_esc(True, outdir)


# EXTRA STUFF
#
# i1 = [i for i,x in enumerate(dtimes) if x > dfoehn1]
# i2 = [i for i,x in enumerate(dtimes) if x < dfoehn2]
# ifoehn = np.intersect1d(i1, i2)
# plt.fill_between(np.array(dtimes)[ifoehn], elev[ifoehn], color='gray', alpha=.1)


# elev = plotSunElevation.get_solar_elevation()

# Create a plot for the day
# figno, hour, elev = plotSunElevation.diurnal()
# plotSunElevation.annotate_diurnal(figno, hour, elev)

# # Extend the plot backwards in time
# plotSunElevation.extend(figno, date1, date2)


# plt.axvspan(0, 21)


# plt.figure(figno)

# Create plot for the month
# plotSunElevation.month()
