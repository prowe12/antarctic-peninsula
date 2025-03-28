#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 08:56:46 2022

@author: prowe

Purpose: Do analysis and make plots for cloud forcing

"""

# Dependencies
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

# My modules
from antarc import params
from antarc.escudero.all.results import get_broadband_monthly_era5
from antarc.escudero.all.results import get_broadband_monthly_meas_only

# Parameter modules
from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import swd_params, lwd_params
from antarc.escudero.all.results.params import SAVEFIGS


get_broadband_monthly_meas = (
    get_broadband_monthly_meas_only.get_broadband_monthly_meas
)
get_monthly_era5 = get_broadband_monthly_era5.get_broadband_monthly_era5


# Run flags
figdir = "antarc/Escudero/all/figures/"
years = range(2017, 2024)
pts_min = 200  # min hours in month to compute stats (of 672 to 744)


# Directories and variables needed from other param files
lwd_base_dir = lwd_params.LWD_DIR
swd_base_dir = swd_params.SWD_DIR
lwd_fmt = lwd_params.LWD_FILEFORMAT
swd_fmt = swd_params.SWD_FILEFORMAT
era_dir = f"{params.MEAS_DIR}Escudero/era5/broadband/"
era5fmt = "era5_esc_broadband_*[0-9]*.nc"
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE


monthly = get_broadband_monthly_meas(
    years,
    lwd_base_dir,
    swd_base_dir,
    lwd_fmt,
    swd_fmt,
    lat,
    lon,
    pts_min,
)

era5 = get_monthly_era5(
    years,
    era_dir,
    era5fmt,
    lat,
    lon,
)


# Setup for figures
mcolor = "blue"  # "darkturquoise"  # measurement color
colors = ["magenta", "purple", "blue", "cyan", "green", "orange", "red"]
colors = ["red", "orange", "green", "cyan", "blue", "purple", "magenta"]

xleg = dt.datetime(2020, 7, 20)
xleg3 = dt.datetime(2020, 4, 20)
xmin = dt.datetime(2020, 1, 1)
xmax = dt.datetime(2020, 12, 1)


# 2-panel plot showing number of measurements

fig, (ax2a, ax2b) = plt.subplots(2, 1, num=4, clear=True)
ax2a.set_position([0.11, 0.57, 0.73, 0.39])
ax2b.set_position([0.11, 0.1, 0.73, 0.39])

mon = np.arange(1, 13)
swd_bot = np.zeros(12)  # [0 for i in range(12)]
lwd_bot = np.zeros(12)  # 0 for i in range(12)]
for i, (year, color) in enumerate(zip(years[::-1], colors[::-1])):
    label = str(year)
    swd_yr = monthly["swd"]["npts"][-(i + 1)]
    lwd_yr = monthly["lwd"]["npts"][-(i + 1)]

    swd_yr[np.isnan(swd_yr)] = 0
    lwd_yr[np.isnan(lwd_yr)] = 0

    ax2a.bar(mon, swd_yr, bottom=swd_bot, color=color, label=label)
    ax2b.bar(mon, lwd_yr, bottom=lwd_bot, color=color)

    ind = np.where(monthly["swd"]["npts"][-(i + 1)] < pts_min)[0]
    ind = np.intersect1d(
        ind, np.where(monthly["swd"]["npts"][-(i + 1)] > 0)[0]
    )
    if np.any(ind):
        ax2a.plot(
            mon[ind],
            swd_bot[ind] + swd_yr[ind] / 2,
            "^",
            color=color,
        )

    ind = np.where(monthly["lwd"]["npts"][-(i + 1)] < pts_min)[0]
    ind = np.intersect1d(
        ind, np.where(monthly["lwd"]["npts"][-(i + 1)] > 0)[0]
    )
    if np.any(ind):
        ax2b.plot(
            mon[ind],
            lwd_bot[ind] + lwd_yr[ind] / 2,
            "^",
            color=color,
        )

    swd_bot += swd_yr
    lwd_bot += lwd_yr


# fmt: off
handles, labels = ax2a.get_legend_handles_labels()
ax2a.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.01, 1.0), loc="upper left")
ax2a.set_ylabel("DSW Flux (Number)")
ax2a.set_xlim([0.5, 12.5])

ax2b.set_ylabel("DLW Flux (Number)")
ax2b.set_xlim([0.5, 12.5])
ax2b.set_xlabel("Month")
# fmt: on

if SAVEFIGS:
    figname = figdir + "broadband_meas_v_era5_npts"
    plt.savefig(figname + ".png", format="png")
    plt.savefig(figname + ".eps", format="eps")


colors = ["red", "orange", "green", "cyan", "blue", "purple", "magenta"]

# 2-panel 1 column plot of monthly mean DLW fluxes
fig, axs = plt.subplots(2, 1, num=10, clear=True)  # , figsize=[5, 8])
ax1a = axs[0]
ax1b = axs[1]

dtm = np.array([dt.datetime(2020, x, 1) for x in range(1, 13)])
off = dt.timedelta(days=3)
wid = 0.87
lft = 0.11
height1 = 0.47
height2 = 0.37
ax1a.set_position([lft, 0.5, wid, height1])
ax1b.set_position([lft, 0.1, wid, height2])


# Require at least a certain number of points
for i, key in enumerate(["swd", "lwd"]):
    monthly["swd"]["mean"][monthly["swd"]["npts"] < pts_min] = np.nan


# Plot the monthly averages over all years
for i, key in enumerate(["swd", "lwd"]):
    ax = axs[i]

    meanval = np.nanmean(monthly[key]["mean"], axis=0)
    # minval = np.nanmax(monthly[key]["min"], axis=0)
    # maxval = (np.nanmin(monthly[key]["max"], axis=0),)
    # ax.errorbar(dtm, meanval, yerr=np.vstack([minval, maxval]), fmt="none")
    ax.plot(dtm, meanval, color=mcolor, linewidth=4, label="Measured")

    meanval = np.nanmean(era5[key + "_clr"]["era5_all"], axis=0)
    ax.plot(dtm, meanval, color="c", linewidth=4, label="ERA5, clear")

i = 0
for year, color in zip(years, colors):
    label = str(year)
    ax1a.plot(dtm, monthly["swd"]["mean"][i], "+", color=color, label=label)
    ax1b.plot(dtm, monthly["lwd"]["mean"][i], "+", color=color)
    i += 1

ax1a.set_ylabel("DSW Flux (W/m$^{2}$)")
ax1a.legend(bbox_to_anchor=(0.31, 1), loc="upper left")
ax1a.xaxis.set_major_formatter(DateFormatter("%m"))
ax1a.text(xmax - dt.timedelta(days=15), 20, "a)")

ax1b.set_ylabel("DLW Flux (W/m$^{2}$)")
ax1b.xaxis.set_major_formatter(DateFormatter("%m"))
ax1b.text(xmax - dt.timedelta(days=15), 245, "b)")

ax1a.set_xticklabels([])
ax1b.set_xlabel("Month")

# Format the dates
date_form = DateFormatter("%m")
ax2b.xaxis.set_major_formatter(date_form)

xrange = [xmin - 2.5 * off, xmax + 2.5 * off]
for ax in axs:
    ax.set_xlim(xrange)

if SAVEFIGS:
    figname = figdir + "broadband_meas_bymonth"
    plt.savefig(figname + ".png", format="png", dpi=600)
    plt.savefig(figname + ".eps", format="eps")
# fmt: on


# fmt: off
dtm = np.array([dt.datetime(2020, x, 1) for x in range(1, 13)])
off = dt.timedelta(days=3)

wid = 0.87
lft = 0.11
height1 = .27
height2 = .2
height3 = .3

# 2-panel 1 column plot of monthly mean flux
fig, axs = plt.subplots(3, 1, num=10, clear=True)  # , figsize=[5, 8])
ax1a = axs[0]
ax1b = axs[1]
ax1c = axs[2]
ax1a.set_position([lft, 0.7, wid, height1])
ax1b.set_position([lft, 0.45, wid, height2])
ax1c.set_position([lft, 0.1, wid, height3])


# Require at least a certain number of points
for i, key in enumerate(["swd", "lwd"]):
    monthly["swd"]["mean"][monthly["swd"]["npts"] < pts_min] = np.nan


# Plot the monthly averages over all years
for i, key in enumerate(["swd", "lwd"]):
    ax = axs[i]
    meanval = np.nanmean(era5[key + "_clr"]["era5_all"], axis=0)
    ax.plot(dtm, meanval, color="c", label="ERA5, clear")
    meanval = np.nanmean(monthly[key]["mean"], axis=0)
    # minval = np.nanmax(monthly[key]["min"], axis=0)
    # maxval = (np.nanmin(monthly[key]["max"], axis=0),)
    # ax.errorbar(dtm, meanval, yerr=np.vstack([minval, maxval]), fmt="none")
    ax.plot(dtm, meanval, color=mcolor, label="Measured")
ax = axs[2]
ax.plot(dtm, np.nanmean(monthly['lwd']["mean"] + monthly['swd']["mean"], axis=0), color=mcolor)
ax.plot(dtm, np.nanmean(era5['lwd_clr']['era5_all'] + era5['swd_clr']["era5_all"], axis=0),'c')



i = 0
for year, color in zip(years, colors):
    label = str(year)
    ax1a.plot(dtm, monthly["swd"]["mean"][i], "+", color=color, label=label)
    ax1b.plot(dtm, monthly["lwd"]["mean"][i], "+", color=color)
    ax1c.plot(dtm, monthly["lwd"]["mean"][i] + monthly["swd"]["mean"][i], "+", color=color)
    i += 1

ax1a.legend(bbox_to_anchor=(0.18, 1), loc="upper left", ncol=3, columnspacing=0.8)
ax1a.xaxis.set_major_formatter(DateFormatter("%m"))

ax1b.xaxis.set_major_formatter(DateFormatter("%m"))

ax1a.set_ylabel("Flux (W/m$^{2}$)")
ax1b.set_ylabel("Flux (W/m$^{2}$)")
ax1c.set_ylabel("Flux (W/m$^{2}$)")

ax1a.text(xmax - dt.timedelta(days=40), 50, "a) DSW")
ax1b.text(xmax - dt.timedelta(days=40), 320, "b) DLW")
ax1c.text(xmax - dt.timedelta(days=75), 600, "c) DSW + DLW")

# fmt: on

ax1a.set_ylim([-10, 400])
ax1b.set_ylim([200, 350])
ax1c.set_ylim([200, 650])

ax1a.set_xticklabels([])
ax1b.set_xticklabels([])
ax1c.set_xlabel("Month")

# Format the dates
date_form = DateFormatter("%m")
ax1c.xaxis.set_major_formatter(date_form)

xrange = [xmin - 2.5 * off, xmax + 2.5 * off]
for ax in axs:
    ax.set_xlim(xrange)

if SAVEFIGS:
    figname = figdir + "broadband_meas_bymonth_b"
    plt.savefig(figname + ".png", format="png", dpi=600)
    plt.savefig(figname + ".eps", format="eps")
