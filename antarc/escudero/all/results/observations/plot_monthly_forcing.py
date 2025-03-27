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
from antarc.escudero.all.results.get_forcing_monthly_meas import (
    get_forcing_monthly_meas,
)

# Parameter modules
from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import swd_params, lwd_params
from antarc.escudero.all.results.params import SAVEFIGS


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


monthly = get_forcing_monthly_meas(
    years,
    lwd_base_dir,
    swd_base_dir,
    lwd_fmt,
    swd_fmt,
    era_dir,
    era5fmt,
    lat,
    lon,
    pts_min,
)


# Setup for figures
colors = ["red", "orange", "green", "cyan", "blue", "purple", "magenta"]
xleg = dt.datetime(2020, 7, 20)
xleg3 = dt.datetime(2020, 4, 20)
xmin = dt.datetime(2020, 1, 1)
xmax = dt.datetime(2020, 12, 1)


mon = list(range(1, 13))
swd_bot = [0 for i in range(12)]
lwd_bot = [0 for i in range(12)]
i = 0
for year, color in zip(years, colors):
    label = str(year)
    swd_yr = monthly["swd"]["npts"][i]
    lwd_yr = monthly["lwd"]["npts"][i]

    swd_bot += swd_yr
    lwd_bot += lwd_yr

    i += 1


# fmt: off
mcolor = "darkturquoise"  # measurement color

dtm = np.array([dt.datetime(2020, x, 1) for x in range(1, 13)])
off = dt.timedelta(days=3)
wid = 0.87
lft = 0.11
height1 = .47
height2 = .37


# fmt: off
wid = 0.87
lft = 0.11
height1 = .47
height2 = .37

# 2-panel 1 column plot of monthly mean cloud forcing without ERA5
fig, axs = plt.subplots(2, 1, num=12, clear=True)  # , figsize=[5, 8])
axs[0].set_position([.12, .55, .86, .42])
axs[1].set_position([.12, .1, .86, .42])
ax1a = axs[0]
ax1b = axs[1]
ax1a.set_position([lft, 0.5, wid, height1])
ax1b.set_position([lft, 0.1, wid, height2])

# Plot the monthly averages over all years
for i, key in enumerate(["swd", "lwd"]):
    ax = axs[i]
    ax.plot(dtm, np.nanmean(monthly[key]["meas"], axis=0), color=mcolor, linewidth=4, label='Mean')

i = 0
for year, color in zip(years, colors):
    label = str(year)
    ax1a.plot(dtm, monthly["swd"]["meas"][i], "+", color=color, label=label)
    ax1b.plot(dtm , monthly["lwd"]["meas"][i], "+", color=color)
    i += 1

ax1a.set_ylim([-210, 10])
ax1b.set_ylim([50, 78])

ax1a.legend(bbox_to_anchor=(0.38, .85), loc="upper left")

ax1a.set_ylabel("DSW Forcing (W/m$^{2}$)")
ax1a.yaxis.set_label_coords(-0.085, 0.54)
ax1b.set_ylabel("DLW Forcing (W/m$^{2}$)")
ax1b.yaxis.set_label_coords(-0.08, 0.5)

ax1b.xaxis.set_major_formatter(DateFormatter("%m"))

ax1a.text(xmax - dt.timedelta(days=15), -20, "a)")
ax1b.text(xmax - dt.timedelta(days=15), 73, "b)")

# Add labels for lines and symbols for Meas and ERA5 fluxes
ax = ax1a
xleg = dt.datetime(2020, 7, 5)
ax.plot([xleg, xleg + dt.timedelta(days=12)], [220, 220], 'k', linewidth=4)
ax.plot([xleg, xleg + dt.timedelta(days=12)], [190, 190], color=mcolor, linewidth=4)
# fmt: on

ax1a.set_xticklabels([])
ax1b.set_xlabel("Month")


xrange = [xmin - 2.5 * off, xmax + 2.5 * off]
for ax in axs:
    ax.set_xlim(xrange)

if SAVEFIGS:
    figname = figdir + "forcing_meas"
    plt.savefig(figname + ".png", format="png", dpi=600)
    plt.savefig(figname + ".eps", format="eps")
# fmt: on


# As above but also including total
# fmt: off
wid = 0.84
lft = 0.13
height1 = .27
height2 = .2
height3 = .3

# 2-panel 1 column plot of monthly mean cloud forcing without ERA5
fig, axs = plt.subplots(3, 1, num=12, clear=True)  # , figsize=[5, 8])
ax1a = axs[0]
ax1b = axs[1]
ax1c = axs[2]
ax1a.set_position([lft, 0.7, wid, height1])
ax1b.set_position([lft, 0.45, wid, height2])
ax1c.set_position([lft, 0.1, wid, height3])

# Plot the monthly averages over all years
for i, key in enumerate(["swd", "lwd"]):
    ax = axs[i]
    ax.plot(dtm, np.nanmean(monthly[key]["meas"], axis=0), color=mcolor, linewidth=4, label='Mean')
axs[2].plot(dtm, np.nanmean(monthly['swd']["meas"], axis=0)+ np.nanmean(monthly['lwd']["meas"], axis=0), color=mcolor, linewidth=4, label='Mean')
axs[2].plot(dtm, 0*np.arange(len(dtm)),'k--')

i = 0
for year, color in zip(years, colors):
    label = str(year)
    ax1a.plot(dtm, monthly["swd"]["meas"][i], "+", color=color, label=label)
    ax1b.plot(dtm, monthly["lwd"]["meas"][i], "+", color=color)
    ax1c.plot(dtm, monthly["swd"]["meas"][i] + monthly["lwd"]["meas"][i], "+", color=color)
    i += 1

ax1a.set_ylim([-210, 10])
ax1b.set_ylim([40, 90])
ax1c.set_ylim([-150, 100])

ax1a.legend(bbox_to_anchor=(0.3, .77), loc="upper left", ncol=2)

ax1a.set_ylabel("Forcing (W/m$^{2}$)")
ax1a.yaxis.set_label_coords(-0.085, 0.54)
ax1b.set_ylabel("Forcing (W/m$^{2}$)")
ax1b.yaxis.set_label_coords(-0.08, 0.5)
ax1c.set_ylabel("Forcing (W/m$^{2}$)")

ax1c.xaxis.set_major_formatter(DateFormatter("%m"))
# ax1c.yaxis.set_minor_locator([-150, -100, -50])

ax1a.text(xmax - dt.timedelta(days=40), -20, "a) DSW")
ax1b.text(xmax - dt.timedelta(days=40), 80, "b) DLW")
ax1c.text(xmax - dt.timedelta(days=65), 70, "c) DSW + DLW")

# Add labels for lines and symbols for Meas and ERA5 fluxes
ax = ax1a
xleg = dt.datetime(2020, 7, 5)
ax.plot([xleg, xleg + dt.timedelta(days=12)], [220, 220], 'k', linewidth=4)
ax.plot([xleg, xleg + dt.timedelta(days=12)], [190, 190], color=mcolor, linewidth=4)
# fmt: on

ax1a.set_xticklabels([])
ax1b.set_xticklabels([])
ax1c.set_xlabel("Month")

xrange = [xmin - 0.8 * off, xmax + 0.8 * off]
for ax in axs:
    ax.set_xlim(xrange)

if SAVEFIGS:
    figname = figdir + "forcing_meas"
    plt.savefig(figname + ".png", format="png", dpi=600)
    plt.savefig(figname + ".eps", format="eps")
# fmt: on
