#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 08:56:46 2022

@author: prowe

Purpose: Plot cloud forcing at Escudero for 2017-2023

"""

# Dependencies
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt


# My modules
from antarc.escudero.get_lwd_flux import get_lwd_clear_date


# Parameter modules specific to the user: the user must create params.py
# in the directory antarc with the variables MEAS_DIR and PROJ_DIR set
from antarc import params

# Parameter modules
from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import lwd_params
from antarc.escudero.all.read_era5_bbnd import read_era5_broadband_down
from antarc.escudero.get_frei_15min import get_frei_15

# TODO: get ERA5 data from web again for 2023/08/01 on b/c it was non-final
# as of 2023/10/25. Get the ERA5 data from the web  (see get_cloud_forcing)


# Directories
MEAS_DIR = params.MEAS_DIR
dir_lwdclear = MEAS_DIR + "Escudero/lwd_clear/"
out_dir = params.PROJECT_DIR + "case_studies/May_2022/figures/"
era_broadband_dir = f"{params.MEAS_DIR}Escudero/era5/broadband/"
surf_dir = f"{params.MEAS_DIR}Escudero/aws/std/"
frei_direc = f"{params.MEAS_DIR}Escudero/frei_met/by_month_every_15min/"


# File formats
erafmt = "era5_esc_broadband_*[0-9]*.nc"
surf_samplefile = "tarp02_met_202305.txt"
frei_sampfname = "esc_met_201803.txt"
frei_pfx = 7


# Select the date range
# start_date = dt.datetime(2017, 1, 1, 0, 0)
# end_date = dt.datetime(2023, 12, 31, 23, 59)
start_date = dt.datetime(2021, 9, 1, 0, 0)
end_date = dt.datetime(2023, 8, 31, 23, 59)


# Params needed from param files
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
era_lat = -62.25
era_lon = -59.0


# Load in the surface temperatures from FREI
(
    fdtime,
    ftemp,
    fpress,
    fslp,
    fwdir,
    fwspd,
    frh,
) = get_frei_15(frei_direc, frei_sampfname, frei_pfx)


# Load in clear scene files
date1 = dt.datetime(2017, 1, 1, 0, 0)
date2 = dt.datetime(2023, 8, 31, 23, 59)
lwd_clear_date, lwd_clear = get_lwd_clear_date(
    lwd_params.LWD_CLEAR_DIR,
    lwd_params.LWD_CLEAR_FILEFORMAT_OUTNU,
    lwd_params.LWD_CLEAR_FILEFORMAT_LOWNU,
    lwd_params.LWD_CLEAR_FILEFORMAT_HINU,
    date1,
    date2,
)

era_keys = [
    "swd",
    "lwd",
    "swd_clr",
    "lwd_clr",
    "sw_net",
    "lw_net",
    "sw_net_clr",
    "lw_net_clr",
]
era = {"date": np.array([])}
for key in era_keys:
    era[key] = np.array([])

years = np.arange(start_date.year, end_date.year + 1)

for year in years:
    # Get ERA5 data for the entire year
    date1 = dt.datetime(year, 1, 1, 0, 0)
    date2 = dt.datetime(year, 12, 31, 23, 59)

    # If we are in the first or last year, use the specified start/end
    if year == start_date.year:
        date1 = start_date
    if year == end_date.year:
        date2 = end_date

    # Useful variables
    year_str = date1.strftime("%Y")

    # QC on input dates
    if year_str != date2.strftime("%Y"):
        raise ValueError("Different years must be processed separately")

    # Get ERA5 data - recommended to use for cloud forcing
    era5rad_outdir = f"{era_broadband_dir}{date1.strftime('%Y')}/"

    # Read in the ERA5 data
    era5bb = read_era5_broadband_down(era5rad_outdir, erafmt, date1, date2)
    if era_lat not in era5bb["lat"].data or era_lon not in era5bb["lon"].data:
        raise ValueError("Bad latititude or longitude")
    # Get closest grid point
    i = np.where(era5bb["lat"] == era_lat)[0][0]
    j = np.where(era5bb["lon"] == era_lon)[0][0]

    era["date"] = np.hstack([era["date"], era5bb["date"].data])
    for key in era_keys:
        era[key] = np.hstack([era[key], era5bb[key][:, i, j].data])


# Loop over era dates and find all hours for which max surface temp > 0 C
if fdtime[0] > era["date"][0]:
    raise ValueError("First surface temp must be before first ERA5 date")

# Add 0.9 K to the Frei temperatures to match Escudero
ftemp += 0.9

# Find the maximum temperature within each hour at Frei
# Because ERA5 accumulates to the hour, for the first hour
# just use the first measurement, while for all the rest
# average over first through last inclusive


i = 0
maxtemp = np.nan * np.ones(len(era["date"]))
for iera, dtime0 in enumerate(era["date"]):
    # Skip to the beginning of the hour
    while fdtime[i] < dtime0 - dt.timedelta(hours=1):
        i += 1
    # Get max time over next hour
    while fdtime[i] <= dtime0:
        i += 1
        maxtemp[iera] = np.nanmax([maxtemp[iera], ftemp[i]])
    # Dial back 1 because we double count the last/first
    i -= 1


# Plot the maximum temperature for each hour
plt.figure(num=1)
plt.plot(fdtime, ftemp, label="Surface temperature, Frei (C)")
plt.plot(era["date"], maxtemp, label="Maximum temperature over hour (C)")
plt.xlabel("Date")
plt.ylabel("Temperature (C)")
plt.legend()

# Subtract 1/2 hour from ERA5 times
era["date"] -= dt.timedelta(minutes=30)
era_date = era["date"]

# Only use times with max temp over some temperature threshold (in C)
iwarm = np.where(maxtemp > 2)[0]
era_date = era["date"][iwarm]
era_warm = {}
for key in era_keys:
    era_warm[key] = era[key][iwarm]
era = era_warm

mon = list(range(1, 13))

# Get ERA5 monthly mean values by year
monthly = {}
for key in era_keys:
    monthly[key] = np.zeros([12, len(years)])
for iyear, year in enumerate(years):
    iyr = [i for i, x in enumerate(era_date) if x.year == year]

    dtime0 = era_date[iyr]
    for key in era_keys:
        val = era[key][iyr]
        for imonth in range(12):
            mon = imonth + 1
            imon = [i for i, x in enumerate(dtime0) if x.month == mon]
            monthly[key][imonth, iyear] = np.nanmean(val[imon])


# fmt: off
# Plot ERA5 broadband components
swu = np.nanmean(monthly["sw_net"] - monthly["swd"], axis=1)
lwu = np.nanmean(monthly["lw_net"] - monthly["lwd"], axis=1)
swu_clr = np.nanmean(monthly["sw_net_clr"] - monthly["swd_clr"], axis=1)
lwu_clr = np.nanmean(monthly["lw_net_clr"] - monthly["lwd_clr"], axis=1)

mon = list(range(1, 13))

fig, (ax1, ax2) = plt.subplots(1, 2, num=2, clear=True, figsize=(9,6)) #,sharex=True)
ax1.set_position([0.58, 0.435, 0.4, 0.546])
ax2.set_position([0.09, 0.08, 0.4, 0.9])
ax = ax1
ax.plot([0,13],[0,0], color='gray')
ax.plot(mon, np.nanmean(monthly['swd_clr'], axis=1), 'c--', label="clear, down")
ax.plot(mon, np.nanmean(monthly['sw_net_clr'], axis=1), 'c', label="clear, net")
ax.plot(mon, np.nanmean(monthly['swd'], axis=1), 'k--', label="scene, down")
ax.plot(mon, np.nanmean(monthly['sw_net'], axis=1), 'r',label="scene, net")
ax.plot(mon, swu, 'k:', label="scene, up")
ax.plot(mon, swu_clr, 'c:', label="clear, up")

ax.legend() #loc='center right', bbox_to_anchor=(1,.25,.5,.5))
ax.set_ylabel("Shortwave Flux (W/m$^{2}$)")
ax.set_xlim([.8, 12.2])
ax.set_ylim([-50, 400])
ax.set_xlabel("Month")

ax = ax2
ax.plot([0,13],[0,0], color='gray')
ax.plot(mon, np.nanmean(monthly['lwd'], axis=1), 'k--', label='scene, down')
ax.plot(mon, np.nanmean(monthly['lwd_clr'], axis=1), 'c--', label='clear, down')
ax.set_prop_cycle(None)
ax.set_ylabel("Longwave Flux (W/m$^{2}$)")
ax.plot(mon, np.nanmean(monthly['lw_net'], axis=1), color='r', label='scene, net')
ax.plot(mon, np.nanmean(monthly['lw_net_clr'], axis=1), 'c', label='clear, net')
ax.plot(mon, lwu, 'k:', label="scene, up")
ax.plot(mon, lwu_clr, 'c:', linewidth=4,label="clear, up")

for imonth in range(12):
    imon = [i for i, x in enumerate(era_date) if x.month == imonth + 1]
    ax.boxplot(era['lwd'][imon], positions=[mon[imonth]],whis=(0, 100))
    ax.boxplot(era['lwd_clr'][imon], positions=[mon[imonth]+.2],whis=(0, 100))
    # ax.boxplot(era['lw_net'][imon]-era['lwd'][imon], positions=[mon[imonth]],whis=(0, 100))
    # ax.violinplot(era['lw_net'][imon]-era['lwd'][imon], positions=[mon[imonth]])
    # ax.boxplot(era['lw_net_clr'][imon]-era['lwd_clr'][imon], positions=[mon[imonth]+.2],whis=(0, 100))
    # ax.boxplot(era['lw_net'][imon], positions=[mon[imonth]],whis=(0, 100))
    ax.boxplot(era['lw_net_clr'][imon], positions=[mon[imonth]+.2],whis=(0, 100))

ax2.legend(loc='lower right', bbox_to_anchor=(1.,0,.5,.5)) 
ax.set_xlabel("Month")
ax.set_ylabel("Longwave Flux (W/m$^{2}$)")
ax.set_ylim([-350, 400])
ax.set_xticks([2,4,6,8,10,12])
ax.set_xticklabels([2,4,6,8,10,12])
ax.set_xlim([.8, 12.2])
# fmt: on


# fmt: off
# Without upwelling, no box and whiskers
fig, (ax1, ax2,ax3) = plt.subplots(1, 3, num=3, clear=True, figsize=(12,6))
ax1.set_position([0.07, 0.08, 0.25, 0.9])
ax2.set_position([0.4, 0.08, 0.25, 0.9])
ax3.set_position([0.72, 0.08, 0.25, 0.9])
ax = ax2
ax.plot([0,13],[0,0], color='gray')
ax.plot(mon, np.nanmean(monthly['swd_clr'], axis=1), 'c--', label="clear, down")
ax.plot(mon, np.nanmean(monthly['sw_net_clr'], axis=1), 'c', label="clear, net")
ax.plot(mon, np.nanmean(monthly['swd'], axis=1), 'k--', label="scene, down")
ax.plot(mon, np.nanmean(monthly['sw_net'], axis=1), 'r',label="scene, net")
ax.legend() #loc='center right', bbox_to_anchor=(1,.25,.5,.5))
ax.set_ylabel("Shortwave Flux (W/m$^{2}$)")
ax.set_xlim([.8, 12.2])
ax.set_ylim([-150, 600])
ax.set_xlabel("Month")

ax = ax1
ax.plot([0,13],[0,0], color='gray')
ax.plot(mon, np.nanmean(monthly['lwd'], axis=1), 'k--', label='scene, down')
ax.plot(mon, np.nanmean(monthly['lwd_clr'], axis=1), 'c--', label='clear, down')
ax.set_prop_cycle(None)
ax.set_ylabel("Longwave Flux (W/m$^{2}$)")
ax.plot(mon, np.nanmean(monthly['lw_net'], axis=1), color='r', label='scene, net')
ax.plot(mon, np.nanmean(monthly['lw_net_clr'], axis=1), 'c', label='clear, net')
ax.legend() #loc='center', bbox_to_anchor=(0.05,0.2,.5,.5)) 
ax.set_xlabel("Month")
ax.set_ylabel("Longwave Flux (W/m$^{2}$)")
ax.set_ylim([-150, 600])
ax.set_xticks([2,4,6,8,10,12])
ax.set_xticklabels([2,4,6,8,10,12])
ax.set_xlim([.8, 12.2])

ax = ax3
ax.plot([0,13],[0,0], color='gray')
ax.plot(mon, np.nanmean(monthly['swd_clr']+monthly['lwd_clr'], axis=1), 'c--', label="clear, down")
ax.plot(mon, np.nanmean(monthly['swd']+monthly['lwd'], axis=1), 'k--', label="scene, down")
ax.plot(mon, np.nanmean(monthly['sw_net_clr']+monthly['lw_net_clr'], axis=1), 'c', label="clear, net")
ax.plot(mon, np.nanmean(monthly['sw_net']+monthly['lw_net'], axis=1), 'r',label="scene, net")
ax.legend() #loc='center right', bbox_to_anchor=(1,.25,.5,.5))
ax.set_ylabel("Total Flux (W/m$^{2}$)")
ax.set_xlim([.8, 12.2])
ax.set_ylim([-150, 600])
ax.set_xlabel("Month")
# fmt: on





# fmt: off
# Without upwelling, with box and whiskers
fig, (ax1, ax2,ax3) = plt.subplots(1, 3, num=4, clear=True, figsize=(12,6))
ax1.set_position([0.07, 0.08, 0.25, 0.9])
ax2.set_position([0.4, 0.08, 0.25, 0.9])
ax3.set_position([0.72, 0.08, 0.25, 0.9])
ax = ax1
ax.plot([0,13],[0,0], color='gray')
ax.plot(mon, np.nanmean(monthly['lwd'], axis=1), 'k--', label='scene, down')
ax.plot(mon, np.nanmean(monthly['lwd_clr'], axis=1), 'c--', label='clear, down')
ax.set_prop_cycle(None)
ax.set_ylabel("Longwave Flux (W/m$^{2}$)")
ax.plot(mon, np.nanmean(monthly['lw_net'], axis=1), color='r', label='scene, net')
ax.plot(mon, np.nanmean(monthly['lw_net_clr'], axis=1), 'c', label='clear, net')
goup = False
for imonth in range(12):
    imon = [i for i, x in enumerate(era_date) if x.month == imonth + 1]
    ax.boxplot(era['lwd'][imon], positions=[mon[imonth]], whis=(0, 100))
    ax.boxplot(era['lwd_clr'][imon], positions=[mon[imonth]+.2], whis=(0, 100))
    ax.boxplot(era['lw_net'][imon], positions=[mon[imonth]], whis=(0, 100))
    if goup:
        ax.text(mon[imonth]-.4, 380, str(len(imon)))
        goup = False
    else:
        ax.text(mon[imonth]-.4, 350, str(len(imon)))
        goup = True
    # ax.boxplot(era['lw_net_clr'][imon], positions=[mon[imonth]+.2],whis=(0, 100))
    
ax.legend() #loc='center', bbox_to_anchor=(0.05,0.2,.5,.5)) 
ax.set_xlabel("Month")
ax.set_ylabel("Longwave Flux (W/m$^{2}$)")
ax1.set_ylim([-150, 860])
ax.set_xticks([2,4,6,8,10,12])
ax.set_xticklabels([2,4,6,8,10,12])
ax.set_xlim([.6, 12.4])

ax = ax2
ax.plot([0,13],[0,0], color='gray')
ax.plot(mon, np.nanmean(monthly['swd_clr'], axis=1), 'c--', label="clear, down")
ax.plot(mon, np.nanmean(monthly['sw_net_clr'], axis=1), 'c', label="clear, net")
ax.plot(mon, np.nanmean(monthly['swd'], axis=1), 'k--', label="scene, down")
ax.plot(mon, np.nanmean(monthly['sw_net'], axis=1), 'r',label="scene, net")
for imonth in range(12):
    imon = [i for i, x in enumerate(era_date) if x.month == imonth + 1]
    ax.boxplot(era['sw_net'][imon], positions=[mon[imonth]], whis=(0, 100))
ax.legend() #loc='center right', bbox_to_anchor=(1,.25,.5,.5))
ax.set_ylabel("Shortwave Flux (W/m$^{2}$)")
ax.set_xlim([.8, 12.2])
ax.set_ylim([-150, 860])
ax.set_xlabel("Month")

ax = ax3
ax.plot([0,13],[0,0], color='gray')
ax.plot(mon, np.nanmean(monthly['swd_clr']+monthly['lwd_clr'], axis=1), 'c--', label="clear, down")
ax.plot(mon, np.nanmean(monthly['swd']+monthly['lwd'], axis=1), 'k--', label="scene, down")
ax.plot(mon, np.nanmean(monthly['sw_net_clr']+monthly['lw_net_clr'], axis=1), 'c', label="clear, net")
ax.plot(mon, np.nanmean(monthly['sw_net']+monthly['lw_net'], axis=1), 'r',label="scene, net")
ax.legend() #loc='center right', bbox_to_anchor=(1,.25,.5,.5))

for imonth in range(12):
    imon = [i for i, x in enumerate(era_date) if x.month == imonth + 1]
    ax.boxplot(era['lw_net'][imon] + era['sw_net'][imon], positions=[mon[imonth]], whis=(0, 100))
    # ax.boxplot(era['lw_net_clr'][imon], positions=[mon[imonth]+.2],whis=(0, 100))

ax.set_ylabel("Total Flux (W/m$^{2}$)")
ax.set_xlim([.8, 12.2])
ax.set_ylim([-150, 860])
ax.set_xlabel("Month")
# fmt: on


# Histogram of Net Fluxes with Season
mon_seasons = [[12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]
binned = []
for imon in mon_seasons:
    iera_mon = []
    for imonth in imon:
        iera_mon += [i for i, x in enumerate(era_date) if x.month == imonth]
    binned.append(era["lw_net"][iera_mon] + era["sw_net"][iera_mon])

plt.figure(num=5)
plt.hist(binned)
plt.xlabel("Net Flux (W/m$^2$)")
plt.ylabel("Number of cases")
plt.legend(["DJF", "MAM", "JJA", "SON"])
