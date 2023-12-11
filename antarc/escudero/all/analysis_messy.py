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
import pytz
import matplotlib.pyplot as plt
import calendar
from matplotlib.dates import DateFormatter

# from os.path import exists
# import pandas as pd


# My modules
# from antarc.run_radtran import get_clear_fluxes
from antarc.load_rad_flux import load_rad_flux, load_rad_flux_utc
from antarc.escudero.get_lwd_flux import get_lwd_clear_date

# from antarc.escudero.all.get_atm_profs import get_atm_profs

# from antarc.escudero.get_swd_flux import get_obs_swd, do_run_libradtran
# from antarc.escudero.get_era5_from_web import get_era5_from_web
# from antarc.escudero.all.get_met_from_web import get_met_15_minutes
# from antarc.escudero.all.get_era5_rad_from_web import get_era5_rad_from_web

# from antarc.escudero.case202205.make_plots import plot_measured_and_clear
# from antarc.escudero.case202205.make_plots import plot_forcings
# from antarc.escudero.case202205.make_plots import plot_meas_clear_forcings
# from antarc.escudero.case202205.make_plots import plot_forcings_and_lidar
# from antarc.escudero.all.read_era5_bbnd import read_era5_broadband_down


# Parameter modules specific to the user: the user must create params.py
# in the directory antarc with the variables MEAS_DIR and PROJ_DIR set
from antarc import params

# Parameter modules
from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import swd_params, lwd_params
from antarc.escudero.all.analysis_rsrc import (
    get_era5_esc_utc,
    ave_over_hour,
    get_stats,
)

#    qc_era5_date,


def get_lwd(date1, date2):
    # Load in measured LWD
    # if date2.year == date1.year + 1:
    #     date2_year1 = dt.datetime(date1.year, 12, 31, 23, 59, 59, tzinfo=pytz.utc)
    #     date1_year2 = dt.datetime(date2.year, 1, 1, tzinfo=pytz.utc)
    #     lwd_date, lwd = load_rad_flux(lwd_dir, lwd_fmt, date1, date2_year1)
    #     lwd_dateb, lwdb = load_rad_flux(lwd_dir, lwd_fmt, date1, date2)
    # else:

    lwd_dir = f"{lwd_base_dir}{date1.strftime('%Y')}/"

    lwd_date, lwd = load_rad_flux_utc(lwd_dir, lwd_fmt, date1, date2)
    # Quality control
    lwd[lwd < 130] = np.nan
    lwd[lwd > 500] = np.nan

    return lwd_date, lwd


def get_lwd_clear(date1, date2):
    # Load in clear-sky LWD calculated with LBLRTM
    lwd_clear_date, lwd_clear = get_lwd_clear_date(
        lwd_params.LWD_CLEAR_DIR,
        lwd_params.LWD_CLEAR_FILEFORMAT_OUTNU,
        lwd_params.LWD_CLEAR_FILEFORMAT_LOWNU,
        lwd_params.LWD_CLEAR_FILEFORMAT_HINU,
        date1,
        date2,
    )
    # lwd_tstamp = [x.timestamp() for x in lwd_date]
    # lwd_clear_tstamp = [x.timestamp() for x in lwd_clear_date]
    lwd_clear = np.array(lwd_clear)
    # lwd_clear_interp = np.interp(lwd_tstamp, lwd_clear_tstamp, lwd_clear)

    # if len(lwd) > 0:
    #     lwd_interp = np.interp(lwd_clear_tstamp, lwd_tstamp, lwd)
    # else:
    #     lwd_interp = []

    return lwd_clear_date, lwd_clear


def get_swd(date1, date2):

    # Load pyranometer data
    swd_dir = f"{swd_params.SWD_DIR}{date1.strftime('%Y')}/"
    swd_fmt = swd_params.SWD_FILEFORMAT
    swd_date0, swd = load_rad_flux_utc(swd_dir, swd_fmt, date1, date2)
    swd_date = np.array(swd_date0)
    # Quality control
    swd[swd < 0] = 0

    # # Get libradtran results for SWD clear
    # swd_clear = []
    # libdate = []
    # this_date = date1
    # while this_date <= date2:
    #     # Set the file to load in and make sure it exists
    #     clrfile = f"{swd_clear_dir}libradtran_{this_date.strftime('%Y%m%d')}.csv"
    #     if not exists(clrfile):
    #         print(f"No clear sky SWD file for {clrfile}")
    #     else:
    #         # Load the libradtran results and get the values
    #         libclear = pd.read_csv(clrfile, delimiter=",")
    #         swd_clear_date0 = pd.to_datetime(libclear["Date"]).to_list()
    #         swd_clear0 = libclear["global"].to_list()

    #         # Append the values to the result for this date
    #         libdate += swd_clear_date0
    #         swd_clear += swd_clear0

    #     # Set date for next time around
    #     this_date += dt.timedelta(days=1)

    # lib_tstamp = [x.timestamp() for x in libdate]
    # swd_clear_tstamp = [x.timestamp() for x in swd_clear_date]
    # swd_tstamp = [x.timestamp() for x in swd_date]
    # # swd_clear = np.interp(swd_tstamp, lib_tstamp, swd_clear)
    # swd_interp = np.interp(swd_clear_tstamp, swd_tstamp, swd)

    # Prove that ERA5 accumulates to the end of the hour
    # Thus we need to subtract 0.5 hours when directly comparing
    # (but not for hour aves, because I did the same in creating them)
    # qc_era5_date(era5, swd_date, swd)

    return swd_date, swd


def get_hourly_meas(lwd_date, lwd, swd_date, swd):

    nlwd_ave, lwd_date_ave, lwd_ave = ave_over_hour(lwd_date, lwd)
    nswd_ave, swd_date_ave, swd_ave = ave_over_hour(swd_date, swd)

    hourly = {
        "lwd_date": lwd_date_ave,
        "lwd": lwd_ave,
        "lwd_pts": nlwd_ave,
        "swd_date": swd_date_ave,
        "swd": swd_ave,
        "swd_pts": nswd_ave,
    }

    return hourly


def get_broadband_stats(era5, meas):

    if not era5 or (len(meas["lwd_date"]) == 0 and len(meas["swd_date"]) == 0):
        return None

    mean_lwd, rms_lwd, pts_lwd, rms_pc_lwd, absdiff_pc_lwd = get_stats(
        era5["date"], era5["lwd"], meas["lwd_date"], meas["lwd"]
    )
    mean_swd, rms_swd, pts_swd, rms_pc_swd, absdiff_pc_swd = get_stats(
        era5["date"], era5["swd"], meas["swd_date"], meas["swd"]
    )

    print(f"For {date1} to {date2}")
    print(f"LWD difference: mean and rms: {mean_lwd}, {rms_lwd}")
    print(f"SWD difference: mean and rms: {mean_swd}, {rms_swd}")

    stats = {
        "swd": {
            "bias": mean_swd,
            "rms": rms_swd,
            "npts": pts_swd,
            "rms_pc": rms_pc_swd,
            "absdiff_pc": absdiff_pc_swd,
        },
        "lwd": {
            "bias": mean_lwd,
            "rms": rms_lwd,
            "npts": pts_lwd,
            "rms_pc": rms_pc_lwd,
            "absdiff_pc": absdiff_pc_lwd,
        },
    }

    return stats


# Run flags
SAVEFIGS = False
PLOTFIGS = False

# # Date ranges for field campaigns
# dt_ranges = [
#     [dt.datetime(2017, 1, 12), dt.datetime(2017, 2, 1)],
#     [dt.datetime(2017, 11, 27), dt.datetime(2018, 2, 21)],
#     [dt.datetime(2018, 7, 1), dt.datetime(2019, 7, 1)],
#     [dt.datetime(2019, 7, 1), dt.datetime(2020, 7, 1)],
#     [dt.datetime(2021, 7, 1), dt.datetime(2023, 1, 1)],
#     [dt.datetime(2023, 1, 1), dt.datetime(2023, 10, 25)],
# ]

# # Date range to process
# # date1 = dt.datetime(2022, 1, 1, 0, 0, tzinfo=pytz.utc)
# # date2 = dt.datetime(2022, 12, 31, 23, 59, tzinfo=pytz.utc)
# i = 3
# date1 = dt_ranges[i][0].replace(tzinfo=pytz.utc)
# date2 = dt_ranges[i][1].replace(tzinfo=pytz.utc)


# Directories and variables needed from main param file
meas_dir = params.MEAS_DIR
# out_dir = params.PROJECT_DIR + "case_studies/May_2022/figures/"

# Directories and variables needed from other param files
# year_str = date1.strftime("%Y")
lwd_base_dir = lwd_params.LWD_DIR
lwd_fmt = lwd_params.LWD_FILEFORMAT
lwd_clear_dir = meas_dir + "Escudero/lwd_clear/"
era_broadband_dir = f"{params.MEAS_DIR}Escudero/era5/broadband/"
era5fmt = "era5_esc_broadband_*[0-9]*.nc"

lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE


# Date range
date1 = dt.datetime(2018, 1, 1)  # , tzinfo=pytz.utc)
date2 = dt.datetime(2018, 2, 1)  # , tzinfo=pytz.utc)

# Load in the measured broadband data
lwd_date, lwd = get_lwd(date1, date2)
swd_date, swd = get_swd(date1, date2)

# Get hourly averages to compare to ERA5
hourly_meas = get_hourly_meas(lwd_date, lwd, swd_date, swd)

# Get ERA5 data
era5 = get_era5_esc_utc(era_broadband_dir, era5fmt, date1, date2, lat, lon)

# Compare ERA5 to hourly
plt.figure()
plt.subplot(211)
plt.plot(era5["date"], era5["lwd"], label="ERA5")
plt.plot(lwd_date, lwd, label="meas")
plt.plot(hourly_meas["lwd_date"], hourly_meas["lwd"], label="meas, ave")
plt.subplot(212)
plt.plot(era5["date"], era5["swd"], label="ERA5")
plt.plot(swd_date, swd, label="meas")
plt.plot(hourly_meas["swd_date"], hourly_meas["swd"], label="meas, ave")
plt.legend()


# Assuming no missing ERA5 data, get indices to matches with measurements
iera_lwd = []
i = 0
for date in hourly_meas["lwd_date"]:
    while era5["date"][i] < date:
        i += 1
    if era5["date"][i] == date:
        iera_lwd += [i]
        i += 1
    elif era5["date"][i] > date:
        raise ValueError("This should not happen!")

# Check out how it worked
diff = [
    (hourly_meas["lwd_date"][i] - era5["date"][iera_lwd][i]).total_seconds()
    for i in range(len(hourly_meas["lwd_date"]))
]

plt.figure()
plt.plot(hourly_meas["lwd_date"], diff)


iera_swd = []
i = 0
for date in hourly_meas["swd_date"]:
    while era5["date"][i] < date:
        i += 1
    if era5["date"][i] == date:
        iera_swd += [i]
        i += 1
    elif era5["date"][i] > date:
        raise ValueError("This should not happen!")


# Now we can plot the differences
plt.figure()
plt.plot(hourly_meas["lwd_date"], era5["lwd"][iera_lwd] - hourly_meas["lwd"])
plt.plot(hourly_meas["swd_date"], era5["swd"][iera_swd] - hourly_meas["swd"])

plt.figure()
plt.plot(
    era5["swd"][iera_swd], era5["swd"][iera_swd] - hourly_meas["swd"], "."
)

plt.figure()
plt.subplot(211)
plt.plot(era5["swd"][iera_swd], hourly_meas["swd"], ".")
plt.plot([0, 850], [0, 850], "k--")
plt.xlabel("ERA5 SWD")
plt.ylabel("Meas SWD")
plt.subplot(212)
plt.plot(era5["lwd"][iera_lwd], hourly_meas["lwd"], ".")
plt.plot([220, 350], [220, 350], "k--")
plt.xlabel("ERA5 LWD")
plt.ylabel("Meas LWD")

# # Annual differences
# years = range(2017, 2024)
# keys = ["mean_lwd", "mean_swd", "rms_lwd", "rms_swd"]
# yearly = {}
# for key in keys:
#     yearly[key] = []

# for year in years:
#     date1 = dt.datetime(year, 1, 1).replace(tzinfo=pytz.utc)
#     date2 = dt.datetime(year, 12, 31, 23, 59).replace(tzinfo=pytz.utc)

#     year_str = date1.strftime("%Y")
#     lwd_dir = lwd_base_dir + year_str + "/"
#     # swd_clear_dir = f"{swd_params.SWD_CLEAR_DIR}{year_str}/"

#     stats = get_stats_daterange(date1, date2)
#     for key in keys:
#         yearly[key].append(stats[key])

# plt.figure()
# plt.plot(years, yearly["rms_swd"], "r", label="RMS Diff: SWD")
# plt.plot(years, yearly["mean_swd"], "r--", label="Bias: ERA5 - Meas; SWD")
# plt.plot(years, yearly["rms_lwd"], "b", label="RMS Diff: LWD")
# plt.plot(years, yearly["mean_lwd"], "b--", label="Bias: ERA5 - Meas; LWD")
# plt.legend()


# Monthly stats

# Monthly differences
years = range(2017, 2024)
months = range(1, 13)


# stats_keys = {
#     "lwd": ["bias", "rms", "npts"],
#     "swd": ["bias", "rms", "npts"],
# }
monthly = {
    "lwd": {
        "bias": np.nan * np.ones(len(years) * len(months)),
        "rms": np.nan * np.ones(len(years) * len(months)),
        "npts": np.zeros(len(years) * len(months)),
        "rms_pc": np.nan * np.ones(len(years) * len(months)),
        "absdiff_pc": np.nan * np.ones(len(years) * len(months)),
    },
    "swd": {
        "bias": np.nan * np.ones(len(years) * len(months)),
        "rms": np.nan * np.ones(len(years) * len(months)),
        "npts": np.zeros(len(years) * len(months)),
        "rms_pc": np.nan * np.ones(len(years) * len(months)),
        "absdiff_pc": np.nan * np.ones(len(years) * len(months)),
    },
}


dates = []
count = 0
for year in years:
    for month in months:
        dtime = dt.datetime(year, month, 1)
        dates.append(dtime.replace(tzinfo=pytz.utc))
        count += 1


pts_min = 7  # minimum number of points to compute stats
i1 = 10  # :84
for i, date1 in enumerate(dates[i1:20]):
    idate = i1 + i
    print(idate)
    year = date1.year
    mon = date1.month
    day = calendar.monthrange(year, mon)[1]
    date2 = dt.datetime(year, mon, day, 23, 59).replace(tzinfo=pytz.utc)

    year_str = date1.strftime("%Y")
    lwd_dir = lwd_base_dir + year_str + "/"
    # swd_clear_dir = f"{swd_params.SWD_CLEAR_DIR}{year_str}/"

    # Load in the measured broadband data
    lwd_date, lwd = get_lwd(date1, date2)
    swd_date, swd = get_swd(date1, date2)

    # Get hourly averages to compare to ERA5
    hourly_meas = get_hourly_meas(lwd_date, lwd, swd_date, swd)

    # Get ERA5 data and convert date to timezone-aware
    era5 = get_era5_esc_utc(era_broadband_dir, era5fmt, date1, date2, lat, lon)

    stats = get_broadband_stats(era5, hourly_meas)

    if not stats:
        continue

    for wave in ("lwd", "swd"):
        if stats[wave]["npts"] >= pts_min:
            for key in stats[wave].keys():
                monthly[wave][key][idate] = stats[wave][key]
        else:
            monthly[wave]["npts"][idate] = stats[wave]["npts"]


# make into numpy arrays
# mon_arr = {}
# for key in keys:
#     mon_arr[key] = np.array(monthly[key])

# Setup for figures
colors = ["red", "orange", "green", "cyan", "blue", "purple", "black"]
months = [x.month for x in dates]
dtimes = np.array(dates)
xleg = dt.datetime(2020, 7, 20)
xmin = dt.datetime(2020, 1, 1)
xmax = dt.datetime(2020, 12, 1)

fig, (ax1a, ax1b) = plt.subplots(2, 1)
fig, (ax2a, ax2b) = plt.subplots(2, 1)
ax1a.set_position([0.11, 0.56, 0.73, 0.41])
ax1b.set_position([0.11, 0.08, 0.73, 0.41])
ax2a.set_position([0.1, 0.56, 0.74, 0.41])
ax2b.set_position([0.1, 0.08, 0.74, 0.41])

for year, color in zip(years, colors):
    i = np.array([j for j, x in enumerate(dates) if x.year == year])
    if len(i) == 0:
        continue
    dtime = [dt.datetime(2020, x.month, 1) for x in dtimes[i]]
    ax1a.plot(dtime, monthly["swd"]["rms"][i], color=color, label=str(year))
    ax1a.plot(dtime, monthly["swd"]["bias"][i], "--", color=color)
    ax1b.plot(dtime, monthly["lwd"]["rms"][i], color=color)
    ax1b.plot(dtime, monthly["lwd"]["bias"][i], "--", color=color)

    ax2a.plot(dtime, monthly["swd"]["npts"][i], color=color, label=str(year))
    ax2b.plot(dtime, monthly["lwd"]["npts"][i], color=color)

ax1a.plot([xleg, xleg + dt.timedelta(days=12)], [148, 148], "k")
ax1a.text(xleg + dt.timedelta(days=15), 144, "RMS diff")
ax1a.plot([xleg, xleg + dt.timedelta(days=12)], [131, 131], "k--")
ax1a.text(xleg + dt.timedelta(days=15), 126, "Bias: ERA - Meas")
ax1a.set_ylabel("SWD Flux (W/m$^{2}$)")
ax1a.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left")
ax1a.xaxis.set_major_formatter(DateFormatter("%m"))
ax1b.set_ylabel("LWD Flux (W/m$^{2}$)")
ax1b.xaxis.set_major_formatter(DateFormatter("%m"))

ax2a.plot([xmin, xmax], [pts_min, pts_min], "k:")
ax2a.set_ylabel("Number of SWD points")
ax2a.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left")
ax2a.xaxis.set_major_formatter(DateFormatter("%m"))
ax2a.set_ylim([0, 800])
ax2a.set_xlim([xmin, xmax])

ax2b.plot([xmin, xmax], [pts_min, pts_min], "k:")
ax2b.set_ylabel("Number of LWD points")
ax2b.xaxis.set_major_formatter(DateFormatter("%m"))
ax2b.set_ylim([0, 800])
ax2b.set_xlim([xmin, xmax])


# TODO: save figures in png and eps formats
if SAVEFIGS:
    figname = "broadband_meas_v_era5"
    plt.savefig(figname + ".png", format="png")
    plt.savefig(figname + ".eps", format="eps")


# fig, (ax5, ax6) = plt.subplots(2, 1)
# ax5.set_position([0.1, 0.56, 0.74, 0.41])
# ax6.set_position([0.1, 0.08, 0.74, 0.41])

# for year, color in zip(years, colors):
#     i = np.array([j for j, x in enumerate(dates) if x.year == year])
#     if len(i) == 0:
#         continue
#     dtime = [dt.datetime(2020, x.month, 1) for x in dtimes[i]]
#     ax5.plot(dtime, monthly["swd"]["rms_pc"][i], color=color, label=str(year))
#     ax5.plot(dtime, monthly["swd"]["absdiff_pc"][i], "--", color=color)
#     ax6.plot(dtime, monthly["lwd"]["rms_pc"][i], color=color)
#     ax6.plot(dtime, monthly["lwd"]["absdiff_pc"][i], "--", color=color)

# ax5.plot([xleg, xleg + dt.timedelta(days=12)], [148, 148], "k")
# ax5.text(xleg + dt.timedelta(days=15), 144, "RMS diff")
# ax5.plot([xleg, xleg + dt.timedelta(days=12)], [131, 131], "k--")
# ax5.text(xleg + dt.timedelta(days=15), 126, "Bias: ERA - Meas")
# ax5.set_ylabel("SWD Flux Diff (%)")
# ax5.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left")
# ax5.xaxis.set_major_formatter(DateFormatter("%m"))
# ax6.set_ylabel("LWD Flux Diff (%)")
# ax6.xaxis.set_major_formatter(DateFormatter("%m"))


# plt.figure()
# plt.plot(months, monthly["mean_swd"], "r--", label="Bias: ERA5 - Meas; SWD")
# plt.plot(months, monthly["rms_lwd"], "b", label="RMS Diff: LWD")
# plt.plot(months, monthly["mean_lwd"], "b--", label="Bias: ERA5 - Meas; LWD")
# plt.legend()

# plt.figure()
# plt.plot(dates, monthly["rms_swd"], "r", label="RMS Diff: SWD")
# plt.plot(dates, monthly["mean_swd"], "r--", label="Bias: ERA5 - Meas; SWD")
# plt.plot(dates, monthly["rms_lwd"], "b", label="RMS Diff: LWD")
# plt.plot(dates, monthly["mean_lwd"], "b--", label="Bias: ERA5 - Meas; LWD")
# plt.legend()


# if PLOTFIGS:
#     # Plot results for fine-resolution measurements
#     halfhour = dt.timedelta(minutes=30)

#     plt.figure(2)
#     plt.clf()

#     plt.subplot(211)
#     plt.plot(swd_date, swd, label="measured")
#     plt.plot(era5["date"] - halfhour, era5["swd"], label="era5")
#     plt.plot(
#         era5["date"] - halfhour, era5["swd_clr"], "k--", label="era5 clear"
#     )
#     plt.ylabel("Downward shortwave (W/m$^{2}$)")
#     plt.legend()

#     plt.subplot(212)
#     plt.plot(lwd_date, lwd)
#     plt.plot(era5["date"] - halfhour, era5["lwd"])
#     plt.plot(era5["date"] - halfhour, era5["lwd_clr"], "k--")
#     plt.plot(lwd_clear_date, lwd_clear, "o", label="this work; clear")
#     plt.ylabel("Downward longwave (W/m$^{2}$)")
#     plt.legend()

#     plt.figure(3)
#     plt.clf()

#     plt.subplot(211)
#     # plt.plot(swd_date, swd, label="measured")
#     plt.plot(swd_date_ave, swd_ave, label="measured, ave")
#     plt.plot(era5["date"], era5["swd"], label="era5")
#     plt.plot(era5["date"], era5["swd_clr"], "k--", label="era5 clear")
#     plt.ylabel("Downward shortwave (W/m$^{2}$)")
#     plt.legend()

#     plt.subplot(212)

#     # plt.plot(lwd_date, lwd, label="measured")
#     plt.plot(lwd_date_ave, lwd_ave, label="measured, ave")
#     plt.plot(era5["date"], era5["lwd"], label="era5")
#     plt.plot(era5["date"], era5["lwd_clr"], "k--", label="era5, clear")
#     plt.plot(lwd_clear_date, lwd_clear, "o", label="this work; clear")
#     plt.ylabel("Downward longwave (W/m$^{2}$)")
#     plt.legend()
