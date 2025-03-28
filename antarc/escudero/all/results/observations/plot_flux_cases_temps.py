#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:21:08 2022

@author: prowe

Notes:
For polar WRF:
    - If it is instantanous by hour, should average measurements from - 1/2 
      to + 1/2 hour, whereas here it is averaged over hour
    - It looks like the shortwave is slightly shifted right. Why?
    - The PWRF DSW is sometimes greater (e.g. 2022/06/02) than the clear. Why?
      This could be because the clear is averaged over hours, while PWRF is
      instaneous, and hits a higher point. If so, then should we *not* average
      the DSW measurements over hour to compare to the shortwave?
"""


import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter
import calendar
import pytz
from netCDF4 import Dataset
import os
import pandas as pd
from scipy.stats import linregress
from dateutil.relativedelta import relativedelta


from antarc import params
from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import swd_params, lwd_params
from antarc.escudero.parameters import esc202202
from antarc.escudero.all.results import get_era5_broadband_monthly_warm
from antarc.escudero.all.results import get_broadband_monthly_meas
from antarc.escudero.all.results.get_hourly_meas import (
    get_hourly_meas_fill as get_hourly_meas,
)
from antarc.escudero.all.analysis_rsrc import (
    get_lwd_qc,
    get_swd_qc,
)
from antarc.escudero.all.get_pwrf_bbnd import get_pwrf_bbnd_daterange
from antarc.escudero.all.results.plot_broadband_esc_pwrf_era5_case_rsrc import (
    get_era5_broadband_month_warm as get_era5_bbnd_mon_warm,
    get_era5_atgridpts,
)

from antarc.escudero.all.results.params import SAVEFIGS, SAVEDIR
from antarc.escudero.all.results.get_event_flux_rsrc import (
    get_event_flux,
    get_event_forcing,
)

get_broadband_monthly_meas = (
    get_broadband_monthly_meas.get_broadband_monthly_meas
)
get_era5_broadband_monthly_warm = (
    get_era5_broadband_monthly_warm.get_era5_broadband_monthly_warm
)


# from antarc.escudero.all.analysis_rsrc import ave_over_every_hour


def get_frei_temps(start, event_hours):

    # The various methods used to get the daily temperature
    temp_methods = ("Min-Max", "4 Sinopticas", "24 Horas", "8 Sinopticas")

    # Drop the start time back to the beginning of the day
    event_start = dt.datetime(start.year, start.month, start.day)

    # The event end should be at the end of the day
    event_end = event_start + dt.timedelta(hours=event_hours)

    # Loop over months
    interpolated_days = []
    dtm = []
    temp = []
    date = event_start
    while date < event_end:
        frei_file = frei_dir + event_start.strftime(temp_fmt)
        print(frei_file)

        # Load in the file
        dfm = pd.read_excel(frei_file, skiprows=1, skipfooter=1)

        day = dfm.iloc[:, 0]
        dfm["date"] = [dt.datetime(date.year, date.month, x) for x in day]
        mask = (dfm["date"] >= event_start) & (dfm["date"] <= event_end)

        event_df = dfm.loc[mask]

        # Check for any missing methods
        val = event_df["Valor"]
        for i, meth in enumerate(event_df["Metodo"]):
            if meth not in temp_methods:
                # Check the type of the value
                print(f"Found method for getting daily temperature: '{meth}'")
                print(f"with corresponding value of type {type(val[i])}")
                # Check if there is a value on the day before and the day
                # after. If so, just interpolate, but make a note
                if i == 0:
                    raise ValueError("i=0, modify code to use previous file")
                elif i >= len(event_df) - 1:
                    raise ValueError("Last i, modify code to use next file")
                elif (
                    event_df["Metodo"][i - 1] in temp_methods
                    and event_df["Metodo"][i + 1] in temp_methods
                ):
                    event_df.loc[i, "Valor"] = (val[i - 1] + val[i + 1]) / 2
                    interpolated_days += [event_df.loc[i, "date"]]
                else:
                    raise ValueError("Possibility not accounted for; see code")

        temp += event_df["Valor"].to_list()
        dtm += event_df["date"].to_list()

        # Increment start for next month
        date += relativedelta(months=1)

    return temp, dtm


def get_hourly_meas_over_years(
    years, mon, lwd_base_dir, lwd_fmt, swd_base_dir, swd_fmt, keys
):
    # For the case month, for all selected years, load and average measured
    # DLW/DSW hourly and append all results to the dictionary of lists meas_hr
    meas_hr0 = {key: [] for key in keys}
    for year in years:
        mdate1 = dt.datetime(year, mon, 1)
        day = calendar.monthrange(year, mon)[1]
        mdate2 = dt.datetime(year, mon, day, 23, 59)

        # Load in the measured broadband data
        lwd_date, lwd = get_lwd_qc(lwd_base_dir, lwd_fmt, mdate1, mdate2)
        swd_date, swd = get_swd_qc(swd_base_dir, swd_fmt, mdate1, mdate2)

        # Get hourly averages to compare to ERA5. Drop back an hour
        # to match ERA5, whose first hour comes from the accumulated hour
        # on the previous day
        dt1 = mdate1 - dt.timedelta(hours=1)
        dt2 = mdate2 - dt.timedelta(hours=1)
        meas0 = get_hourly_meas(lwd_date, lwd, swd_date, swd, dt1, dt2)

        # Only keep non-nan results
        ikeep = np.nonzero(~np.isnan(meas0["swd"]))[0]
        ikeep = np.intersect1d(ikeep, np.nonzero(~np.isnan(meas0["lwd"]))[0])
        for key in keys:
            meas_hr0[key] += [
                x for i, x in enumerate(meas0[key]) if i in ikeep
            ]

    # Drop the dates back 1/2 hour and get the total
    meas_hr0["date"] = [x - dt.timedelta(minutes=30) for x in meas_hr0["date"]]
    meas_hr0["tot"] = [x + y for x, y in zip(meas_hr0["lwd"], meas_hr0["swd"])]

    # Return numpy arrays
    for key in keys:
        meas_hr0[key] = np.array(meas_hr0[key])
    meas_hr0["tot"] = np.array(meas_hr0["tot"])

    return meas_hr0


# Directories and variables needed from other param files
lwd_base_dir = lwd_params.LWD_DIR
swd_base_dir = swd_params.SWD_DIR
lwd_fmt = lwd_params.LWD_FILEFORMAT
swd_fmt = swd_params.SWD_FILEFORMAT
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
era5fmt = "era5_esc_broadband_*[0-9]*.nc"
era_fmt = "era5_esc_broadband_*[0-9]*.nc"
era_lat = -62.25
era_lon = -59.0
file1 = f"{SAVEDIR}broadband_means.csv"
file2 = f"{SAVEDIR}broadband_diffs.csv"


# Directories
esc_dir = params.MEAS_DIR + "Escudero/"
era_dir = f"{esc_dir}/era5/broadband/"
dir_pwrf = params.MEAS_DIR + "Escudero/pwrf/broadband/"
frei_dir = f"{esc_dir}frei_met/daily_ave_temp/"

# Formats
fmt_pwrf_lwd = "wrfout_LWD_d03_%Y%m%d.nc"
fmt_pwrf_swd = "wrfout_SWD_d03_%Y%m%d.nc"
# temp_fmt = "esc_met_%Y%m.txt"
temp_fmt = "TemperaturaMediaDiariaClimatológica_%Y-%m.xlsx"
temp_fmt2 = "esc_met_%Y%m_new.txt"


# Get all polar wrf files
files_pwrf = os.listdir(dir_pwrf)
files_pwrf_lwd = [x for x in files_pwrf if "wrfout_LWD_d03_" in x]
files_pwrf_swd = [x for x in files_pwrf if "wrfout_SWD_d03_" in x]
files_pwrf_lwd.sort()
files_pwrf_swd.sort()

if len(files_pwrf_lwd) != len(files_pwrf_lwd):
    raise ValueError("Different number of LWD and SWD files")

pwrf_dates = [dt.datetime.strptime(x, fmt_pwrf_lwd) for x in files_pwrf_lwd]


# Parameters
redo_output_files = True
make_figs = True
make_qc_figs = False
years = range(2017, 2024)
pts_min = 200  # min hours in month to compute stats (of 672 to 744)

plot_forcings = True
season = "winter"  # "summer"

summer_months = [1, 2, 12]
winter_months = [5, 6, 7, 8]

date1_long = dt.datetime(2017, 1, 1, 0, 0)
date2_long = dt.datetime(2023, 8, 31, 23, 59)
# date1_long = dt.datetime(2021, 9, 1, 0, 0)
# date2_long = dt.datetime(2023, 8, 31, 23, 59)

keys = ["date", "lwd", "swd", "lwd_pts", "swd_pts"]


# Get surface temperatures, measurements and ERA5 for the month
# for all years and for the event
print("Looping over months")
era_2c = [[] for i in range(12)]
era = [[] for i in range(12)]
meas_hr = [[] for i in range(12)]
meas_hr_2c = [[] for i in range(12)]
temp_mon = [[] for i in range(12)]
dtm_temp_mon = [[] for i in range(12)]
if season == "summer":
    months = summer_months
else:
    months = winter_months
for mon in months:
    imon = mon - 1

    # Get surface temperatures at Frei
    temp_mon[imon] = []
    dtm_temp_mon[imon] = []
    for year in years:
        nday = calendar.monthrange(year, mon)[1]
        event_start = dt.datetime(year, mon, 1, 1)
        event_hours = nday * 24
        print(f"Getting Frei temps for {event_start} ", end="")
        print(f"for {event_hours/24} days: ", end="")
        temp0, date0 = get_frei_temps(event_start, event_hours)
        temp_mon[imon] += temp0
        dtm_temp_mon[imon] += date0

    meas0 = get_hourly_meas_over_years(
        years, mon, lwd_base_dir, lwd_fmt, swd_base_dir, swd_fmt, keys
    )
    meas0_2c = {}
    era0_2c, era0 = get_era5_bbnd_mon_warm(date1_long, date2_long, mon, 2)

    imeas = np.nonzero(np.isin(meas0["date"], era0["date"]))[0]
    date_meas = [meas0["date"][i] for i in imeas]
    imeas2c = np.nonzero(np.isin(meas0["date"], era0_2c["date"]))[0]
    date_meas2c = [meas0["date"][i] for i in imeas2c]
    iera2c = np.nonzero(np.isin(era0_2c["date"], date_meas2c))[0]
    # date_era = [era0['date'][i] for i in iera]
    # np.sum([(x - y).total_seconds() for x,y in zip(date_era, date_meas)])
    for key in ["lwd", "swd", "tot"]:
        meas0_2c[key] = meas0[key][imeas2c]
        meas0[key] = meas0[key][imeas]
        era0_2c[key] = era0_2c[key][iera2c]
        if plot_forcings:
            iera = np.nonzero(np.isin(era0["date"], date_meas))[0]
            meas0_2c[key] -= era0_2c[key + "_clr"][iera2c]
            meas0[key] -= era0[key + "_clr"][iera]
            era0[key] -= era0[key + "_clr"]
            era0_2c[key] -= era0_2c[key + "_clr"][iera2c]
    meas_hr[imon] = meas0
    meas_hr_2c[imon] = meas0_2c
    era[imon] = era0
    era_2c[imon] = era0_2c


# Get seasonal values
era_all = {}
era_2c_all = {}
meas_hr_all = {}
meas_2c_all = {}
for i, key in enumerate(["swd", "lwd", "tot", "date"]):
    era_all[key] = []
    era_all[key + "_clr"] = []
    meas_hr_all[key] = []
    meas_2c_all[key] = []
    era_2c_all[key] = []
    era_2c_all[key + "_clr"] = []
    for mon in months:
        meas_hr_all[key] += list(meas_hr[mon - 1][key])
        era_all[key] += list(era[mon - 1][key])
        era_2c_all[key] += list(era_2c[mon - 1][key])
        if key != "date":
            meas_2c_all[key] += list(meas_hr_2c[mon - 1][key])
            era_all[key + "_clr"] += list(era[mon - 1][key + "_clr"])
            era_2c_all[key + "_clr"] += list(era_2c[mon - 1][key + "_clr"])


temp_all = []
dtm_temp_all = []
for mon in months:
    temp_all += temp_mon[mon - 1]
    dtm_temp_all += dtm_temp_mon[mon - 1]


plt.figure()
plt.plot(dtm_temp_all, temp_all, ".")
plt.plot(meas_hr_all["date"], meas_hr_all["lwd"], ".")


# Match up temperature and radiation measurements
# First fast forward to first match
j = 0
while dtm_temp_all[j] < meas_hr_all["date"][0]:
    j += 1

nave_day = 0
meas_match_day = {"lwd": 0, "swd": 0, "tot": 0}
meas_match = {"lwd": [], "swd": [], "tot": []}
temp_match = []
date_match = []
for i, dtm in enumerate(meas_hr_all["date"]):
    meas_day = dt.datetime(dtm.year, dtm.month, dtm.day, 0)
    if meas_day == dtm_temp_all[j]:
        # We have a match
        nave_day += 1
        for key in ["lwd", "swd", "tot"]:
            meas_match_day[key] += meas_hr_all[key][i]
    elif meas_day > dtm_temp_all[j]:
        # We have passed the matching day, finish it
        for key in ["lwd", "swd", "tot"]:
            meas_match[key].append(meas_match_day[key] / nave_day)
            meas_match_day[key] = 0

        # and prepare for next
        nave_day = 0
        date_match.append(meas_day)
        temp_match.append(temp_all[j])
        j += 1

        if meas_day > dtm_temp_all[j]:
            # Fast forward to the matching day
            while dtm_temp_all[j] < meas_day:
                j += 1

        if meas_day == dtm_temp_all[j]:
            # We have a match
            nave_day += 1
            for key in ["lwd", "swd", "tot"]:
                meas_match_day[key] += meas_hr_all[key][i]
        else:
            raise ValueError("problem")


# There should not be any nans
if np.any(np.isnan(temp_all)):
    raise ValueError("There should not be any NaNs in the temperatures!")

# For winter, there is a string in temp_all - remove it
for i, t0 in enumerate(temp_all):
    if type(t0) == str:
        temp_all[i] = np.nan

# Get cases with surface temperature over 2C
temp_2c = [x for x in temp_all if x > 2]

# QC
if make_qc_figs:
    imon = 0
    fig, (axa, axb, axc) = plt.subplots(3, 1, num=1, clear=True)
    axa.plot(
        meas_hr[imon]["date"], meas_hr[imon]["lwd"], ".", label="meas lwd"
    )
    axa.plot(era[imon]["date"], era[imon]["lwd"], "+", label="ERA5 lwd")
    axb.plot(meas_hr[imon]["date"], meas_hr[imon]["swd"], "+", label="swd")
    axb.plot(era[imon]["date"], era[imon]["swd"], "+", label="ERA5 swd")
    axc.plot(meas_hr[imon]["date"], meas_hr[imon]["tot"], "s", label="tot")
    axc.plot(era[imon]["date"], era[imon]["tot"], "+", label="ERA5 tot")
    axc.legend()


if season == "summer":
    lft1 = 0.1
    lft2 = 0.77
    wid1 = 0.64
    wid2 = 0.2
    height_lw = 0.11
    height_sw = 0.39
    height_tot = 0.39
    bot3 = 0.08
    bot2 = 0.48
    bot1 = 0.6
    fignum = 1
else:
    lft1 = 0.1
    lft2 = 0.77
    wid1 = 0.64
    wid2 = 0.2
    height_lw = 0.33 / 420 * 220
    height_sw = 0.33
    height_tot = 0.33
    bot3 = 0.08
    bot2 = 0.44
    bot1 = 0.64
    fignum = 2


fig, axs = plt.subplots(4, 2, figsize=(6, 6), num=2, clear=True)

moff = 0
mall_off = 1  # 1.9  # offset for measurement
m2coff = 2  # 2.6  # offset for ERA5 box-and-whiskers when surface T above 2C
# eoff = 2.1  # offset for ERA5 box-and-whiskers
# e2coff = 2.8  # offset for ERA5 box-and-whiskers when surface T above 2C

ax1 = axs[0, 0]
ax2 = axs[1, 0]
ax3 = axs[2, 0]
ax4 = axs[0, 1]
ax5 = axs[1, 1]
ax6 = axs[2, 1]
ax7 = axs[3, 0]
ax8 = axs[3, 1]

if plot_forcings:
    lft1 = 0.1
    lft2 = 0.81
    wid1 = 0.69
    wid2 = 0.18
    height_sw = 0.23
    height_lw = 0.20
    height_tot = 0.22
    height_temps = 0.22
    bot1 = 0.76  # 0.68
    bot2 = 0.55  # 0.39
    bot3 = 0.32  # 0.08
    bot4 = 0.09
    ylab = "Forcing (W/m$^2$)"
else:
    # lft1 = 0.1
    # lft2 = 0.81
    # wid1 = 0.69
    # wid2 = 0.18
    # height_lw = 0.11
    # height_sw = 0.39
    # height_tot = 0.39
    # bot3 = 0.08
    # bot2 = 0.48
    # bot1 = 0.6
    lft1 = 0.1
    lft2 = 0.81
    wid1 = 0.69
    wid2 = 0.18
    height_sw = 0.3 * 3 / 4
    height_lw = 0.26 * 3 / 4
    height_tot = 0.28 * 3 / 4
    height_temps = 0.28 * 3 / 4
    bot1 = 0.72  # 0.68
    bot2 = 0.53  # 0.39
    bot3 = 0.33  # 0.08
    bot4 = 0.08
    ylab = "Flux (W/m$^2$)"

ax1.set_position([lft1, bot1, wid1, height_sw])
ax2.set_position([lft1, bot2, wid1, height_lw])
ax3.set_position([lft1, bot3, wid1, height_tot])
ax4.set_position([lft2, bot1, wid2, height_sw])
ax5.set_position([lft2, bot2, wid2, height_lw])
ax6.set_position([lft2, bot3, wid2, height_tot])
ax7.set_position([lft1, bot4, wid1, height_tot])
ax8.set_position([lft2, bot4, wid2, height_tot])

date_form = DateFormatter("%m/%d")
orng = "orange"
mcolor = "green"  # "darkturquoise"  # measurement color
acolor = "blue"
wcolor = "orange"

# Plot box-and-whiskers for the seasonal values for all cases and cases > 2c
for i, key in enumerate(["swd", "lwd", "tot"]):
    axb = axs[i][1]

    axb.boxplot(
        meas_hr_all[key],
        whis=(0, 100),
        positions=[mall_off],
        boxprops=dict(color=acolor),
        medianprops=dict(color=acolor),
        whiskerprops=dict(color=acolor),
        capprops=dict(color=acolor),
    )
    axb.plot(mall_off, np.nanmean(meas_hr_all[key]), "*", color=acolor)

    # Boxplots for surface temp > 2C according to ERA5
    # (looks very similar to results based on Frei temps)
    # axb.boxplot(
    #     meas_2c_all[key],
    #     whis=(0, 100),
    #     positions=[m2coff],
    #     boxprops=dict(color=wcolor),
    #     medianprops=dict(color=wcolor),
    #     whiskerprops=dict(color=wcolor),
    #     capprops=dict(color=wcolor),
    # )
    # axb.plot(m2coff, np.nanmean(meas_2c_all[key]), "*", color=wcolor)

# Add on the temperatures
temp_all = [x for x in temp_all if np.isnan(x) == False]
i += 1
axs[i][1].boxplot(
    temp_all,
    whis=(0, 100),
    positions=[mall_off],
    boxprops=dict(color=wcolor),
    medianprops=dict(color=wcolor),
    whiskerprops=dict(color=wcolor),
    capprops=dict(color=wcolor),
)
axs[i][1].plot(mall_off, np.nanmean(temp_all), "*", color=wcolor)

axs[i][1].boxplot(
    temp_2c,
    whis=(0, 100),
    positions=[m2coff],
    boxprops=dict(color=wcolor),
    medianprops=dict(color=wcolor),
    whiskerprops=dict(color=wcolor),
    capprops=dict(color=wcolor),
)
axs[i][1].plot(m2coff, np.nanmean(temp_2c), "*", color=wcolor)

# axb.boxplot(
#     era_all[key],
#     whis=(0, 100),
#     positions=[eoff],
#     boxprops=dict(color="k"),
#     medianprops=dict(color="k"),
# )
# axb.plot(eoff, np.nanmean(era_all[key]), "k*")

# axb.boxplot(
#     era_2c_all[key],
#     whis=(0, 100),
#     positions=[e2coff],
#     boxprops=dict(color="k"),
#     medianprops=dict(color="k"),
# )
# axb.plot(e2coff, np.nanmean(era_2c_all[key]), "k*")


# Loop over PWRF dates and get starting and ending date indices of AR cases
icase_beg = [0]
icase_end = []
for i in range(1, len(pwrf_dates)):
    if (pwrf_dates[i] - pwrf_dates[i - 1]).days > 1:
        icase_end.append(i - 1)
        icase_beg.append(i)
# icase_beg.pop()
icase_end.append(len(pwrf_dates) - 1)

# Choose subset for selected months
icase_beg0 = []
icase_end0 = []
for ibeg, iend in zip(icase_beg, icase_end):
    if pwrf_dates[ibeg].month in months:
        icase_beg0.append(ibeg)
        icase_end0.append(iend)

# Sort them by mean surface temperaure
if season == "summer":
    # fmt: off
    #         0  1 12/8 2/20  1/19  2/7  2/4  2/8  2/27  2/4 2/21
    iorder = [6, 8, 5,   10,  7,     9,   1,   3,  4,     2,  0]
    # fmt: on
else:
    iorder = [4, 3, 2, 1, 0]
icase_beg = [icase_beg0[i] for i in iorder]
icase_end = [icase_end0[i] for i in iorder]

# # Sort them by month/day
# pwrf_event_dates = [pwrf_dates[i] for i in icase_beg0]
# event_days = [x.timetuple().tm_yday for x in pwrf_event_dates]
# isort = np.argsort(event_days)
# icase_beg = [icase_beg0[i] for i in isort]
# icase_end = [icase_end0[i] for i in isort]

# if 12 in months:
#     # Put December first
#     idec = [
#         i
#         for i in range(len(icase_beg))
#         if pwrf_dates[icase_beg[i]].month == 12
#     ]
#     idec = idec[0]
#     icase_beg = icase_beg[idec:] + icase_beg[:idec]
#     icase_end = icase_end[idec:] + icase_end[:idec]


era_events = {}
pwrf_events = {}
meas_events = {}
for key in ["swd", "lwd", "tot"]:
    era_events[key] = []
    pwrf_events[key] = []
    meas_events[key] = []
    era_events[key + "_clr"] = []

# TOOD: Loop over AR events
# Loop over AR events
print("\nLooping over AR events")
temp_mean = []
meas_mean = {}
frei_temps_events = []
for key in ["swd", "lwd", "tot"]:
    meas_mean[key] = []
count = 0
events = []
fmt = "%Y/%m/%d"
for ibeg, iend in zip(icase_beg, icase_end):
    print(
        f"\n{pwrf_dates[ibeg].strftime(fmt)} - {pwrf_dates[iend].strftime(fmt)}"
    )

    event = pwrf_dates[ibeg].strftime("%Y%m%d")
    events.append(f"{event[:4]}/{event[4:6]}/{event[6:8]}")
    case_year = pwrf_dates[ibeg].year
    case_mon = pwrf_dates[ibeg].month
    case_day = pwrf_dates[ibeg].day
    case_mon_end = pwrf_dates[iend].month
    case_day_end = pwrf_dates[iend].day
    duration = (pwrf_dates[iend] - pwrf_dates[ibeg]).days + 1

    # Parameters
    if event == "20220207":
        # This uses data in a different folder, from prior PWRF runs
        event = "20220206"
        dir_pwrf_20220206 = f"{esc_dir}case_studies/case_202202/"
        pwrf_file = "PWRF_D03_Escudero_20220206_09_LWSW_hrly.nc"
        case_year = 2022
        case_mon = 2
        duration = 1 + 21 / 24
        date1 = dt.datetime(case_year, case_mon, 6, 0, 0)
        date2 = dt.datetime(case_year, case_mon, 9, 23, 0)
        pwrf_date1 = dt.datetime(case_year, case_mon, 6, 0, 0)
        pwrf_date2 = dt.datetime(case_year, case_mon, 9, 23, 0)
        date1_utc = esc202202.DATE1 - dt.timedelta(days=5)
        date2_utc = esc202202.DATE2
        plot_date1 = dt.datetime(case_year, case_mon, 4)
        plot_date2 = dt.datetime(case_year, case_mon, 12)

        event_start = dt.datetime(case_year, 2, 6, 0, 30, 30)
        # Start date for computing seconds for interpolating to created plots
        start_date_for_s = dt.datetime(case_year, 1, 1)

    else:
        date1 = dt.datetime(case_year, case_mon, case_day, 0, 0)
        date2 = dt.datetime(case_year, case_mon_end, case_day_end, 23, 59)
        pwrf_date1 = date1
        pwrf_date2 = date2
        date1_utc = date1.replace(tzinfo=pytz.utc) - dt.timedelta(days=2)
        date2_utc = date2.replace(tzinfo=pytz.utc) + dt.timedelta(days=2)
        plot_date1 = date1_utc.replace(tzinfo=None)
        plot_date2 = date2_utc.replace(tzinfo=None)
        event_start = dt.datetime(case_year, case_mon, case_day, 0, 30, 30)
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    imon = case_mon - 1
    print(f"case date start {date1}")
    print(f"case date end   {date2}")
    print(f"case date duration {duration}")

    # Directories
    era_dir = f"{esc_dir}/era5/broadband/"
    dir_swd = f"{esc_dir}pyranometer/stand/{case_year}/"
    dir_lwd = f"{esc_dir}pyrgeometer/stand/{case_year}/"

    # Date range for data for the month for all years for box-and_whiskers
    # The start date for computing seconds for interpolating to created plots
    start_date_for_s = dt.datetime(case_year, case_mon, 1)
    event_hours = int(duration * 24)

    frei_temps, _ = get_frei_temps(event_start, event_hours)

    # Load in the polar wrf data
    if event == "20220206":
        with Dataset(dir_pwrf_20220206 + pwrf_file) as nci:
            ntimes = len(nci["LWD_tosta"])
            pwrf = {
                "swd": nci["SWD_tosta"][:].data,
                "lwd": nci["LWD_tosta"][:].data,
            }
        dates = np.linspace(pwrf_date1.day, pwrf_date2.day + 23 / 24, ntimes)
        pwrf_dt = []
        for i, d in enumerate(dates):
            h = int(d % int(d) * 24)
            pwrf_dt.append(
                dt.datetime(pwrf_date1.year, pwrf_date1.month, int(d), h)
            )
        pwrf["date"] = pwrf_dt
        pwrf["tot"] = pwrf["swd"] + pwrf["lwd"]
    else:
        pwrf = get_pwrf_bbnd_daterange(pwrf_date1, pwrf_date2)

    # Get ERA5 data for this case (downloaded by get_era5_from_web):
    direc = era_dir + str(case_year) + "/"
    era5 = get_era5_atgridpts(
        direc, era_fmt, date1_utc, date2_utc, era_lat, era_lon
    )
    era5["date"] -= dt.timedelta(minutes=30)

    if plot_forcings:
        meas_event, era_event, pwrf_event = get_event_forcing(
            date1_utc,
            date2_utc,
            dir_swd,
            dir_lwd,
            event_start,
            event_hours,
            era5,
            pwrf,
        )
    else:
        meas_event, era_event, pwrf_event = get_event_flux(
            date1_utc,
            date2_utc,
            dir_swd,
            dir_lwd,
            event_start,
            event_hours,
            era5,
            pwrf,
        )

    if make_figs:
        # moff = 0.0
        # eoff = 0.2
        # poff = 0.4
        # Add box-and-whiskers for the AR event
        for i, key in enumerate(["swd", "lwd", "tot"]):
            axa = axs[i][0]

            pos = count
            # fmt: off
            if count == 0 and i == 0:
                axa.plot(pos + moff, np.nanmean(meas_event[key]), "*", color=mcolor, label="Meas")

            # axa.boxplot(era_event[key], whis=(0, 100), positions=[pos + eoff],
            #     boxprops=dict(color="k"),
            #     medianprops=dict(color="k"),
            #     whiskerprops=dict(color="k"),
            #     capprops=dict(color="k"),
            # )            

            axa.boxplot(meas_event[key], whis=(0, 100), positions=[pos + moff],
                boxprops=dict(color=mcolor),
                medianprops=dict(color=mcolor),
                whiskerprops=dict(color=mcolor),
                capprops=dict(color=mcolor),
            )

            axa.plot(pos + moff, np.nanmean(meas_event[key]), "*", color=mcolor)
            # axa.plot(pos + eoff, np.nanmean(era_event[key]), "k*")
            # fmt: on

            # combine up the events
            meas_events[key] += list(meas_event[key])
            # era_events[key + "_clr"] += list(era_event[key + "_clr"])

            meas_mean[key].append(np.nanmean(meas_event[key]))

        i += 1
        # Temperatures
        # Remove nans, but warn
        if season == "summer":
            yasterisk = 9
        else:
            yasterisk = 4
        if np.any(np.isnan(frei_temps)):
            raise ValueError(f"Nans found for AR beginning {date1}")
            frei_temps = [x for x in frei_temps if np.isnan(x) == False]
            axs[i][0].plot(pos + moff, yasterisk, "k*")

        axs[i][0].boxplot(
            frei_temps,
            whis=(0, 100),
            positions=[pos + moff],
            boxprops=dict(color=mcolor),
            medianprops=dict(color=mcolor),
            whiskerprops=dict(color=mcolor),
            capprops=dict(color=mcolor),
        )
        axs[i][0].plot(pos + moff, np.nanmean(frei_temps), "*", color=mcolor)
        temp_mean.append(np.nanmean(frei_temps))
        frei_temps_events += frei_temps

        count += 1


# Add the averages for all events
# if 12 in months:
#     # moff = 1.0
#     # eoff = 1.2
#     # poff = 1.4
# else:
#     # moff = 1.0
#     # eoff = 1.0
#     # poff = 1.2
for i, key in enumerate(["swd", "lwd", "tot"]):
    axa = axs[i][1]
    axa.boxplot(
        meas_events[key],
        whis=(0, 100),
        positions=[moff],
        boxprops=dict(color=mcolor),
        medianprops=dict(color=mcolor),
        whiskerprops=dict(color=mcolor),
        capprops=dict(color=mcolor),
    )

    # axa.boxplot(
    #     era_events[key],
    #     whis=(0, 100),
    #     positions=[eoff],
    #     boxprops=dict(color="k"),
    #     medianprops=dict(color="k"),
    #     whiskerprops=dict(color="k"),
    #     capprops=dict(color="k"),
    # )

    axa.plot(moff, np.nanmean(meas_events[key]), "*", color=mcolor)
    # axa.plot(eoff, np.nanmean(era_events[key]), "k*")

i += 1
axa = axs[i][1]
axa.boxplot(
    frei_temps_events,
    whis=(0, 100),
    positions=[moff],
    boxprops=dict(color=mcolor),
    medianprops=dict(color=mcolor),
    whiskerprops=dict(color=mcolor),
    capprops=dict(color=mcolor),
)
axa.plot(moff, np.nanmean(frei_temps_events), "*", color=mcolor)

# Add labels
# Set the ylimits to be the same as the panel to the left
# axb.set_xticks([1, 2, 3])
# axb.set_xticks([1, 2, 3])
# axb.set_xticklabels(["Event", ">2C", "all"], rotation=10)


eventlabels = [
    x[5:7].lstrip("0") + "/" + x[8:10].lstrip("0") + "\n" + x[:4]
    for x in events
]
axs[3][0].set_xticks(np.arange(len(events)))  #  + 0.2)
axs[3][0].set_xticklabels(eventlabels, rotation=0)
axs[3][0].set_xlabel("Date of AR Case (day/month year)")

for i in range(4):
    axs[i][1].set_xlim([-0.4, 2.3])
    axs[i][1].set_xticks([0, 1, 2])
axs[3][1].set_xticklabels(["Events\nall", "DJF\nall", "DJF\nT>2C"])

axs[0][0].set_xticklabels([])
axs[1][0].set_xticklabels([])
axs[0][1].set_xticklabels([])
axs[1][1].set_xticklabels([])
axs[2][0].set_xticklabels([])
axs[2][1].set_xticklabels([])

# Comment out to see y tick labels on panels to right
axs[0][1].set_yticklabels([])
axs[1][1].set_yticklabels([])
axs[2][1].set_yticklabels([])
axs[3][1].set_yticklabels([])

axs[0][0].set_ylabel(ylab)  # "SWD (W/m$^2$)")
axs[1][0].set_ylabel(ylab)  # ("LWD (W/m$^2$)")
axs[2][0].set_ylabel(ylab)  # ("LWD+SWD (W/m$^2$)")
axs[3][0].set_ylabel("Temperature (C)")

if plot_forcings:
    if season == "summer":
        for i in range(2):
            axs[0][i].set_ylim([-820, 300])
            axs[1][i].set_ylim([-10, 120])
            axs[2][i].set_ylim([-750, 499])
            axs[3][i].set_ylim([-6, 12])
        axs[0][0].text(-0.4, 190, "a) DSW")
        axs[0][1].text(-0.3, 190, "b) DSW")
        axs[1][0].text(-0.4, 100, "c) DLW")
        axs[1][1].text(-0.3, 105, "d) DLW")
        axs[2][0].text(-0.4, 370, "e) DSW + DLW")
        axs[2][1].text(-0.3, 370, "f) DSW + DLW")
        axs[3][0].text(-0.4, 10, "g) Temperature")
        axs[3][1].text(-0.3, 10, "h) Temperature")
        # axs[0][0].legend(
        #     ["Meas", "ERA5"], loc="lower right", ncol=3, borderpad=0.4
        # )
        fig_file = f"{SAVEDIR}esc_forcing_cases_summerb.png"
    else:
        for i in range(2):
            axs[0][i].set_ylim([-300, 120])
            axs[1][i].set_ylim([-20, 120])
            axs[2][i].set_ylim([-300, 200])
            axs[3][i].set_ylim([-21, 8.5])
        axs[0][0].text(-0.4, 80, "a) DSW")
        axs[1][0].text(-0.4, 103, "c) DLW")
        axs[2][0].text(-0.4, 150, "e) DSW + DLW")
        axs[0][1].text(-0.3, 80, "b) DSW")
        axs[1][1].text(-0.3, 103, "d) DLW")
        axs[2][1].text(-0.3, 150, "f) DSW + DLW")
        axs[3][0].text(-0.4, 6, "g) Temperature")
        axs[3][1].text(-0.3, 6, "h) Temperature")
        # axs[0][0].legend(["Meas", "ERA5"], loc="lower center", borderpad=0.4)
        fig_file = f"{SAVEDIR}esc_forcing_cases_winterb.png"
else:
    if season == "summer":
        for i in range(2):
            axs[0][i].set_ylim([-50, 1100])
            axs[1][i].set_ylim([190, 380])
            axs[2][i].set_ylim([150, 1500])
            axs[3][i].set_ylim([-6, 10.5])

        axs[0][0].text(-0.4, 1000, "a) DSW")
        axs[0][1].text(-0.3, 1000, "b) DSW")
        axs[1][0].text(-0.4, 358, "c) DLW")
        axs[1][1].text(-0.3, 358, "d) DLW")
        axs[2][0].text(-0.4, 1350, "e) DSW + DLW")
        axs[2][1].text(-0.3, 1350, "f) DSW + DLW")

        axs[3][0].text(-0.4, 9, "g) Temperature")
        axs[3][1].text(-0.3, 9, "h) Temperature")
        fig_file = f"{SAVEDIR}esc_broadband_cases_summer.png"
    else:
        for i in range(2):
            axs[0][i].set_ylim([-20, 400])
            axs[1][i].set_ylim([130, 350])
            axs[2][i].set_ylim([130, 550])
            axs[3][i].set_ylim([-2.5, 10.5])
        axs[0][0].text(4.6, 350, "a)")
        axs[1][0].text(4.6, 280, "c)")
        axs[2][0].text(4.6, 500, "e)")
        axs[0][1].text(0.4, 350, "b)")
        axs[1][1].text(0.4, 280, "d)")
        axs[2][1].text(0.4, 500, "f)")
        axs[3][0].text(-0.4, 9, "g) Temperature")
        axs[3][1].text(-0.3, 9, "h) Temperature")

        fig_file = f"{SAVEDIR}esc_broadband_cases_winter.png"
    # axs[0][0].legend(["Meas", "ERA5", "PWRF"])


# Add the matched measurements to the plot
for i, key in enumerate(["swd", "lwd", "tot"]):
    axb = axs[i][1]

    axb.boxplot(
        meas_match[key],
        whis=(0, 100),
        positions=[m2coff],
        boxprops=dict(color=wcolor),
        medianprops=dict(color=wcolor),
        whiskerprops=dict(color=wcolor),
        capprops=dict(color=wcolor),
    )
    axb.plot(m2coff, np.nanmean(meas_match[key]), "*", color=wcolor)

if SAVEFIGS:
    if make_qc_figs:
        raise ValueError("Set SAVEFIGS to False for this option!")
    plt.savefig(fig_file, format="png", dpi=600)


temp_match_day = temp_match
meas_match_day = meas_match


# fmt: off
panel_label = ['a) DSW', 'b) DLW', 'c) DSW + DLW']
if season == 'summer':
    fig, ax = plt.subplots(3, 1, num=3, figsize=(7, 5.8), clear=True)

    ytxt1 = [-160, 60, -100]
    ytxt2 = [-220, 48, -150]
    ytxt3 = [-620, 10, -620]
    ytxt4 = [-50, 80, -5]
    temp_mean = np.array(temp_mean)
    temp_match_day = np.array(temp_match_day)

    for i, key in enumerate(["swd", "lwd", "tot"]):
        # temp_match_day = np.array(temp_match_day)
        meas_mean[key] = np.array(meas_mean[key])
        meas_match_day[key] = np.array(meas_match_day[key])
        itemp = np.nonzero(temp_match_day > -1)[0]
        itemp = np.intersect1d(itemp, np.nonzero(temp_match_day < 3.5)[0])
        pfit = np.polyfit(temp_match_day[itemp], meas_match_day[key][itemp], 1)
        yall = np.polyval(pfit, temp_match_day[itemp])
        m, b, r, p, se = linregress(temp_match_day, meas_match_day[key])
        if i==1:
            ax[i].text(-4, ytxt1[i], f"m= {m:.0f}; r={r:.2f}; p={p:.2f}; se={se:.1f}", color="b")
        else:
            ax[i].text(-4, ytxt1[i], f"m= {m:.0f}; r={r:.2f}; p={p:.1f}; se={se:.0f}", color="b")
    
        idx = np.isfinite(temp_mean) & np.isfinite(meas_mean[key])
        p = np.polyfit(temp_mean[idx], meas_mean[key][idx], 1)
        xevent = [np.nanmin(temp_mean), np.nanmax(temp_mean)]
        yevent = np.polyval(p, xevent)
        m, b, r, p, se = linregress(temp_mean[idx], meas_mean[key][idx])
        if i==1:
            ax[i].text(-4, ytxt2[i], f"m={m:.0f}; r={r:.1f}; p={p:.1f}; se={se:.0f}", color="darkgreen")
        else:
            ax[i].text(-4, ytxt2[i], f"m={m:.0f}; r={r:.1f}; p={p:.2f}; se={se:.0f}", color="darkgreen")
        ax[i].text(-4, ytxt4[i], panel_label[i])
    
        # ax[i].plot(-100,0, "b.", label=lab)
        lab1 = f"Summer all: {np.nanmean(meas_hr_all[key]):.0f}"
        lab2 = f"Summer ARs: {np.nanmean(meas_events[key]):.0f}"
        ax[i].plot(temp_match_day, meas_match_day[key], "b.", alpha=0.2, label=lab1)
        ax[i].plot(temp_match_day[itemp], yall, "b", alpha=0.5)
        ax[i].plot(temp_mean, meas_mean[key], "go", label=lab2)
        ax[i].plot(xevent, yevent, "g") #, label='Best fit')
        ax[i].legend(loc='lower left')    
        ax[i].set_xlim([-4.2, 5])

    fig_file_corr = f"{SAVEDIR}esc_forcing_corr_summer.png"
        
else:
    lab = 'Fall/Winter, all'
    ytxt1 = [-28, 75, 17]
    ytxt2 = [-32, 65, 7]
    ytxt3 = [-15, 30, 30]
    ytxt4 = [0, 85, 80]
    temp_mean = np.array(temp_mean)
    fig, ax = plt.subplots(3, 1, num=4, figsize=(7, 6), clear=True)
    for i, key in enumerate(["swd", "lwd", "tot"]):
        meas_mean[key] = np.array(meas_mean[key])
        pfit = np.polyfit(temp_match, meas_match[key], 1)
        yall = np.polyval(pfit, temp_match_day)
        m, b, r, p, se = linregress(temp_match, meas_match[key])
        rstr = f"{r:.2f}".lstrip('0')
        pstr = f"{p:.9f}".lstrip('0')
        if i==0:
            ax[i].text(-13.5, ytxt1[i], f"m={m:.1f}; r={r:.1f}; p<.01; se={se:.2f}", color="b")
        else:
            ax[i].text(-13.5, ytxt1[i], f"m={m:.1f}; r={r:.1f}; p<.01; se={se:.1f}", color="b")

        idx = np.isfinite(temp_mean) & np.isfinite(meas_mean[key])
        p = np.polyfit(temp_mean[idx], meas_mean[key][idx], 1)
        xevent = [np.nanmin(temp_mean), np.nanmax(temp_mean)]
        yevent = np.polyval(p, xevent)
        m, b, r, p, se = linregress(temp_mean[idx], meas_mean[key][idx])
        rstr = f"{r:.1f}".lstrip('0')
        pstr = f"{p:.1f}".lstrip('0')        
        ax[i].text(-13.5, ytxt2[i], f"m={m:.0f}; r={r:.1f}; p={p:.1f}; se={se:.0f}", color="darkgreen")            
        ax[i].text(-13.5, ytxt4[i], panel_label[i])
    
        # ax[i].plot(-100,0, "b.", label=lab)
        lab1 = f"Fall/winter all: {np.nanmean(meas_hr_all[key]):.0f}"
        lab2 = f"Fall/winter ARs: {np.nanmean(meas_events[key]):.0f}"
        ax[i].plot(temp_match_day, meas_match_day[key], "b.", alpha=0.2, label=lab1)
        ax[i].plot(temp_mean, meas_mean[key], "go", label=lab2)
        ax[i].plot(temp_match_day, yall, "b", alpha=0.5, label='Best fit, all')
        ax[i].plot(xevent, yevent, "g", label='Best fit, ARs')
        ax[i].set_xlim([-14, 2])
        if i == 1:
            ax[i].legend(loc='lower left')
            ax[i].set_ylim([0, 95])

        fig_file_corr = f"{SAVEDIR}esc_forcing_corr_winter.png"

    

for i, key in enumerate(["swd", "lwd", "tot"]):
    ax[i].set_ylabel(ylab)
ax[i].set_xlabel("Surface Temperature ($\degree$C)")
# fmt: on


# labelx = -0.3  # axes coords
# for j in range(2):
#     ax[j].yaxis.set_label_coords(labelx, 0.5)
fig.align_ylabels(ax)
fig.tight_layout()


if SAVEFIGS:
    if make_qc_figs:
        raise ValueError("Set SAVEFIGS to False for this option!")
    plt.savefig(fig_file_corr, format="png", dpi=600)
    plt.savefig(fig_file_corr.rstrip("png") + "eps", format="eps")


# Print out the results
print("\n\n")
print("Means for all events: ")
print("DSW ERA5/PWRF/Meas, DLW ERA5/PWRF/Meas, Tot DSW ERA5/PWRF/Meas")
meas_mean_events = dict()
for key in ["swd", "lwd", "tot"]:
    print(f"{np.nanmean(era_events[key]):.0f}", end="  ")
    print(f"{np.nanmean(pwrf_events[key]):.0f}", end="  ")
    print(f"{np.nanmean(meas_events[key]):.0f}", end="    ")
    meas_mean_events[key] = np.nanmean(meas_events[key])
print()


meas_mean_all = dict()
print("\nMeans for season: DSW ERA5/Meas, DLW ERA5/Meas, Tot DSW ERA5/Meas")
for i, key in enumerate(["swd", "lwd", "tot"]):
    print(f"{np.nanmean(era_all[key]):.1f}", end="\t")
    print(f"{np.nanmean(meas_hr_all[key]):.1f}", end="\t")
    meas_mean_all[key] = np.nanmean(meas_hr_all[key])
print()

print("\nMeans for season for surface T > 2C")
for i, key in enumerate(["swd", "lwd", "tot"]):
    print(f"{np.nanmean(era_2c_all[key]):.0f}", end="\t")
    print(f"{np.nanmean(meas_2c_all[key]):.0f}", end="\t")


print("\nERA5 event clear sky means: DSW/DLW/tot")
for i, key in enumerate(["swd", "lwd", "tot"]):
    val = np.nanmean(era_events[key + "_clr"])
    print(f"{val:.1f}", end="\t")

print("\nERA5 season clear sky means: DSW/DLW/tot")
for i, key in enumerate(["swd", "lwd", "tot"]):
    val = np.nanmean(era_all[key + "_clr"])
    print(f"{val:.1f}", end="\t")

print()
