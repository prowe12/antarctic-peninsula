#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:21:08 2022

@author: prowe
"""


import datetime as dt
import numpy as np
import pandas as pd
import csv


from antarc.era5_read_station_broadband import read_era5_broadband_down_by_file
from antarc.escudero.all.results import get_era5_broadband_monthly_warm
from antarc.escudero.all.results import get_broadband_monthly_meas


get_broadband_monthly_meas = (
    get_broadband_monthly_meas.get_broadband_monthly_meas
)
get_era5_broadband_monthly_warm = (
    get_era5_broadband_monthly_warm.get_era5_broadband_monthly_warm
)


def get_era5_broadband_month_warm(start_date, end_date, month, tempthresh):
    era_warm, era_all = get_era5_broadband_monthly_warm(
        start_date, end_date, tempthresh
    )

    # Get just values for month
    imon = [i for i, x in enumerate(era_all["date"]) if x.month == month]
    for key in era_all.keys():
        era_all[key] = era_all[key][imon]

    imon = [i for i, x in enumerate(era_warm["date"]) if x.month == month]
    for key in era_warm.keys():
        era_warm[key] = era_warm[key][imon]

    # Add on the totals
    era_all["tot"] = era_all["swd"] + era_all["lwd"]
    era_all["tot_clr"] = era_all["swd_clr"] + era_all["lwd_clr"]
    era_warm["tot"] = era_warm["swd"] + era_warm["lwd"]
    era_warm["tot_clr"] = era_warm["swd_clr"] + era_warm["lwd_clr"]

    return era_warm, era_all


def get_broadband(fnames):
    rad = []
    mydates = []
    for filename in fnames:

        if (
            filename[filename.rfind("/") :] == "/esc_lwd20220209_2035.csv"
            or filename[filename.rfind("/") :] == "/esc_swd20220209_2033.csv"
        ):
            # Add nans in the missing region
            nanvec = np.array([np.nan, np.nan])
            dtimevec = np.array(
                [dt.datetime(2022, 2, 7, 21, 24), dt.datetime(2022, 2, 9, 20)]
            )
            rad = np.hstack([rad, nanvec])
            mydates = np.hstack([mydates, dtimevec])

        # Read in the csv file for ifile
        df = pd.read_csv(filename, delimiter=",")

        # Get dates, skipping the 0th element, which is the lower part of the header
        dates_str = df["Date"] + " " + df["Time"]

        # Get rad and samples as float
        rad = np.hstack([rad, np.array([float(x) for x in df["Radiation"]])])

        # Get the date
        dt0 = [dt.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in dates_str]
        mydates = np.hstack([mydates, dt0])

    # Get the hourly averages
    radsum = 0
    naves = 0
    try:
        hour = mydates[0].hour
    except:
        print("pause")
    dateave = []
    radave = []
    for i, mydate in enumerate(mydates):
        if mydate.hour == hour:
            radsum += rad[i]
            naves += 1
        else:
            dt0 = mydates[i - 1]
            thisdate = dt.datetime(
                dt0.year, dt0.month, dt0.day, dt0.hour, 30, 30
            )
            dateave.append(thisdate)
            radave.append(radsum / naves)
            hour = mydate.hour
            radsum = rad[i]
            naves = 1

    return mydates, rad, dateave, radave


def get_era5_atgridpts(direc, fmt, dt1, dt2, lat, lon):
    era5a = read_era5_broadband_down_by_file(direc, fmt, dt1, dt2)

    # Get the closest grid point
    i = np.where(era5a["lat"] == lat)[0][0]
    j = np.where(era5a["lon"] == lon)[0][0]

    return {
        "date": era5a["dates"],
        "lwd": era5a["lwd"][:, i, j],
        "swd": era5a["swd"][:, i, j],
        "tot": era5a["lwd"][:, i, j] + era5a["swd"][:, i, j],
        "lwd_clr": era5a["lwd_clr"][:, i, j],
        "swd_clr": era5a["swd_clr"][:, i, j],
        "tot_clr": era5a["lwd_clr"][:, i, j] + era5a["swd_clr"][:, i, j],
    }


def write_means_tofile(
    fname, case_date1, event_days, era_event, pwrf_event, meas_event
):
    """Write the results to files"""

    header1 = [
        "Date",
        "Duration (days)",
        "Mean SWD (W/m2)",
        "",
        "",
        "Mean LWD (W/m2)",
        "",
        "",
        "Mean Net down (W/m2)",
        "",
        "",
    ]
    header2 = [
        "",
        "",
        "ERA5",
        "PWRF",
        "Meas",
        "ERA5",
        "PWRF",
        "Meas",
        "ERA5",
        "PWRF",
        "Meas",
    ]

    data = [case_date1.strftime("%Y/%m/%d"), event_days]
    for key in ["swd", "lwd", "tot"]:
        data.append(f"{np.nanmean(era_event[key]):.0f}")
        data.append(f"{np.nanmean(pwrf_event[key]):.0f}")
        data.append(f"{np.nanmean(meas_event[key]):.0f}")

    with open(fname, "w", encoding="UTF8", newline="") as f:
        writer = csv.writer(f)

        # write the header
        writer.writerow(header1)
        writer.writerow(header2)

        # write the data
        writer.writerow(data)


def append_means_tofile(
    fname, case_date1, event_days, era_event, pwrf_event, meas_event
):
    """Append the results to files"""

    data = [case_date1.strftime("%Y/%m/%d"), event_days]
    for key in ["swd", "lwd", "tot"]:
        data.append(f"{np.nanmean(era_event[key]):.0f}")
        data.append(f"{np.nanmean(pwrf_event[key]):.0f}")
        data.append(f"{np.nanmean(meas_event[key]):.0f}")

    with open(fname, "a", encoding="UTF8", newline="") as f:
        # append the data
        writer = csv.writer(f)
        writer.writerow(data)


def write_diffs_tofile(fname, date1, days, era, pwrf, meas):
    """Write the results to files"""

    header1 = [
        "Date",
        "Duration (days)",
        "Mean SWD (W/m2)",
        "",
        "Mean LWD (W/m2)",
        "",
        "Mean Net Down (W/m2)",
        "",
    ]

    header2 = [
        "",
        "",
        "ERA5 - Meas",
        "PWRF - Meas",
        "ERA5 - Meas",
        "PWRF - Meas",
        "ERA5 - Meas",
        "PWRF - Meas",
    ]

    data = [date1.strftime("%Y/%m/%d"), days]

    with open(fname, "w", encoding="UTF8", newline="") as f:
        writer = csv.writer(f)

        # write the header
        writer.writerow(header1)
        writer.writerow(header2)

        # write the data
        writer.writerow(data)

    append_diffs_tofile(fname, date1, days, era, pwrf, meas)


def append_diffs_tofile(fname, date1, days, era, pwrf, meas):
    """Append the results to files"""

    data = [date1.strftime("%Y/%m/%d"), days]
    for key in ["swd", "lwd", "tot"]:
        ediff = np.nanmean(era[key]) - np.nanmean(meas[key])
        pdiff = np.nanmean(pwrf[key]) - np.nanmean(meas[key])
        data.append(f"{ediff:.1f}")
        data.append(f"{pdiff:.1f}")
        # if ediff < 2:
        #     data.append(f"{ediff:.1f}")
        # else:
        #     data.append(f"{ediff:.0f}")
        # if pdiff < 2:
        #     data.append(f"{pdiff:.1f}")
        # else:
        #     data.append(f"{pdiff:.0f}")

    with open(fname, "a", encoding="UTF8", newline="") as f:
        # append the data
        writer = csv.writer(f)
        writer.writerow(data)
