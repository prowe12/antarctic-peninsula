#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:13:40 2024

@author: prowe
"""

import datetime as dt
import matplotlib.pyplot as plt
import numpy as np


from antarc.getfilenames import getfilenames_daterange
from antarc.escudero.all.results import get_era5_broadband_monthly_warm
from antarc.escudero.all.results import get_broadband_monthly_meas
from antarc.escudero.all.results.plot_broadband_esc_pwrf_era5_case_rsrc import (
    get_broadband,
)

get_broadband_monthly_meas = (
    get_broadband_monthly_meas.get_broadband_monthly_meas
)
get_era5_broadband_monthly_warm = (
    get_era5_broadband_monthly_warm.get_era5_broadband_monthly_warm
)

from antarc.era5_read_station_broadband import read_era5_broadband_down_by_file

make_qc_figs = False


# def getfilenames_daterange(
#     direc: str, flen: int, pfix: str, dt1: dt.datetime, dt2: dt.datetime
# ) -> list:
#     """
#     Get a list of filenames for dates between start and end date
#     @param direc  Directory where files are
#     @param flen  Length of file
#     @param pfix  File prefix
#     @param dt1  Start date
#     @param dt2  End date
#     """
#     plen = len(pfix)
#     date1 = int(dt1.strftime("%Y%m%d"))
#     date2 = int(dt2.strftime("%Y%m%d"))

#     fnames = []
#     for year in range(date1.year, date2.year + 1):
#         fnames0 = getfilenames(direc + "year/", flen, pfix)
#     fnames.append(fnames0)
#     dates = np.array([float(fname[plen : plen + 8]) for fname in fnames])
#     ikeep = np.intersect1d(
#         np.where(dates >= date1)[0], np.where(dates <= date2)[0]
#     )
#     return [direc + fnames[x] for x in ikeep]


def get_meas(dir_swd, dir_lwd, date1, date2, event_start, event_hours):
    # # # # # # #   Get measurements   # # # # # # # # # # #
    flen = len("esc_swd20220206_0301.csv")
    sfiles = getfilenames_daterange(dir_swd, flen, "esc_swd", date1, date2)
    lfiles = getfilenames_daterange(dir_lwd, flen, "esc_lwd", date1, date2)
    sw_date, sw_flux, dateave, sw_fluxave = get_broadband(sfiles)
    sw_flux[sw_flux < 0] = 0
    sw_fluxave = np.array(sw_fluxave)
    sw_fluxave[sw_fluxave < 0] = 0
    lw_date, lw_flux, lw_dateave, lw_fluxave = get_broadband(lfiles)
    if (
        sum([(x - y).total_seconds() for x, y, in zip(lw_dateave, dateave)])
        > 0
    ):
        raise ValueError("dates differ")

    # Get the measurements for the event
    dta = [event_start + dt.timedelta(hours=x) for x in range(event_hours)]
    dta_sec = [(x - event_start).total_seconds() for x in dta]
    # Meas SWD
    dateave_sec = [(x - event_start).total_seconds() for x in dateave]
    sw_event = np.interp(dta_sec, dateave_sec, sw_fluxave)
    # Meas LWD
    dateave_sec = [(x - event_start).total_seconds() for x in lw_dateave]
    lw_event = np.interp(dta_sec, dateave_sec, lw_fluxave)
    meas_event = {"swd": sw_event, "lwd": lw_event, "tot": sw_event + lw_event}

    return meas_event, dta


def get_event_flux(
    date1,
    date2,
    dir_swd,
    dir_lwd,
    event_start,
    event_hours,
    era5,
    pwrf,
):
    edates = era5["date"]

    # # # # # # #   Get measurements   # # # # # # # # # # #
    flen = len("esc_swd20220206_0301.csv")
    sfiles = getfilenames_daterange(dir_swd, flen, "esc_swd", date1, date2)
    lfiles = getfilenames_daterange(dir_lwd, flen, "esc_lwd", date1, date2)
    sw_date, sw_flux, dateave, sw_fluxave = get_broadband(sfiles)
    sw_flux[sw_flux < 0] = 0
    sw_fluxave = np.array(sw_fluxave)
    sw_fluxave[sw_fluxave < 0] = 0
    lw_date, lw_flux, lw_dateave, lw_fluxave = get_broadband(lfiles)
    if (
        sum([(x - y).total_seconds() for x, y, in zip(lw_dateave, dateave)])
        > 0
    ):
        raise ValueError("dates differ")

    # Get the measurements for the event
    dta = [event_start + dt.timedelta(hours=x) for x in range(event_hours)]
    dta_sec = [(x - event_start).total_seconds() for x in dta]
    # Meas SWD
    dateave_sec = [(x - event_start).total_seconds() for x in dateave]
    sw_event = np.interp(dta_sec, dateave_sec, sw_fluxave)
    # Meas LWD
    dateave_sec = [(x - event_start).total_seconds() for x in lw_dateave]
    lw_event = np.interp(dta_sec, dateave_sec, lw_fluxave)
    meas_event = {"swd": sw_event, "lwd": lw_event, "tot": sw_event + lw_event}

    # Get ERA5 for the month for all years and for the event
    # era_2c, era = get_era5_bbnd_mon_warm(date1_long, date2_long, case_mon, 2)
    # era_dt_sec = [(x - event_start).total_seconds() for x in era[imon]["date"]]
    era_dt_sec = [(x - event_start).total_seconds() for x in edates]
    era_event = {
        "swd": np.interp(dta_sec, era_dt_sec, era5["swd"]),
        "lwd": np.interp(dta_sec, era_dt_sec, era5["lwd"]),
        "swd_clr": np.interp(dta_sec, era_dt_sec, era5["swd_clr"]),
        "lwd_clr": np.interp(dta_sec, era_dt_sec, era5["lwd_clr"]),
    }
    era_event["tot"] = era_event["swd"] + era_event["lwd"]
    era_event["tot_clr"] = era_event["swd_clr"] + era_event["lwd_clr"]

    # For PWRF, used dates on the hour
    # Since the Polar WRF GLW values are instantaneous, assume they vary
    # linearly over the hour and interpolate
    estart = event_start.replace(minute=0, second=0)
    dtb = [estart + dt.timedelta(hours=x) for x in range(event_hours)]
    dtb_sec = [(x - event_start).total_seconds() for x in dtb]
    pwrf_sec = [(x - event_start).total_seconds() for x in pwrf["date"]]
    pwrf_event = {
        "lwd": np.interp(dtb_sec, pwrf_sec, np.array(pwrf["lwd"])),
        "swd": np.interp(dtb_sec, pwrf_sec, np.array(pwrf["swd"])),
    }
    pwrf_event["tot"] = pwrf_event["lwd"] + pwrf_event["swd"]

    if make_qc_figs:
        # QC
        plt.figure(num=5, clear=True)
        plt.subplot(2, 1, 1)
        plt.plot(dateave, sw_fluxave, label="meas")
        plt.plot(pwrf["date"], pwrf["swd"], label="PWRF")
        plt.plot(dta, sw_event, ".", label="meas; event")
        plt.plot(dta, era_event["swd"], ".", label="ERA5; event")
        plt.plot(dta, pwrf_event["swd"], ".", label="PWRF; event")
        plt.subplot(2, 1, 2)
        plt.plot(dateave, lw_fluxave, label="meas")
        plt.plot(pwrf["date"], pwrf["lwd"], label="PWRF")
        plt.plot(dta, lw_event, ".", label="meas; event")
        plt.plot(dta, era_event["lwd"], ".", label="ERA5; event")
        plt.plot(dta, pwrf_event["lwd"], ".", label="PWRF; event")
        plt.legend()

    return meas_event, era_event, pwrf_event


def get_event_forcing(
    date1,
    date2,
    dir_swd,
    dir_lwd,
    event_start,
    event_hours,
    era5,
    pwrf,
):

    meas_event, era_event, pwrf_event = get_event_flux(
        date1,
        date2,
        dir_swd,
        dir_lwd,
        event_start,
        event_hours,
        era5,
        pwrf,
    )

    for key in ["swd", "lwd", "tot"]:
        meas_event[key] -= era_event[key + "_clr"]
        pwrf_event[key] -= era_event[key + "_clr"]
        era_event[key] -= era_event[key + "_clr"]

    return meas_event, era_event, pwrf_event


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


# def get_broadband(fnames):
#     """
#     Get hourly averaged broadband fluxes
#     """
#     rad = []
#     mydates = []
#     dateave = []
#     radave = []

#     for filename in fnames:

#         if (
#             filename[filename.rfind("/") :] == "/esc_lwd20220209_2035.csv"
#             or filename[filename.rfind("/") :] == "/esc_swd20220209_2033.csv"
#         ):
#             # Add nans in the missing region
#             nanvec = np.array([np.nan, np.nan])
#             dtimevec = np.array(
#                 [dt.datetime(2022, 2, 7, 21, 24), dt.datetime(2022, 2, 9, 20)]
#             )
#             rad = np.hstack([rad, nanvec])
#             mydates = np.hstack([mydates, dtimevec])

#         # Read in the csv file for ifile
#         df = pd.read_csv(filename, delimiter=",")

#         # Get dates, skipping the 0th element, which is the lower part of the header
#         dates_str = df["Date"] + " " + df["Time"]

#         # Get the date
#         dt0 = [dt.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in dates_str]

#         # Attach rad and date for this file to all rads/dates (rad as float)
#         mydates = np.hstack([mydates, dt0])
#         rad = np.hstack([rad, np.array([float(x) for x in df["Radiation"]])])

# # Get the hourly averages
# radsum = 0
# naves = 0
# hour = mydates[0].hour
# dateave = []
# radave = []
# for i, mydate in enumerate(mydates):
#     if mydate.hour == hour:
#         radsum += rad[i]
#         naves += 1
#     else:
#         dt0 = mydates[i - 1]
#         thisdate = dt.datetime(
#             dt0.year, dt0.month, dt0.day, dt0.hour, 30, 30
#         )
#         dateave.append(thisdate)
#         radave.append(radsum / naves)
#         hour = mydate.hour
#         radsum = rad[i]
#         naves = 1

#     return mydates, rad, dateave, radave


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
