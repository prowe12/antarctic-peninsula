#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 10:47:58 2024

@author: prowe
"""

import numpy as np
import datetime as dt
import calendar

from antarc.escudero.all.analysis_rsrc import (
    get_lwd_qc,
    get_swd_qc,
    ave_over_every_hour,
)


def get_stats(date, val):
    """
    Get the mean ERA5 and mean measurement at measurement times and get
    the bias and rms differences between the ERA5 and measured value:
    - mean ERA5
    - mean measured
    - mean difference
    - rms difference
    - mean difference (%) (removed)
    - rms difference (%) (removed)
    - npts
    """

    npts_meas = len(val[~np.isnan(val)])

    if npts_meas == 0:
        return {
            "mean": np.nan,
            "max": np.nan,
            "min": np.nan,
            "rms": np.nan,
            "npts": np.nan,
        }

    return {
        "mean": np.nanmean(val),
        "max": np.nanmax(val),
        "min": np.nanmin(val),
        "rms": np.sqrt(np.mean((val[np.isnan(val) == False]) ** 2)),
        "npts": npts_meas,
    }


# def get_broadband_stats(era5, meas, npts):
#     """
#     Get broadband stats
#     """
#     if not era5 or (len(meas["lwd_date"]) == 0 and len(meas["swd_date"]) == 0):
#         return None, None
#     lwd = get_broadband_stats_key(era5, meas, "lwd", npts)
#     swd = get_broadband_stats_key(era5, meas, "swd", npts)
#     return lwd, swd


# def get_broadband_stats_key(era5, meas, key, npts):
#     """
#     Get broadband stats for lwd or swd
#     """
#     datestr = key + "_date"
#     return get_stats(era5["date"], era5[key], meas[datestr], meas[key])


def get_hourly_meas(lwd_date, lwd, swd_date, swd, date1, date2):
    """
    Average measurements to an hour
    """
    nlwd, lwd_date_a, lwd_a = ave_over_every_hour(lwd_date, lwd, date1, date2)
    nswd, swd_date_a, swd_a = ave_over_every_hour(swd_date, swd, date1, date2)

    if len(lwd_date_a) == 0:
        return {
            "date": swd_date_a,
            "lwd": lwd_a,
            "swd": swd_a,
            "lwd_date": lwd_date_a,
            "swd_date": swd_date_a,
            "lwd_pts": nlwd,
            "swd_pts": nswd,
        }
    elif len(lwd_date_a) != len(swd_date_a):
        raise ValueError("Problem")
    else:
        return {
            "date": swd_date_a,
            "lwd": lwd_a,
            "swd": swd_a,
            "lwd_date": lwd_date_a,
            "swd_date": swd_date_a,
            "lwd_pts": nlwd,
            "swd_pts": nswd,
        }

    return None


def get_broadband_monthly_meas(
    years: list[int],
    lwd_base_dir: str,
    swd_base_dir: str,
    lwd_fmt: str,
    swd_fmt: str,
    lat: float,
    lon: float,
    pts_min: int,
):

    mon_keys = ["mean", "max", "min", "rms", "npts"]
    mon_keys_all = mon_keys
    months = range(1, 13)
    ones = np.ones([len(years), len(months)])
    monthly = {
        "lwd": {key: np.nan * ones for key in mon_keys_all},
        "swd": {key: np.nan * ones for key in mon_keys_all},
    }
    dates = [[[] for i in range(len(months))] for j in range(len(years))]
    count = 0
    for iyr, year in enumerate(years):
        for imo, mon in enumerate(months):
            dtime = dt.datetime(year, mon, 1)
            dates[iyr][imo] = dtime
            count += 1

            date1 = dt.datetime(year, mon, 1)
            print(f"On {date1.year}/{date1.month}")

            # Get the date corresponding to the last hour of the month
            day = calendar.monthrange(year, mon)[1]
            date2 = dt.datetime(year, mon, day, 23, 59)

            # Load in the measured broadband data
            lwd_date, lwd = get_lwd_qc(lwd_base_dir, lwd_fmt, date1, date2)
            swd_date, swd = get_swd_qc(swd_base_dir, swd_fmt, date1, date2)

            # Get hourly averages to compare to ERA5. Drop back an hour
            # to match ERA5, whose first hour comes from the accumulated hour
            # on the previous day
            dt1 = date1 - dt.timedelta(hours=1)
            dt2 = date2 - dt.timedelta(hours=1)
            meas_hr = get_hourly_meas(lwd_date, lwd, swd_date, swd, dt1, dt2)

            # Get the statistics for this date range
            # lwd_stats, swd_stats = get_broadband_stats(meas_hr, pts_min)
            lwd_stats = get_stats(meas_hr["lwd_date"], meas_hr["lwd"])
            swd_stats = get_stats(meas_hr["swd_date"], meas_hr["swd"])

            for key in mon_keys:
                if lwd_stats:
                    monthly["lwd"][key][iyr, imo] = lwd_stats[key]
                if swd_stats:
                    monthly["swd"][key][iyr, imo] = swd_stats[key]

    return monthly
