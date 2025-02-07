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
    get_era5_esc_utc,
    get_stats,
    get_lwd_qc,
    get_swd_qc,
    ave_over_every_hour,
)


def get_broadband_stats(era5, meas, npts):
    """
    Get broadband stats
    """
    if not era5 or (len(meas["lwd_date"]) == 0 and len(meas["swd_date"]) == 0):
        return None, None
    lwd = get_broadband_stats_key(era5, meas, "lwd", npts)
    swd = get_broadband_stats_key(era5, meas, "swd", npts)
    return lwd, swd


def get_broadband_stats_key(era5, meas, key, npts):
    """
    Get broadband stats for lwd or swd
    """
    datestr = key + "_date"
    return get_stats(era5["date"], era5[key], meas[datestr], meas[key])


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


def get_broadband_monthly_era5(
    years: list[int],
    era_dir: str,
    era5fmt: str,
    lat: float,
    lon: float,
):

    mon_keys = ["era5"]
    mon_keys_all = mon_keys + ["era5_all"]
    months = range(1, 13)
    ones = np.ones([len(years), len(months)])
    monthly = {
        "lwd": {key: np.nan * ones for key in mon_keys_all},
        "swd": {key: np.nan * ones for key in mon_keys_all},
        "lwd_clr": {key: np.nan * ones for key in mon_keys_all},
        "swd_clr": {key: np.nan * ones for key in mon_keys_all},
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

            # Get ERA5 data
            era5 = get_era5_esc_utc(era_dir, era5fmt, date1, date2, lat, lon)

            # Get ERA5 monthly all (not just matched)
            if era5:
                monthly["lwd"]["era5_all"][iyr, imo] = np.nanmean(era5["lwd"])
                monthly["swd"]["era5_all"][iyr, imo] = np.nanmean(era5["swd"])
                monthly["lwd_clr"]["era5_all"][iyr, imo] = np.nanmean(
                    era5["lwd_clr"]
                )
                monthly["swd_clr"]["era5_all"][iyr, imo] = np.nanmean(
                    era5["swd_clr"]
                )

    return monthly
