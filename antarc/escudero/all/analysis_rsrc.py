#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 15:33:27 2023

@author: prowe
"""

# Dependencies
import datetime as dt
import numpy as np
import pytz
from matplotlib import pyplot as plt
import os


# My modules
from antarc.escudero.all.read_era5_bbnd import (
    read_era5_broadband_down,
    read_era5_broadband_down_utc,
)
from antarc.load_rad_flux import load_rad_flux_utc


def get_files(direc, fmt, date1, date2):
    samplefile = date1.strftime(fmt)
    allfiles = os.listdir(direc)
    allfiles.sort()
    fnames = [
        x
        for x in allfiles
        if len(x) == len(samplefile)
        and (
            dt.datetime.strptime(x, fmt) >= date1
            and dt.datetime.strptime(x, fmt) < date2
        )
    ]
    return fnames


def get_swd_qc(base_dir: str, fmt: str, date1, date2):
    """
    Load, QC, and return pyranometer data
    """
    direc = f"{base_dir}{date1.strftime('%Y')}/"
    swd_date0, swd = load_rad_flux_utc(direc, fmt, date1, date2)
    swd = np.array(swd)
    swd_date = np.array(swd_date0)
    # Quality control
    swd[swd < 0] = 0

    return swd_date, swd


def get_lwd_qc(base_dir, fmt, date1, date2):
    """
    Get the LWD flux where the innermost directory is the year,
    and quality control the results
    """
    direc = f"{base_dir}{date1.strftime('%Y')}/"

    lwd_date, lwd = load_rad_flux_utc(direc, fmt, date1, date2)

    # Quality control
    lwd = np.array(lwd)
    lwd[lwd < 130] = np.nan
    lwd[lwd > 500] = np.nan

    return lwd_date, lwd


def get_hourly_bbnd_meas(lwd_date, lwd, swd_date, swd):

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


def get_era5_escudero_byyear(
    direc: str, fmt: str, date1, date2, lat: float, lon: float
):
    """
    Get ERA5 results at Escudero
    @param lat  The Escudero latitude
    @param lon  The Escudero longitude
    @returns  Dictionary with ERA5 fields
    """
    era5bb = read_era5_broadband_down(direc, fmt, date1, date2)

    if not era5bb:
        return None

    # Get closest grid point
    i = np.argmin(abs(era5bb["lat"].data - lat))
    j = np.argmin(abs(era5bb["lon"].data - lon))

    era5 = {
        "date": era5bb["date"],
        "swd": era5bb["swd"][:, i, j].data,
        "lwd": era5bb["lwd"][:, i, j].data,
        "swd_clr": era5bb["swd_clr"][:, i, j].data,
        "lwd_clr": era5bb["lwd_clr"][:, i, j].data,
    }

    return era5


def get_era5_escudero_byyear_utc(
    direc: str, fmt: str, date1, date2, lat: float, lon: float
):
    """
    Get ERA5 results at Escudero
    @param lat  The Escudero latitude
    @param lon  The Escudero longitude
    @returns  Dictionary with ERA5 fields
    """
    era5bb = read_era5_broadband_down_utc(direc, fmt, date1, date2)

    if not era5bb:
        return None

    # Get closest grid point
    i = np.argmin(abs(era5bb["lat"].data - lat))
    j = np.argmin(abs(era5bb["lon"].data - lon))

    era5 = {
        "date": era5bb["date"],
        "swd": era5bb["swd"][:, i, j].data,
        "lwd": era5bb["lwd"][:, i, j].data,
        "swd_clr": era5bb["swd_clr"][:, i, j].data,
        "lwd_clr": era5bb["lwd_clr"][:, i, j].data,
    }

    return era5


def get_era5_esc_utc(
    direc: str, fmt: str, date1, date2, lat: float, lon: float
):
    """
    Get ERA5 results at Escudero
    @param lat  The Escudero latitude
    @param lon  The Escudero longitude
    @returns  Dictionary with ERA5 fields
    """
    dir_yr = direc + f"{date1.strftime('%Y')}/"
    era5 = get_era5_escudero_byyear_utc(dir_yr, fmt, date1, date2, lat, lon)

    if not era5:
        return None

    if date2.year == date1.year:
        return era5

    elif date2.year != date1.year + 1:
        raise ValueError("date2 year cannot be more than 1 year later")

    # print(f"get next year, starting on day after last: {era5['date'][-1]}")
    date1 = era5["date"][-1] + dt.timedelta(days=1)
    # Crop off hours, to start at beginning of day, because files are per day
    date1 = dt.datetime(date1.year, date1.month, date1.day)
    dir_yr = direc + f"{date2.strftime('%Y')}/"
    era5b = get_era5_escudero_byyear_utc(dir_yr, fmt, date1, date2, lat, lon)

    if not era5b:
        return era5

    for key in era5.keys():
        era5[key] = np.hstack([era5[key], era5b[key]])

    return era5


def get_era5_escudero(
    direc: str, fmt: str, date1, date2, lat: float, lon: float
):
    """
    Get ERA5 results at Escudero
    @param lat  The Escudero latitude
    @param lon  The Escudero longitude
    @returns  Dictionary with ERA5 fields
    """
    direc_year = direc + f"{date1.strftime('%Y')}/"
    era5 = get_era5_escudero_byyear(direc_year, fmt, date1, date2, lat, lon)

    if not era5:
        return None

    if date2.year == date1.year:
        # era5["date"] = np.array(
        #     [x.replace(tzinfo=pytz.utc) for x in era5["date"]]
        # )
        era5["date"] = np.array(era5["date"])
        return era5

    elif date2.year != date1.year + 1:
        raise ValueError("date2 year cannot be more than 1 year later")

    # print(f"get next year, starting on day after last: {era5['date'][-1]}")
    date1 = era5["date"][-1] + dt.timedelta(days=1)
    # Crop off hours, to start at beginning of day, because files are per day
    date1 = dt.datetime(date1.year, date1.month, date1.day, tzinfo=pytz.utc)
    direc_year = direc + f"{date2.strftime('%Y')}/"
    era5b = get_era5_escudero_byyear(direc_year, fmt, date1, date2, lat, lon)

    if not era5b:
        return era5

    for key in era5.keys():
        era5[key] = np.hstack([era5[key], era5b[key]])

    era5["date"] = np.array([x.replace(tzinfo=pytz.utc) for x in era5["date"]])

    return era5


def get_end_of_hour(date):
    """
    Return the date falling at the end of the current hour; that is,
    the beginning of the next hour
    """
    new = dt.datetime(date.year, date.month, date.day, date.hour)
    # new = new.replace(tzinfo=pytz.utc)
    return new + dt.timedelta(hours=1)


def ave_over_hour(dates, vals):
    """
    Average over an hour, and assign the date at the end of the hour
    """
    if len(dates) == 0:
        return [], [], []

    current_hour = dates[0].strftime("%Y%m%d%H")
    current_date = dates[0]
    naves = 0
    tot = 0
    val_ave = []
    date_ave = []  # Do not include the first hour
    nave = []
    for i, date in enumerate(dates):
        if date.strftime("%Y%m%d%H") == current_hour:
            # Add to average
            naves += 1
            tot += vals[i]

        else:
            # Finish up last
            nave.append(naves)
            val_ave.append(tot / naves)
            # assign the end of the hour
            date_ave.append(get_end_of_hour(current_date))
            if date_ave[-1].minute > 0 or date_ave[-1].second > 0:
                print()
            # Prep next
            current_date = date
            current_hour = date.strftime("%Y%m%d%H")
            naves = 1
            tot = vals[i]

    # Finish off the last
    nave.append(naves)
    nave.append(naves)
    val_ave.append(tot / naves)
    date_ave.append(get_end_of_hour(current_date))

    return nave, date_ave, val_ave


def ave_over_every_hour_fxd(mdates, vals, date1, date2):
    """
    Average over an hour from the start date to the end date, without replacing
    """
    # Set up a vector of hours from date1 to date2
    first_hour = date1.replace(minute=0, second=0) + dt.timedelta(hours=1)
    last_hour = date2.replace(minute=0, second=0) + dt.timedelta(hours=1)
    hour = first_hour
    hours = []
    while hour <= last_hour:
        hours.append(hour)
        hour += dt.timedelta(hours=1)

    # Quality control times
    tdiff = [
        (hours[i + 1] - y).total_seconds() for i, y in enumerate(hours[:-1])
    ]
    if not np.allclose(tdiff, 3600):
        raise ValueError("one or more times is bad")

    nave = np.zeros(len(hours))
    tot = np.zeros(len(hours))
    val_ave = np.nan + np.zeros(len(hours))

    if len(mdates) == 0:
        return nave, hours, val_ave

    # Skip to where the measured data starts
    ihour = 0
    while mdates[0] > hours[ihour] and ihour < len(hours):
        ihour += 1

    i = 0
    date = mdates[0]
    for ihour in range(ihour, len(hours)):
        hour = hours[ihour]
        while date <= hour and i < len(mdates):
            # Still in the current hour: add to average
            nave[ihour] += 1
            try:
                tot[ihour] += vals[i]
            except:
                print("pause")
            # Prep for next round
            i += 1
            if i >= len(mdates):
                break
            date = mdates[i]

    val_ave[nave > 0] = tot[nave > 0] / nave[nave > 0]

    # current_date = dates[0]
    # current_hour_str = dates[0].strftime("%Y%m%d%H")
    # naves = 0
    # tot = 0
    # val_ave = []
    # date_ave = []  # Do not include the first hour
    # nave = []
    # for i, date in enumerate(dates):
    #     if date.strftime("%Y%m%d%H") == current_hour_str:
    #         # Still in the current hour: add to average
    #         naves += 1
    #         tot += vals[i]

    #     else:
    #         # Started next hour. Finish up last
    #         nave.append(naves)
    #         val_ave.append(tot / naves)
    #         # assign the end of the hour
    #         date_ave.append(get_end_of_hour(current_date))
    #         if date_ave[-1].minute > 0 or date_ave[-1].second > 0:
    #             raise ValueError("what is this?")

    #         # Check if the next hour is skipped
    #         if current_hour + dt.timedelta(hours=1) > date:
    #             raise ValueError("An hour was skipped!")

    #         # Prep next
    #         current_date = date
    #         current_hour = date.replace(minute=0, second=0)
    #         current_hour_str = date.strftime("%Y%m%d%H")
    #         naves = 1
    #         tot = vals[i]

    # Finish off the last
    # nave.append(naves)
    # nave.append(naves)
    # val_ave.append(tot / naves)
    # date_ave.append(get_end_of_hour(current_date))

    # # Quality control times
    # tdiff = [
    #     (date_ave[i + 1] - y).total_seconds()
    #     for i, y in enumerate(date_ave[:-1])
    # ]

    # if np.any(np.array(tdiff) != 3600):
    #     raise ValueError("Times have a problem")

    return nave, hours, val_ave


def ave_over_every_hour(mdates, vals, date1, date2):
    """
    Average over an hour, and assign the date at the end of the hour
    """
    # Set up a vector of hours from date1 to date2
    first_hour = date1.replace(minute=0, second=0) + dt.timedelta(hours=1)
    last_hour = date2.replace(minute=0, second=0) + dt.timedelta(hours=1)
    hour = first_hour
    hours = []
    while hour <= last_hour:
        hours.append(hour)
        hour += dt.timedelta(hours=1)

    # # Quality control times: make sure they are about 1 hour apart (3600 s)
    # tdiff = [
    #     (hours[i + 1] - y).total_seconds() for i, y in enumerate(hours[:-1])
    # ]
    # if not np.allclose(tdiff, 3600):
    #     raise ValueError("one or more times is bad")

    nave = np.zeros(len(hours))
    tot = np.zeros(len(hours))
    val_ave = np.nan + np.zeros(len(hours))

    if len(mdates) == 0:
        return nave, hours, val_ave

    # Skip to where the measured data starts
    ihour = 0
    while mdates[0] > hours[ihour] and ihour < len(hours):
        ihour += 1

    # Average all points with an hour, leaving values nan if the entire
    # hour is missing. This will assume missing data has the same average
    # as the average number of points for the data that exists
    i = 0
    date = mdates[0]
    for ihour in range(ihour, len(hours)):
        hour = hours[ihour]
        while date <= hour and i < len(mdates):
            # Still in the current hour: add to average
            nave[ihour] += 1
            tot[ihour] += vals[i]

            # TODO: this is new!
            # If we are on the hour, then we are about to advance the hour
            # Include this point in the average for the upcoming hour
            if date == hour and len(tot) <= ihour:
                tot[ihour + 1] += vals[i]
                nave[ihour + 1] += 1

            # Advance index to date and date for next time through
            i += 1
            if i >= len(mdates):
                break
            date = mdates[i]

    val_ave[nave > 0] = tot[nave > 0] / nave[nave > 0]

    # current_date = dates[0]
    # current_hour_str = dates[0].strftime("%Y%m%d%H")
    # naves = 0
    # tot = 0
    # val_ave = []
    # date_ave = []  # Do not include the first hour
    # nave = []
    # for i, date in enumerate(dates):
    #     if date.strftime("%Y%m%d%H") == current_hour_str:
    #         # Still in the current hour: add to average
    #         naves += 1
    #         tot += vals[i]

    #     else:
    #         # Started next hour. Finish up last
    #         nave.append(naves)
    #         val_ave.append(tot / naves)
    #         # assign the end of the hour
    #         date_ave.append(get_end_of_hour(current_date))
    #         if date_ave[-1].minute > 0 or date_ave[-1].second > 0:
    #             raise ValueError("what is this?")

    #         # Check if the next hour is skipped
    #         if current_hour + dt.timedelta(hours=1) > date:
    #             raise ValueError("An hour was skipped!")

    #         # Prep next
    #         current_date = date
    #         current_hour = date.replace(minute=0, second=0)
    #         current_hour_str = date.strftime("%Y%m%d%H")
    #         naves = 1
    #         tot = vals[i]

    # Finish off the last
    # nave.append(naves)
    # nave.append(naves)
    # val_ave.append(tot / naves)
    # date_ave.append(get_end_of_hour(current_date))

    # # Quality control times
    # tdiff = [
    #     (date_ave[i + 1] - y).total_seconds()
    #     for i, y in enumerate(date_ave[:-1])
    # ]

    # if np.any(np.array(tdiff) != 3600):
    #     raise ValueError("Times have a problem")

    return nave, hours, val_ave


def ave_over_time_offset(dates, vals, toff):
    # apply the time offset
    dates = np.array([x + toff for x in dates])
    return ave_over_hour(dates, vals)


def get_stats(era5_date, era5_val, date, val):
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

    # Make sure there are no nans in era5
    npts_era5 = len(era5_val[~np.isnan(era5_val)])
    if npts_era5 != len(era5_val):
        raise ValueError("One or more ERA5 values is nan")
    if npts_era5 == 0:
        raise ValueError(f"No ERA5 data for {era5_date[0] - era5_date[-1]}")

    # We should have a 1:1 correspondence with dates
    if len(era5_date) == len(date):
        tdiff = np.array(
            [(x - y).total_seconds() for x, y in zip(era5_date, date)]
        )
        if np.any(np.abs(tdiff) > 1):
            raise ValueError("ERA5 and meas dates are not same")
    else:
        # Remove era5 dates that are before the measurement dates
        i1 = 0
        while era5_date[i1] < date[0]:
            i1 += 1
            # Quit if we run out of dates
            if i1 > len(era5_date) - 1:
                return None
        era5_val = era5_val[i1:]
        era5_date = era5_date[i1:]
        # Keep all dates that are the same
        i2 = 0
        im2 = 0
        while era5_date[i2] == date[im2]:
            i2 += 1
            im2 += 1
            if i2 > len(era5_date) - 1 or im2 > len(date) - 1:
                break
        era5_val = era5_val[: i2 + 1]
        era5_date = era5_date[: i2 + 1]
        val = val[: im2 + 1]
        date = date[: im2 + 1]

    npts_meas = len(val[~np.isnan(val)])

    if npts_meas == 0:
        return {
            "era5": np.nan,
            "meas": np.nan,
            "bias": np.nan,
            "rms": np.nan,
            "npts": np.nan,
        }

    return {
        "era5": np.nanmean(era5_val[~np.isnan(val)]),
        "meas": np.nanmean(val),
        "bias": np.nanmean(era5_val - val),
        "rms": np.sqrt(np.nanmean((era5_val - val) ** 2)),
        "npts": npts_meas,
    }

    # squarediff = 0
    # absdiff = 0
    # squarediff_fract = 0
    # absdiff_fract = 0
    # npts = 0
    # iera = 0
    # imeas = 0
    # while iera < len(era5_date) and imeas < len(date):
    #     if era5_date[iera] == date[imeas]:
    #         if not np.isnan(val[imeas]):
    #             npts += 1
    #             absdiff0 = era5_val[iera] - val[imeas]
    #             if val[imeas] == 0 and abs(era5_val[iera]) <= 0.1:
    #                 absdiff0_fract = 0
    #             elif val[imeas] == 0 and abs(era5_val[iera]) > 0.1:
    #                 absdiff0_fract = absdiff0 / era5_val[iera]
    #                 print("warning: unexpectedly large era5 value")
    #             else:
    #                 absdiff0_fract = absdiff0 / val[imeas]
    #             absdiff += absdiff0
    #             squarediff += absdiff0**2
    #             absdiff_fract += absdiff0_fract
    #             squarediff_fract += absdiff0_fract**2
    #         iera += 1
    #         imeas += 1
    #     elif era5_date[iera] < date[imeas]:
    #         iera += 1
    #     elif date[imeas] < era5_date[iera]:
    #         imeas += 1
    #     else:
    #         raise ValueError("something is wrong")
    # if npts == 0:
    #     return np.nan, np.nan, npts, np.nan, np.nan

    # rmsdiff = np.sqrt(squarediff / npts)
    # meandiff = absdiff / npts
    # rmsdiff_pc = 100 * np.sqrt(squarediff_fract / npts)
    # absdiff_pc = 100 * absdiff_fract / npts

    # return meandiff, rmsdiff, npts, rmsdiff_pc, absdiff_pc


def get_forcing_stats(era5_date, era5_val, era5_clr, date, val):
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

    # Make sure there are no nans in era5
    npts_era5 = len(era5_val[~np.isnan(era5_val)])
    if npts_era5 != len(era5_val):
        raise ValueError("One or more ERA5 values is nan")
    if npts_era5 == 0:
        raise ValueError(f"No ERA5 data for {era5_date[0] - era5_date[-1]}")

    # We should have a 1:1 correspondence with dates
    if len(era5_date) == len(date):
        tdiff = np.array(
            [(x - y).total_seconds() for x, y in zip(era5_date, date)]
        )
        if np.any(np.abs(tdiff) > 1):
            raise ValueError("ERA5 and meas dates are not same")
    else:
        # Remove era5 dates that are before the measurement dates
        i1 = 0
        while era5_date[i1] < date[0]:
            i1 += 1
            # Quit if we run out of dates
            if i1 > len(era5_date) - 1:
                return None
        era5_val = era5_val[i1:]
        era5_date = era5_date[i1:]
        # Keep all dates that are the same
        i2 = 0
        im2 = 0
        while era5_date[i2] == date[im2]:
            i2 += 1
            im2 += 1
            if i2 > len(era5_date) - 1 or im2 > len(date) - 1:
                break
        era5_val = era5_val[: i2 + 1]
        era5_date = era5_date[: i2 + 1]
        val = val[: im2 + 1]
        date = date[: im2 + 1]

    npts_meas = len(val[~np.isnan(val)])

    if npts_meas == 0:
        return {
            "era5": np.nan,
            "meas": np.nan,
            "bias": np.nan,
            "rms": np.nan,
            "npts": np.nan,
        }

    return {
        "era5": np.nanmean((era5_val - era5_clr)[~np.isnan(val)]),
        "meas": np.nanmean(val - era5_clr),
        "bias": np.nanmean(era5_val - val),
        "rms": np.sqrt(np.nanmean((era5_val - val) ** 2)),
        "npts": npts_meas,
    }


def get_rmsdiff(era5_date, era5_val, date, val):
    """
    Get the rms difference between ERA5 and measured value
    """

    _, rmsdiff = get_stats(era5_date, era5_val, date, val)

    return rmsdiff


def get_rmsdiff_toff(era5_date, era5_val, date, val, toff):
    """
    Get the rms difference for various time offsets
    """
    nave, date, val = ave_over_time_offset(date, val, toff)
    return get_rmsdiff(era5_date, era5_val, date, val)


def qc_era5_date(era5, swd_date, swd):
    """Prove to myself that ERA5 accumulates to the end of the hour"""
    min_offsets = [-60, -45, -30, -15, 0, 15, 30, 45, 60]
    toffs = [dt.timedelta(minutes=x) for x in min_offsets]
    rms = [
        get_rmsdiff_toff(era5["date"], era5["swd"], swd_date, swd, toff)
        for toff in toffs
    ]

    plt.figure()
    plt.plot(min_offsets, rms)
