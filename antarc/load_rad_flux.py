#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 11:15:15 2022

@author: prowe
"""

# # # # # # # #
import numpy as np
import datetime as dt
import pytz
import pandas as pd

# My modules
from antarc.getfilenames import getfiles_for_daterange, getfnames_daterange_utc


def load_rad_flux(
    direc: str, fmt: str, start_dt: dt.datetime, end_dt: dt.datetime
) -> tuple[np.ndarray, np.ndarray]:
    """
    Load all downward radiative flux data between two dates
    @param direc  Directory
    @param fmt  File format
    @param start_dt  Starting date
    @param end_dt  Ending date
    """

    if direc[-5:] != str(start_dt.year) + "/":
        raise ValueError("Directory format is incorrect")

    files = []
    files_wo_direc = []
    for year in range(start_dt.year, end_dt.year + 1):
        direc = direc[:-5] + str(year) + "/"
        files0, _ = getfiles_for_daterange(direc, fmt, start_dt, end_dt)
        files += [direc + x for x in files0]
        files_wo_direc += files0

    # QC that the same file does not exist in different year directories
    if len(files_wo_direc) != len(set(files_wo_direc)):
        raise ValueError("Same file found in two different directories")

    # Load the pyranometer data
    dtimes = np.array([])
    flux = np.array([])
    for i, file in enumerate(files):
        pir = pd.read_csv(file)

        # Get dates from pyranometer file data
        ntimes = len(pir)
        fmt = "%Y-%m-%d%H:%M:%S"
        dtimes0 = [
            dt.datetime.strptime(pir["Date"][i] + pir["Time"][i], fmt).replace(
                tzinfo=pytz.utc
            )
            for i in range(ntimes)
        ]
        dtimes = np.hstack([dtimes, dtimes0])
        flux = np.hstack([flux, pir["Radiation"].values])

    inds = np.where(np.logical_and(dtimes >= start_dt, dtimes <= end_dt))[0]

    return dtimes[inds], flux[inds]


def load_rad_flux_utc(
    direc: str, fmt: str, start_dt: dt.datetime, end_dt: dt.datetime
) -> tuple[np.ndarray, np.ndarray]:
    """
    Load all downward radiative flux data between two dates
    assuming all dates are already in UTC, and thus we can use
    timezone unaware dates
    @param direc  Directory
    @param fmt  File format
    @param start_dt  Starting date (in UTC; timezone unaware)
    @param end_dt  Ending date (in UTC; timezone unaware)
    """

    if direc[-5:] != str(start_dt.year) + "/":
        raise ValueError("Directory format is incorrect")

    files = []
    files_wo_direc = []
    for year in range(start_dt.year, end_dt.year + 1):
        direc = direc[:-5] + str(year) + "/"
        files0, _ = getfnames_daterange_utc(direc, fmt, start_dt, end_dt)
        files += [direc + x for x in files0]
        files_wo_direc += files0

    # QC that the same file does not exist in different year directories
    if len(files_wo_direc) != len(set(files_wo_direc)):
        raise ValueError("Same file found in two different directories")

    # Load the pyranometer data
    dtimes = np.array([])
    flux = np.array([])
    for i, file in enumerate(files):
        pir = pd.read_csv(file)

        # Get dates from pyranometer file data
        ntimes = len(pir)
        fmt = "%Y-%m-%d%H:%M:%S"
        dtimes0 = [
            dt.datetime.strptime(pir["Date"][i] + pir["Time"][i], fmt)
            for i in range(ntimes)
        ]
        dtimes = np.hstack([dtimes, dtimes0])
        flux = np.hstack([flux, pir["Radiation"].values])

    inds = np.where(np.logical_and(dtimes >= start_dt, dtimes <= end_dt))[0]

    return dtimes[inds], flux[inds]
