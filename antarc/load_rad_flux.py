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
from antarc.getfilenames import getfiles_for_daterange


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
    files, _ = getfiles_for_daterange(direc, fmt, start_dt, end_dt)

    # Load the pyranometer data
    dtimes = np.array([])
    flux = np.array([])
    for i, file in enumerate(files):
        pir = pd.read_csv(direc + file)

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
