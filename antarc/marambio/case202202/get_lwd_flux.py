#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 10:16:21 2022

@author: prowe
"""

# Dependencies
import numpy as np
import datetime as dt

# My modules
from antarc.getfilenames import getfiles_for_daterange
from antarc.load_rad_flux import load_rad_flux

# Parameter modules
from antarc.escudero.parameters import esc202202, lwd_params


def get_lwd_clear():
    """
    Get clear sky simulation of longwave radiative flux
    @return  Datetimes
    @retrun  longwave downward radiative flux
    """
    date1 = esc202202.DATE1
    date2 = esc202202.DATE2

    direc = lwd_params.LWD_CLEAR_DIR
    fmt_lo = lwd_params.LWD_CLEAR_FILEFORMAT_LOWNU
    fmt_hi = lwd_params.LWD_CLEAR_FILEFORMAT_HINU

    # Extend date range by 1 day on either side
    utc = dt.timezone.utc
    date1 = dt.datetime(date1.year, date1.month, date1.day - 4, tzinfo=utc)
    date2 = dt.datetime(date2.year, date2.month, date2.day + 2, tzinfo=utc)

    # Load all pyranometer data between two dates
    lo_files, lo_dates = getfiles_for_daterange(direc, fmt_lo, date1, date2)
    hi_files, hi_dates = getfiles_for_daterange(direc, fmt_hi, date1, date2)

    # QC
    if len(lo_dates) != len(hi_dates) or np.any(lo_dates != hi_dates):
        raise ValueError("Low and high wavenumber file dates do not match")

    # Load the calculated clear-sky broadband data
    lwd_clear = []
    for lofile, hifile in zip(lo_files, hi_files):
        lwd0 = np.loadtxt(direc + lofile) + np.loadtxt(direc + hifile)
        lwd_clear.append(lwd0)

    return lo_dates, lwd_clear


def get_lwd(
    start_date=esc202202.DATE1, end_date=esc202202.DATE2
) -> tuple[np.ndarray, np.ndarray]:
    """
    Get measured longwave downward flux
    @return  Datetimes corresponding to flux
    @return  Longwave downward radiative flux
    """

    pyr_dir = lwd_params.LWD_DIR
    pyr_fmt = lwd_params.LWD_FILEFORMAT

    pyr_dir += str(start_date.year) + "/"

    dates, lwd = load_rad_flux(pyr_dir, pyr_fmt, start_date, end_date)
    return dates, lwd
