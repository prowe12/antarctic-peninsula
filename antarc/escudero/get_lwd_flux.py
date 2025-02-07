#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 10:16:21 2022

@author: prowe
"""

# Dependencies
import datetime as dt
import numpy as np

# My modules
from antarc.getfilenames import getfiles_for_daterange
from antarc.load_rad_flux import load_rad_flux


def get_lwd_clear_date(direc, fmt_out, fmt_lo, fmt_hi, date1, date2):
    """
    Get previously run clear sky calculation of longwave downward (LWD)
    radiative fluxes within the wavenumber range of the Escudero pyrgeometer.
    @params direc
    @params fmt_lor  Format for even-lower wavenumber files
    @params fmt_lo  Format for low wavenumber files
    @params fmt_hi  Format for high wavenumber files
    @params date1
    @param date2
    @return  The dates
    @retrun  The total longwave downward radiative fluxes for each date.
    """

    # Get the filenames for all calculated clear-sky LWD fluxes between
    # the date1 and date2 for each set of wavenumbers (low and high, giving
    # lo_files and hi_files).
    # lor_files, lor_dates = getfiles_for_daterange(direc, fmt_lor, date1, date2)
    lo_files, lo_dates = getfiles_for_daterange(direc, fmt_lo, date1, date2)
    hi_files, hi_dates = getfiles_for_daterange(direc, fmt_hi, date1, date2)
    out_files, out_dates = getfiles_for_daterange(direc, fmt_out, date1, date2)

    # QC: warn if no files are found and throw an error if there is not a
    # one-to-one correspondence between files for the low and high sets
    # of wavenumbers, since the fluxes need to be added.
    if len(lo_dates) == 0 and len(hi_dates) == 0:
        return np.array([]), np.array([])
    if len(lo_dates) != len(hi_dates) or np.any(lo_dates != hi_dates):
        for i, lo_date in enumerate(lo_dates):
            if hi_dates[i] != lo_date:
                print(f"Matching file missing for {hi_dates[i] } or {lo_date}")
        raise ValueError("Low and high wavenumber file dates do not match")

    if len(lo_dates) != len(out_dates) or np.any(lo_dates != out_dates):
        for i, lo_date in enumerate(lo_dates):
            if out_dates[i] != lo_date:
                print(f"Matching file missing for {hi_dates[i] } or {lo_date}")
        raise ValueError("Low and out-band wavenumber file dates do not match")

    # Load the calculated clear-sky longwave broadband data for the low
    # and high wavenumber ranges and add them together to get the total
    # flux over the same wavenumber range as the pyrgeometer.
    lwd_clear = []
    for lofile, hifile, outfile in zip(lo_files, hi_files, out_files):
        lwd0 = (
            np.loadtxt(direc + lofile)
            + np.loadtxt(direc + hifile)
            + np.loadtxt(direc + outfile)
        )
        lwd_clear.append(lwd0)

    return lo_dates, np.array(lwd_clear)


def get_lwd_clear(esc_case, lwd_params):
    """
    Get previously run clear sky calculation of longwave downward (LWD)
    radiative fluxes within the wavenumber range of the Escudero pyrgeometer.
    @return  The dates
    @retrun  The total longwave downward radiative fluxes for each date.
    """

    # Get the date range for this particular case study, as the starting date
    # (date1) and the ending date (date2)
    date1 = esc_case.DATE1
    date2 = esc_case.DATE2

    # Get parameters that depend on how the LWD was calculated and stored,
    # including the directory and the format for file names for the set of
    # low wavenumbers (fmt_lo) and high wavenumbers (fmt_hi)
    direc = lwd_params.LWD_CLEAR_DIR
    fmt_out = lwd_params.LWD_CLEAR_FILEFORMAT_OUTNU
    fmt_lo = lwd_params.LWD_CLEAR_FILEFORMAT_LOWNU
    fmt_hi = lwd_params.LWD_CLEAR_FILEFORMAT_HINU

    # Extend date range by 1 day on either side
    # utc = dt.timezone.utc
    date1 = dt.datetime(date1.year, date1.month, date1.day - 4)
    date2 = dt.datetime(date2.year, date2.month, date2.day + 2)
    # , tzinfo=utc)

    # Get the filenames for all calculated clear-sky LWD fluxes between
    # the date1 and date2 for each set of wavenumbers (low and high, giving
    # lo_files and hi_files).
    _, out_dates = getfiles_for_daterange(direc, fmt_out, date1, date2)
    lo_files, lo_dates = getfiles_for_daterange(direc, fmt_lo, date1, date2)
    hi_files, hi_dates = getfiles_for_daterange(direc, fmt_hi, date1, date2)

    # QC: warn if no files are found and throw an error if there is not a
    # one-to-one correspondence between files for the low and high sets
    # of wavenumbers, since the fluxes need to be added.
    if len(lo_dates) == 0 or len(hi_dates) == 0 or len(out_dates) == 0:
        print("Warning: no pyrgeometer files found!")
    if len(lo_dates) != len(hi_dates) or np.any(lo_dates != hi_dates):
        raise ValueError("Low and high wavenumber file dates do not match")
    if len(out_dates) != len(hi_dates) or np.any(out_dates != hi_dates):
        raise ValueError("out and low wavenumber file dates do not match")
    if (
        len(lo_dates) == 0
        or len(hi_dates) == 0
        or len(out_dates) == 0
        or len(out_dates) != len(hi_dates)
        or np.any(out_dates != hi_dates)
    ):
        raise ValueError("probably missing the out-of-band file")

    # Load the calculated clear-sky longwave broadband data for the low
    # and high wavenumber ranges and add them together to get the total
    # flux over the same wavenumber range as the pyrgeometer.
    lwd_clear = []
    for lofile, hifile in zip(lo_files, hi_files):
        lwd0 = np.loadtxt(direc + lofile) + np.loadtxt(direc + hifile)
        outfile = lofile[: lo_files[0].find("238")] + "0_238.txt"
        lwd0 += np.loadtxt(direc + outfile)
        lwd_clear.append(lwd0)

    return lo_dates, lwd_clear


def get_lwd(lwd_params, start_date, end_date) -> tuple[np.ndarray, np.ndarray]:
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
