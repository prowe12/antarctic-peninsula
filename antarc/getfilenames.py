#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 11:34:07 2022

@author: prowe
"""

import numpy as np
import os
import datetime as dt
from typing import Any


def getfilenames(direc: str, flen: int, pfix: str) -> list:
    """
    Get all filenames having the specified length and prefix
    @param direc  Desired directory to get file names from
    @param flen  Length of file names
    @param pfix  File prefix
    @return  A list of file names
    """
    allfnames = np.sort(os.listdir(direc))
    return [x for x in allfnames if len(x) == flen and x[: len(pfix)] == pfix]


def get_filenames(direc: str, fmt: str) -> list:
    """
    Get all filenames having the specified length and prefix
    @param direc  Desired directory to get file names from
    @param fmt  Format of desired files from direc
    @return  A list of file names
    """
    matchlen = fmt.find("%Y")
    match = fmt[:matchlen]
    files = np.sort(os.listdir(direc))
    samplefile = dt.datetime.now().strftime(fmt)
    flen = len(samplefile)
    return [x for x in files if len(x) == flen and x[:matchlen] == match]


def get_filenames_3(direc: str, fmt: str) -> list:
    """
    Get all filenames having the specified length and prefix
    @param direc  Desired directory to get file names from
    @param fmt  Format of desired files from direc
    @return  A list of file names
    """
    files = np.sort(os.listdir(direc))
    samplefile = dt.datetime.now().strftime(fmt)
    flen = len(samplefile)
    return [x for x in files if len(x) == flen and x[:3] == samplefile[:3]]


def get_filenames_dates(direc: str, fmt: str):
    """
    Get all filenames having the specified length and prefix, and datetimes
    @param direc  Desired directory to get file names from
    @param fmt  Format of desired files from direc
    @return  A tuple of file names and datetimes
    """
    files = get_filenames(direc, fmt)
    file_n_date = map(lambda f: (f, dt.datetime.strptime(f, fmt)), files)
    return tuple(file_n_date)


def get_filenames_and_dates(direc: str, fmt: str) -> tuple[list, list]:
    """
    Get all filenames having the specified length and prefix, and datetimes
    @param direc  Desired directory to get file names from
    @param fmt  Format of desired files from direc
    @return  A list of file names
    @return  A list of dates corresponding to the file names
    """
    files = get_filenames(direc, fmt)
    dates = [dt.datetime.strptime(f, fmt) for f in files]
    return files, dates


def getfilenames_daterange(
    direc: str, flen: int, pfix: str, dt1: dt.datetime, dt2: dt.datetime
) -> list:
    """
    Get a list of filenames for dates between start and end date
    @param direc  Directory where files are
    @param flen  Length of file
    @param pfix  File prefix
    @param dt1  Start date
    @param dt2  End date
    """
    plen = len(pfix)
    date1 = int(dt1.strftime("%Y%m%d"))
    date2 = int(dt2.strftime("%Y%m%d"))

    fnames = getfilenames(direc, flen, pfix)
    dates = np.array([float(fname[plen : plen + 8]) for fname in fnames])
    ikeep = np.intersect1d(
        np.where(dates >= date1)[0], np.where(dates <= date2)[0]
    )
    return [direc + fnames[x] for x in ikeep]


def getfiles_for_daterange(
    direc: str, fmt: str, dt1: dt.datetime, dt2: dt.datetime
):
    """
    Get a list of filenames and dates for dates between two values
    @param direc  Directory where files are
    @param fmt  Format of file name
    @param dt1  Desired starting date
    @param dt2  Desired ending date
    @return List of files between start and end dates
    @return List of dates corresponding to files
    """
    filesl = get_filenames(direc, fmt)
    datesl = [
        dt.datetime.strptime(f, fmt).replace(tzinfo=dt.timezone.utc)
        for f in filesl
    ]
    dates = np.array(datesl)
    files = np.array(filesl)
    ikeep = np.where(np.logical_and(dates >= dt1, dates <= dt2))[0]
    return [files[x] for x in ikeep], dates[ikeep]


def getfnames_daterange_utc(
    direc: str, fmt: str, dt1: dt.datetime, dt2: dt.datetime
):
    """
    Get a list of filenames and dates for dates between two values
    assuming UNAWARE datetimes
    @param direc  Directory where files are
    @param fmt  Format of file name
    @param dt1  Desired starting date
    @param dt2  Desired ending date
    @return List of files between start and end dates
    @return List of dates corresponding to files
    """
    filesl = get_filenames(direc, fmt)
    datesl = [dt.datetime.strptime(f, fmt) for f in filesl]
    dates = np.array(datesl)
    files = np.array(filesl)
    ikeep = np.where(np.logical_and(dates >= dt1, dates <= dt2))[0]
    return [files[x] for x in ikeep], dates[ikeep]
