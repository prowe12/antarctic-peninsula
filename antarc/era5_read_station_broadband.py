#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:44:48 2022

@author: prowe
"""

from netCDF4 import Dataset, num2date
import numpy as np
import datetime as dt
import glob


def read_era5_broadband_down(
    filename: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray,]:
    """
    Read ERA5 downwelling-to-surface broadband fluxes from files created
    by Xun Zou
    @param filename  Name of netcdf file
    @return datetimes corresponding to swd and lwd
    @return shortwave downard radiation (swd), W/m2
    @return longwave downward radiation (lwd), W/m2
    Reminders:
    ncid.variables.keys()
    ncid["ssrd_tosta"]
    """

    ncid = Dataset(filename, "r")

    date = num2date(
        ncid["time"][:].data,
        ncid["time"].units,
        only_use_python_datetimes=True,
        only_use_cftime_datetimes=False,
    )

    return date, ncid["ssrd_tosta"][:].data, ncid["strd_tosta"][:].data


def get_dates_from_filenames(files, dfmt, date1, date2):
    # dates = get_dates_from_filenames(files, '%Y%m%d', date1, date2)
    iyear = files[0].find(str(date1.year))
    idot = files[0].find(".nc")
    dates = [dt.datetime.strptime(f[iyear:idot], dfmt) for f in files]
    dates = [d for d in dates if date1 <= d <= date2]
    return dates


def read_era5_broadband_down_by_file(
    direc: str, fmt: str, date1: dt.datetime, date2: dt.datetime
) -> tuple[np.ndarray, np.ndarray, np.ndarray,]:
    """
    Read ERA5 downwelling-to-surface broadband fluxes from files download
    using scripts in this code repo
    @param filename  Name of netcdf file
    @return datetimes corresponding to swd and lwd
    @return shortwave downard radiation (swd), W/m2
    @return longwave downward radiation (lwd), W/m2
    Reminders:
    ncid.variables.keys()
    ncid["ssrd_tosta"]
    """

    if date1.year != date2.year:
        raise ValueError("ERA5 dates must all fall within a year")

    # Get all filenames in between date1 and date2
    fmt = direc + "era5_esc_broadband_*[0-9]*.nc"
    files = glob.glob(fmt)

    if len(files) == 0:
        return None

    islash = files[0].rfind("/")
    files = [f[islash + 1 :] for f in files]

    iyear = files[0].find(str(date1.year))
    idot = files[0].find(".nc")
    dfmt = "%Y%m%d"

    erafiles = []
    for fname in np.sort(files):
        thisdate = dt.datetime.strptime(fname[iyear:idot], dfmt)
        thisdate = thisdate.replace(tzinfo=dt.timezone.utc)
        if date1 <= thisdate <= date2:
            erafiles.append(fname)

    era5 = {}
    init = True
    for fname in erafiles:
        with Dataset(direc + fname, "r") as nci:
            era_date = num2date(
                nci["time"][:].data,
                nci["time"].units,
                only_use_python_datetimes=True,
                only_use_cftime_datetimes=False,
            )
            if init:
                era5["dates"] = era_date
                era5["swd"] = nci["ssrd"][:]
                era5["lwd"] = nci["strd"][:]
                era5["swd_clr"] = nci["ssrdc"][:]
                era5["lwd_clr"] = nci["strdc"][:]
                era5["lat"] = nci["latitude"][:]
                era5["lon"] = nci["longitude"][:]
                init = False
            else:
                era5["dates"] = np.hstack([era5["dates"], era_date])
                era5["swd"] = np.vstack([era5["swd"], nci["ssrd"][:]])
                era5["lwd"] = np.vstack([era5["lwd"], nci["strd"][:]])
                era5["swd_clr"] = np.vstack([era5["swd_clr"], nci["ssrdc"][:]])
                era5["lwd_clr"] = np.vstack([era5["lwd_clr"], nci["strdc"][:]])

            # if "dates" not in era5:
            #     era5["dates"] = era_date
            # else:
            #     era5["dates"] = np.hstack([era5["dates"], era_date])
            # if "swd" not in era5:
            #     era5["swd"] = nci["ssrd"][:]
            # else:
            #     era5["swd"] = np.vstack([era5["swd"], nci["ssrd"][:]])
            # if "lwd" not in era5:
            #     era5["lwd"] = nci["strd"][:]
            # else:
            #     era5["lwd"] = np.vstack([era5["lwd"], nci["strd"][:]])

    era5["swd"] = era5["swd"] / 3600
    era5["lwd"] = era5["lwd"] / 3600
    era5["swd_clr"] /= 3600
    era5["lwd_clr"] /= 3600
    return era5
