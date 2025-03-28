#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:44:48 2022

@author: prowe
"""

# Modules
from netCDF4 import Dataset, num2date
import numpy as np
import datetime as dt
import glob

# My modules
from antarc.map_interp import map_interp


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
    """
    Assuming the files have the format dfmt, return the dates
    that are part of each file name within date1 and date2
    @param files  Filenames
    @param dfmt  Format of each filename
    @param date1  Starting datetime for returning dates
    @param date2  Ending datetime for returning dates
    @returns  A list of datetimes
    """
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
    # Get all filenames in between date1 and date2
    fmt = direc + "era5_esc_broadband_*[0-9]*.nc"
    files = glob.glob(fmt)
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

    era5 = {"dates": []}
    for fname in erafiles:
        with Dataset(direc + fname, "r") as nci:
            era_date = num2date(
                nci["time"][:].data,
                nci["time"].units,
                only_use_python_datetimes=True,
                only_use_cftime_datetimes=False,
            )
            if "dates" not in era5:
                era5["dates"] = era_date
            else:
                era5["dates"] = np.hstack([era5["dates"], era_date])
            if "lat" not in era5:
                era5["lat"] = nci["latitude"][:]
            if "lon" not in era5:
                era5["lon"] = nci["longitude"][:]
            if "swd" not in era5:
                era5["swd"] = nci["ssrd"][:]
            else:
                era5["swd"] = np.vstack([era5["swd"], nci["ssrd"][:]])
            if "lwd" not in era5:
                era5["lwd"] = nci["strd"][:]
            else:
                era5["lwd"] = np.vstack([era5["lwd"], nci["strd"][:]])

    era5["swd"] = era5["swd"] / 3600
    era5["lwd"] = era5["lwd"] / 3600
    return era5


def read_era5_profiles(fname: str) -> dict[float]:
    """
    Read ERA5 profiles from downloaded file
    @param fname  Name of netcdf file
    @returns  A dictionary containing the ERA5 data
    """
    era5 = {}
    with Dataset(fname, "r") as nci:
        era5["date"] = num2date(
            nci["time"][:].data,
            nci["time"].units,
            only_use_python_datetimes=True,
            only_use_cftime_datetimes=False,
        )
        era5["lat"] = nci["latitude"][:]
        era5["lon"] = nci["longitude"][:]
        era5["level"] = nci["level"][:]
        era5["time"] = nci["time"][:]
        era5["alt"] = nci["z"][:]
        era5["temp"] = nci["t"][:]
        era5["rhw"] = nci["r"][:]
        era5["ozone"] = nci["o3"][:]
        era5["uwind"] = nci["u"][:]
        era5["vwind"] = nci["v"][:]

    return era5


def read_era5_profiles_interp(fname: str, lat: float, lon: float, date):
    """
    Interpolate the ERA data to the given lat, lon, and datetime
    @param era  A dictionary of the era data
    @param lat  The latitude
    @param lon  The longitude
    @returns  A new dictionary containing the interpolated data
    """
    lats, lons, dates, press, zmat, tmat, rhmat, o3mat, umat, vmat = load_era(
        fname
    )

    era = {"press": press}

    ozone = map_interp(o3mat, date, lat, lon, dates, lats, lons)
    era["temp"] = map_interp(tmat, date, lat, lon, dates, lats, lons)
    era["rhw"] = map_interp(rhmat, date, lat, lon, dates, lats, lons)
    era["alt"] = map_interp(zmat, date, lat, lon, dates, lats, lons)
    era["ozone"] = ozone * 1000  # kg/kg * 1000 g/kg => g/kg

    # For now, allow for missing winds
    if "uwind" in era.keys() and "vwind" in era.keys():
        era["uwind"] = map_interp(umat, date, lat, lon, dates, lats, lons)
        era["vwind"] = map_interp(vmat, date, lat, lon, dates, lats, lons)
    else:
        era["uswind"] = np.nan + np.zeros(len(era["temp"]))
        era["vwind"] = np.nan + np.zeros(len(era["temp"]))

    return era


def load_era(fname: str):
    """
    Load ERA5 data
    @param  fname  ERA5 file name
    """
    with Dataset(fname, "r", format="NETCDF4") as nci:
        # Variables, for convenience
        # 'longitude', 'latitude', 'level', 'time', 't', 'r'
        # Dimensions for t, r, o3, z:
        # (time, level, latitude, longitude)
        lats = nci["latitude"][:]
        lons = nci["longitude"][:]
        dates = num2date(
            nci["time"][:], nci["time"].units, nci["time"].calendar
        )

        # Check units
        # Expected variables and units
        expectedvars = {
            "level": "millibars",
            "t": "K",
            "r": "%",
            "z": "m**2 s**-2",
            "o3": "kg kg**-1",
            "u": "m s**-1",
            "v": "m s**-1",
        }

        # For now, allow for case where u/v are missing
        # (but print file names so we can replace them)
        for exptd in expectedvars:
            if exptd not in nci.variables.keys() and exptd in ["u", "v"]:
                print(f"{fname}: {exptd} missing")
            elif exptd not in nci.variables.keys() and exptd in ["u", "v"]:
                raise ValueError(f"{fname}: variable {exptd} missing")
            elif nci[exptd].units != expectedvars[exptd]:
                msg1 = "Expected units of " + expectedvars[exptd] + " for "
                msg2 = exptd + ", but got " + nci[exptd].units
                raise ValueError(msg1 + msg2)

        press = nci["level"][:]
        zmat = nci["z"][:] / 9.80665 / 1000
        tmat = nci["t"][:]
        rhmat = nci["r"][:]
        o3mat = nci["o3"][:]

        # Catch for missing u/v
        if "u" in nci.variables.keys():
            umat = nci["u"][:]
            vmat = nci["v"][:]
        else:
            umat = None
            vmat = None

        return lats, lons, dates, press, zmat, tmat, rhmat, o3mat, umat, vmat


def read_era5_profiles_monthly(fname: str, lat: float, lon: float, date):
    """
    Interpolate the ERA data to the given lat, lon, and datetime
    @param era  A dictionary of the era data
    @param lat  The latitude
    @param lon  The longitude
    @returns  A new dictionary containing the interpolated data
    """
    lats, lons, dates, press, tmat, rhmat = load_era_monthly(fname)

    era = {"press": press}

    era["temp"] = map_interp(tmat, date, lat, lon, dates, lats, lons)
    era["rhw"] = map_interp(rhmat, date, lat, lon, dates, lats, lons)

    return era


def load_era_monthly(fname: str):
    """
    Load ERA5 data
    @param  fname  ERA5 file name
    """
    with Dataset(fname, "r", format="NETCDF4") as nci:
        # Variables, for convenience
        # 'longitude', 'latitude', 'level', 'time', 't', 'r'
        # Dimensions for t, r, o3, z:
        # (time, level, latitude, longitude)
        lats = nci["latitude"][:]
        lons = nci["longitude"][:]
        dates = num2date(
            nci["time"][:], nci["time"].units, nci["time"].calendar
        )

        # Check units
        # Expected variables and units
        expectedvars = {
            "level": "millibars",
            "t": "K",
            "r": "%",
        }

        # For now, allow for case where u/v are missing
        # (but print file names so we can replace them)
        for exptd in expectedvars:
            if nci[exptd].units != expectedvars[exptd]:
                msg1 = "Expected units of " + expectedvars[exptd] + " for "
                msg2 = exptd + ", but got " + nci[exptd].units
                raise ValueError(msg1 + msg2)

        press = nci["level"][:]
        tmat = nci["t"][:]
        rhmat = nci["r"][:]

        return lats, lons, dates, press, tmat, rhmat
