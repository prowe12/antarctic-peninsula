#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  26 09:33:02 2023

@author: prowe

ssrd
-----
This parameter is the amount of solar radiation (shortwave radiation) that 
reaches a horizontal plane at the surface of the Earth. This parameter 
comprises both direct and diffuse solar radiation.

Radiation from the Sun (solar, or shortwave, radiation) is partly reflected 
back to space by clouds and particles in the atmosphere (aerosols) and some 
of it is absorbed. The rest is incident on the Earth's surface (represented 
by this parameter). See further documentation.

To a reasonably good approximation, this parameter is the model equivalent 
of what would be measured by a pyranometer (an instrument used for measuring 
solar radiation) at the surface. However, care should be taken when comparing 
model parameters with observations, because observations are often local to 
a particular point in space and time, rather than representing averages over 
a model grid box.

This parameter is accumulated over a particular time period which depends 
on the data extracted. The units are joules per square metre (J m-2). 
To convert to watts per square metre (W m-2), the accumulated values should 
be divided by the accumulation period expressed in seconds. The ECMWF 
convention for vertical fluxes is positive downwards.
"""

import glob
import datetime as dt
from netCDF4 import Dataset, num2date
import numpy as np


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

    if iyear == -1:
        raise ValueError(f"First file does not include {date1.year}")

    erafiles = []
    for fname in np.sort(files):
        thisdate = dt.datetime.strptime(fname[iyear:idot], dfmt)
        thisdate = thisdate.replace(tzinfo=dt.timezone.utc)
        if date1 <= thisdate <= date2:
            erafiles.append(fname)

    # ssrdc: Surface solar radiation downward clear-sky
    # ssrc: Surface net short-wave (solar) radiation, clear sky
    # strd: Surface long-wave (thermal) radiation downwards
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

            if "swd_clr" not in era5:
                era5["swd_clr"] = nci["ssrdc"][:]
            else:
                era5["swd_clr"] = np.vstack([era5["swd_clr"], nci["ssrdc"][:]])

            if "lwd_clr" not in era5:
                era5["lwd_clr"] = nci["strdc"][:]
            else:
                era5["lwd_clr"] = np.vstack([era5["lwd_clr"], nci["strdc"][:]])

    era5["swd"] = era5["swd"] / 3600
    era5["lwd"] = era5["lwd"] / 3600
    era5["swd_clr"] = era5["swd_clr"] / 3600
    era5["lwd_clr"] = era5["lwd_clr"] / 3600

    return era5
