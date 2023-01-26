#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:44:48 2022

@author: prowe
"""

from netCDF4 import Dataset, num2date
import numpy as np


def read_era5_broadband_down(
    filename: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray,]:
    """
    @param filename  Name of netcdf file
    @return datetimes corresponding to swd and lwd
    @return shortwave downard radiation (swd), W/m2
    @return longwave downward radiation (lwd), W/2
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
