#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:03:39 2024

@author: prowe

Get ERA5 cloud base heights

ERA5 cloud base height (https://codes.ecmwf.int/grib/param-db/228023)
Name: Cloud base height
Short name: cbh
Unit: m
ID: 228023
Description: The height above the Earth's surface of the base of the lowest 
cloud layer, at the specified time. This parameter is calculated by searching 
from the second lowest model level upwards, to the height of the level where 
cloud fraction becomes greater than 1% and condensate content greater than 
1.E-6 kg kg-1. Fog (i.e., cloud in the lowest model layer) is not considered 
when defining cloud base height.

"""

from netCDF4 import Dataset, num2date
import numpy as np


from antarc.params import MEAS_DIR
from antarc.escudero.parameters.esc_params import LATITUDE, LONGITUDE


def get_era5_cbh(season: str):
    """
    Get cloud base height (cbh) for ERA5 data for season, at the nearest
    gridpoint to Escudero station.
    season = "All seasons"  # "winter"  # "All seasons"
    Return the date and the cbh, with NaN for cases with no cbh specified.
    """
    measdir = MEAS_DIR + "Escudero/"
    era5dir = measdir + "era5/clouds/"

    # Select years and season desired
    years = [2017, 2018, 2019, 2020, 2021, 2022]
    months_by_season = {
        "summer": [11, 0, 1],
        "fall": [2, 3, 4],
        "winter": [5, 6, 7],
        "spring": [8, 9, 10],
    }

    if season == "All seasons":
        for_months = False
    else:
        for_months = True

    esc_lat = LATITUDE  # negative for south
    esc_lon = LONGITUDE  # negative for west

    if for_months:
        months = months_by_season[season]

    era_cbh = np.array([])
    dtime_cbh = np.array([])
    for year in years:

        era5file = "cloud_" + str(year) + ".nc"

        with Dataset(era5dir + era5file) as era5:
            dtm = num2date(
                era5["time"][:].data,
                era5["time"].units,
                only_use_python_datetimes=True,
                only_use_cftime_datetimes=False,
            )
            iera5_lat = np.argmin(np.abs(era5["latitude"][:].data - esc_lat))
            iera5_lon = np.argmin(np.abs(era5["longitude"][:].data - esc_lon))
            cbh = era5["cbh"][:, iera5_lat, iera5_lon]
            cbh[cbh.mask] = np.nan

            # Include only the months of interest
            if for_months:
                cbh0 = [x for i, x in enumerate(cbh) if dtm[i].month in months]
                dtm0 = [x for i, x in enumerate(dtm) if dtm[i].month in months]
            else:
                cbh0 = [x for x in cbh]
                dtm0 = [x for x in dtm]

        dtime_cbh = np.hstack([dtime_cbh, dtm0])
        era_cbh = np.hstack([era_cbh, cbh0])

    return dtime_cbh, era_cbh
