#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 08:56:46 2022

@author: prowe

Purpose: Get the ERA5 monthly broadband

"""

# Dependencies
import datetime as dt
import numpy as np

# My modules
from antarc.escudero.all.read_era5_bbnd import read_era5_broadband_down

# Parameter modules specific to the user: the user must create params.py
# in the directory antarc with the variables MEAS_DIR and PROJ_DIR set
from antarc import params

# Note: ERA5 data from web again for 2023/08/01 was non-final
# as of 2023/10/25 (to get again, see get_cloud_forcing).


def get_era5_broadband_monthly(
    start_date: dt.datetime,
    end_date: dt.datetime,
):
    """
    Get broadband monthly stats between dates for all cases with surface
    temperataures above a threshold
    start_date: starting date
    end_date: ending date
    """

    # Directories
    era_broadband_dir = f"{params.MEAS_DIR}Escudero/era5/broadband/"

    # File formats
    erafmt = "era5_esc_broadband_*[0-9]*.nc"
    era_lat = -62.25
    era_lon = -59.0

    # Load in clear scene files
    date1 = dt.datetime(2017, 1, 1, 0, 0)
    date2 = dt.datetime(2023, 8, 31, 23, 59)

    era_keys = [
        "swd",
        "lwd",
        "swd_clr",
        "lwd_clr",
        "sw_net",
        "lw_net",
        "sw_net_clr",
        "lw_net_clr",
    ]
    era = {"date": np.array([])}
    for key in era_keys:
        era[key] = np.array([])

    years = np.arange(start_date.year, end_date.year + 1)
    for year in years:
        # Get ERA5 data for the entire year
        date1 = dt.datetime(year, 1, 1, 0, 0)
        date2 = dt.datetime(year, 12, 31, 23, 59)

        # If we are in the first or last year, use the specified start/end
        if year == start_date.year:
            date1 = start_date
        if year == end_date.year:
            date2 = end_date

        # Useful variables
        year_str = date1.strftime("%Y")

        # QC on input dates
        if year_str != date2.strftime("%Y"):
            raise ValueError("Different years must be processed separately")

        # Get ERA5 data - recommended to use for cloud forcing
        era5rad_outdir = f"{era_broadband_dir}{date1.strftime('%Y')}/"

        # Read in the ERA5 data
        era5bb = read_era5_broadband_down(era5rad_outdir, erafmt, date1, date2)
        if (
            era_lat not in era5bb["lat"].data
            or era_lon not in era5bb["lon"].data
        ):
            raise ValueError("Bad latititude or longitude")
        # Get closest grid point
        i = np.where(era5bb["lat"] == era_lat)[0][0]
        j = np.where(era5bb["lon"] == era_lon)[0][0]

        era["date"] = np.hstack([era["date"], era5bb["date"].data])
        for key in era_keys:
            era[key] = np.hstack([era[key], era5bb[key][:, i, j].data])

    # Subtract 1/2 hour from ERA5 times
    era["date"] -= dt.timedelta(minutes=30)

    return era
