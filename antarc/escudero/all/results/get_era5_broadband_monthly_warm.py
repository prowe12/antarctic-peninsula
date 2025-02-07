#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 08:56:46 2022

@author: prowe

Purpose: Plot cloud forcing at Escudero for 2017-2023

"""

# Dependencies
import datetime as dt
import numpy as np

# My modules
# Parameter modules specific to the user: the user must create params.py
# in the directory antarc with the variables MEAS_DIR and PROJ_DIR set
from antarc.escudero.get_frei_15min import get_frei_15
from antarc.escudero.all.results.get_era5_broadband_monthly import (
    get_era5_broadband_monthly,
)

from antarc import params


def get_era5_broadband_monthly_warm(
    start_date: dt.datetime,
    end_date: dt.datetime,
    temperature_threshold: float,
):
    """
    start_date = dt.datetime(2021, 9, 1, 0, 0)
    end_date = dt.datetime(2023, 8, 31, 23, 59)
    temperature_threshold = 2
    Return  ERA5 monthly values for surface temperatures above the threshold
    Return  ERA5 monthly values for all data within date range
    """

    # Directories
    frei_direc = f"{params.MEAS_DIR}Escudero/frei_met/by_month_every_15min/"
    frei_sampfname = "esc_met_201803.txt"
    frei_pfx = 7

    # Load in the surface temperatures from FREI
    (
        fdtime,
        ftemp,
        fpress,
        fslp,
        fwdir,
        fwspd,
        frh,
    ) = get_frei_15(frei_direc, frei_sampfname, frei_pfx)

    era = get_era5_broadband_monthly(start_date, end_date)

    # Add the 1/2 hour back to era date in prep for analysis
    era["date"] += dt.timedelta(minutes=30)

    # Loop over era dates and find all hours for which max surface temp > 0 C
    if fdtime[0] > era["date"][0]:
        raise ValueError("First surface temp must be before first ERA5 date")

    # Add 0.9 K to the Frei temperatures to account for the lower altitude
    # of Escudero
    ftemp += 0.9

    # Find the maximum temperature within each hour at Frei
    # Because ERA5 accumulates to the hour, for the first hour
    # just use the first measurement, while for all the rest
    # average over first through last inclusive
    i = 0
    maxtemp = np.nan * np.ones(len(era["date"]))
    for iera, dtime0 in enumerate(era["date"]):
        # Skip to the beginning of the hour
        while fdtime[i] < dtime0 - dt.timedelta(hours=1):
            i += 1
        # Get max time over next hour
        while fdtime[i] <= dtime0:
            i += 1
            maxtemp[iera] = np.nanmax([maxtemp[iera], ftemp[i]])
        # Dial back 1 because we double count the last/first
        i -= 1

    # Subtract 1/2 hour from ERA5 times
    era["date"] -= dt.timedelta(minutes=30)

    # Only use times with max temp over some temperature threshold (in C)
    iwarm = np.where(maxtemp > temperature_threshold)[0]
    era_warm = {}
    for key in era.keys():
        era_warm[key] = era[key][iwarm]

    return era_warm, era
