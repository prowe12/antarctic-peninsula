#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 13:36:00 2022

@author: prowe
"""

import cdsapi


def get_params(datestr: list[str], times, latlon, press_levs):
    """
    Get the ERA5 data as indicated in params from the web
    @returns  A list of the parameter structures
    @returns  The output file names
    """
    # api parameters:
    # area: [north, west, south, east]; west is negative
    return {
        "format": "netcdf",
        "product_type": "reanalysis",
        "variable": [
            "pressure",
            "Geopotential",
            "Temperature",
            "Relative humidity",
            "Ozone mass mixing ratio",
        ],
        "pressure_level": press_levs,
        "year": [datestr[0]],
        "month": [datestr[1]],
        "day": [datestr[2]],
        "time": times,
        "grid": [0.25, 0.25],
        "area": [
            latlon["north"],
            latlon["west"],
            latlon["south"],
            latlon["east"],
        ],
    }


def get_outfile(prefix: str, datestring: str, ext: str = ".nc"):
    """
    Get the name of the output file
    @param  prefix  Desired beginging of file name
    @param  datestring  Date string
    @param  ext  Extension, beginning with '.'; e.g. '.nc'
    """
    return prefix + datestring + ext


def get_paramsets(
    prefix,
    yearstr,
    monthstr,
    days,
    times,
    latlon,
    press_levs,
):
    """
    Get dictionaries with the output file names and the parameters
    dictionary of parameters needed for retrieving ERA5 data
    and pack it together into
    @returns  A list of the parameter structures
    @returns  The output file names
    """

    paramsets = []
    for day in days:

        paramset = get_paramset(
            prefix,
            yearstr,
            monthstr,
            day,
            times,
            latlon,
            press_levs,
        )

        paramsets.append(paramset)
    return paramsets


def get_paramset(
    prefix,
    yearstr,
    monthstr,
    day,
    times,
    latlon,
    press_levs,
):
    """
    Get dictionaries with the output file names and the parameters
    dictionary of parameters needed for retrieving ERA5 data
    and pack it together into
    @returns  A list of the parameter structures
    @returns  The output file names
    """

    daystr = str(day).zfill(2)
    datestring = daystr
    outfile = get_outfile(prefix, datestring, ".nc")

    # api parameters:
    # area: [north, west, south, east]; west is negative
    params = get_params(
        [yearstr, monthstr, daystr],
        times,
        latlon,
        press_levs,
    )

    return {"outfile": outfile, "params": params}


def get_era5(outdir, paramsets, dataset):
    """
    Get the ERA5 data as indicated in params from the web
    """
    c_era = cdsapi.Client()
    for paramset in paramsets:

        # retrieve the path to the file
        fid = c_era.retrieve(dataset, paramset["params"])

        # download the file and save it to the output directory
        fid.download(outdir + paramset["outfile"])
