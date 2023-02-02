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
    # TODO: remove [:1] below
    for day in days[:1]:

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


# LATITUDE = esc_params.LATITUDE
# LONGITUDE = esc_params.LONGITUDE
# ALTITUDE = esc_params.ALTITUDE

# DATE1 = esc202202.DATE1
# DATE2 = esc202202.DATE2


# outdir = "/Users/prowe/Sync/measurements/Escudero/era5/"

# north = np.round(LATITUDE / 0.25) * 0.25 + 0.25
# east = np.round(LONGITUDE / 0.25) * 0.25 + 0.25
# latlon = {
#     "north": north,
#     "south": north - 0.25,
#     "east": east,
#     "west": east - 0.5,
# }

# # Pressure levels
# press_levs = np.hstack(
#     [
#         np.arange(1000, 925 - 5, -25),
#         np.arange(900, 200, -50),
#         np.arange(200, 100 - 5, -25),
#     ]
# )
# press_levs_ints = list(press_levs) + [70, 50, 30, 20, 10, 7, 5, 3, 2]
# press_levs = [str(p) for p in press_levs_ints]

# # Times
# times = [f"0{i}:00" for i in range(10)] + [f"{i}:00" for i in range(10, 24)]
# dataset = "reanalysis-era5-pressure-levels"
# download_flag = True

# yearstr = str(DATE1.year)
# monthstr = str(DATE1.month).zfill(2)
# days = list(range(DATE1.day - 1, DATE2.day + 2))  # pad out a day on each side

# prefix = "era5_esc_202202"

# paramsets = get_paramsets(
#     prefix,
#     yearstr,
#     monthstr,
#     days,
#     times,
#     latlon,
#     press_levs,
# )
# get_era5(outdir, paramsets, dataset)
