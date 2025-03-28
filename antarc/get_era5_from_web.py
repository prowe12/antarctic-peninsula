#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 13:36:00 2022

@author: prowe
"""

from os.path import exists
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
            "u_component_of_wind",
            "v_component_of_wind",
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
    # Uncomment to add cloud info
    # "specific_cloud_liquid_water_content",
    # "specific_cloud_ice_water_content",
    # "fraction_of_cloud_cover",


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

        outdir0 = f"{outdir}{paramset['params']['year']}/"

        if exists(outdir0 + paramset["outfile"]):
            continue
        # retrieve the path to the file
        fid = c_era.retrieve(dataset, paramset["params"])

        # download the file and save it to the output directory
        fid.download(outdir0 + paramset["outfile"])


def get_paramsets_monthly(
    prefix,
    yearstr,
    monthstr,
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

    paramset = get_paramset_monthly(
        prefix,
        yearstr,
        monthstr,
        latlon,
        press_levs,
    )

    paramsets.append(paramset)
    return paramsets


def get_paramset_monthly(
    prefix,
    yearstr,
    monthstr,
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

    outfile = prefix + ".nc"

    # api parameters:
    # area: [north, west, south, east]; west is negative

    params = {
        "format": "netcdf",
        "product_type": "monthly_averaged_reanalysis",
        "variable": [
            "temperature",
            "relative_humidity",
        ],
        "pressure_level": press_levs,
        "year": yearstr,
        "month": monthstr,
        "time": "00:00",
        "grid": [0.25, 0.25],
        "area": [
            latlon["north"],
            latlon["west"],
            latlon["south"],
            latlon["east"],
        ],
    }
    # Uncomment to add cloud info
    # "specific_cloud_liquid_water_content",
    # "specific_cloud_ice_water_content",
    # "fraction_of_cloud_cover",

    # "variable": [
    #     "pressure",
    #     "Temperature",
    #     "Relative humidity",
    #     "Ozone mass mixing ratio",
    #     "u_component_of_wind",
    #     "v_component_of_wind",
    # ],

    return {"outfile": outfile, "params": params}


# c.retrieve(
#     'reanalysis-era5-pressure-levels-monthly-means',
#     {
#         'format': 'netcdf',
#         'product_type': 'monthly_averaged_reanalysis',
#         'variable': [
#             'fraction_of_cloud_cover', 'relative_humidity', 'specific_cloud_ice_water_content',
#             'specific_cloud_liquid_water_content', 'temperature', 'u_component_of_wind',
#             'v_component_of_wind',
#         ],
#         'pressure_level': [
#             '30', '50', '70',
#             '100', '125', '150',
#             '175', '200', '225',
#             '250', '300', '350',
#             '400', '450', '500',
#             '550', '600', '650',
#             '700', '750', '775',
#             '800', '825', '850',
#             '875', '900', '925',
#             '950', '975', '1000',
#         ],
#         'year': '2018',
#         'month': '01',
#         'time': '00:00',
#         'area': [
#             90, -180, -90,
#             180,
#         ],
#     },
#     'download.nc')
