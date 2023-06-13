#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 16:22:16 2019

@author: prowe
"""

# Dependencies
import numpy as np
import pysolar
from os.path import exists
import pandas as pd
import datetime as dt

# My modules
from antarc.load_rad_flux import load_rad_flux
from antarc.run_libradtran import run_libradtran


def get_obs_swd(swd_params, esc_case) -> tuple[np.ndarray, np.ndarray]:
    """
    Get shortwave downward radiation and forcing betwen start and end dates
    @return   Dates between start and end dates
    @return   Measured shortwave downward radiation
    """
    # Location and time-dependent parameters
    pyr_dir = swd_params.SWD_DIR
    pyr_fmt = swd_params.SWD_FILEFORMAT
    start_date = esc_case.DATE1
    end_date = esc_case.DATE2

    # Get measured swd and dates
    pyr_dir += str(start_date.year) + "/"
    dtimes, rad = load_rad_flux(pyr_dir, pyr_fmt, start_date, end_date)
    dtimes = np.array(dtimes)

    return dtimes, rad


def get_pysolar_swd(
    swd_params, esc_params, dtimes
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get shortwave downward radiation and forcing betwen start and end dates
    @param  dtimes: Dates between start and end dates
    @return   Simulated clear sky downward radiation
    @return   Diffuse radiation according to pysolar
    @return   Factor for diffuse radiation
    """

    # Location and time-dependent parameters
    lat = esc_params.LATITUDE
    lon = esc_params.LONGITUDE
    pfac = swd_params.SW_FACTOR

    # Get diffuse radiation from pysolar
    ndates = len(dtimes)
    diffuse = 0 * np.zeros(ndates)
    sunalt = 0 * np.zeros(ndates)
    for i in range(ndates):
        thisdate = dtimes[i]
        altitude_deg = pysolar.solar.get_altitude(lat, lon, thisdate)
        if altitude_deg > 0:
            diffuse[i] = pysolar.util.diffuse_underclear(lat, lon, thisdate)
        sunalt[i] = altitude_deg

    # Get factor for scaling pysolar diffuse radiance to get clear sky rad
    diffuse_fac = np.polyval(pfac, sunalt)

    return diffuse * diffuse_fac, diffuse, diffuse_fac


def get_albedo(dayofyear) -> float:
    """
    Get the albedo based on the current datetime
    @param dtime  The current datetime
    @return  The albedo
    """
    # fmt: off
    day = [
        1,    13,  25,  37,  49,  61,  73,  85,  97, 109, 121, 133, 145, 157, 
        169, 181, 193, 205, 217, 229, 241, 253, 265, 277, 289, 301, 313, 325, 
        337, 349, 361, 365, 
    ]
    albedo = [
        0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.25, 0.3, 0.4, 
        0.5, 0.55, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.55, 0.5, 0.4, 
        0.3, 0.25, 0.2, 0.2, 0.2, 0.2, 
    ]
    # fmt: on
    return np.interp(dayofyear, day, albedo)


def get_ozone(dtime) -> float:
    """
    Get the ozone amount based on the current datetime
    @param dtime  The current datetime
    @return  The ozone amount
    From https://www.temis.nl/protocols/o3col/overpass_omi.php
    Notes: average for 2021/02/05 - 2021/02/08 OMI given by:
        (268.5+260.9+240.9+243.3+275.5+256.1)/6

    2020/05/10-20:
    20200510 181240 -62.22  -59.10  11 84.55  279.2    3.8  21019  0.77  411.0 16  2  -71.1  24.9  31.6
    20200511 171721 -62.16  -59.13  12 81.90  289.2    4.0  20320  1.00  665.0  3  2  -73.1  13.2  63.3
    20200512 180024 -62.19  -58.89   1 84.31  288.3    3.0  22033  0.92  668.0 12  2  -70.3  12.7  40.6
    20200513 170500 -62.22  -58.30  32 82.10  259.3    4.3  18891  1.00  793.0  2  2  -75.3   7.9  66.4
    20200514 174806 -62.18  -58.48  22 84.17  270.5    4.3  20753  0.69  829.0  9  2  -69.5   6.3  47.7
    20200515 165238 -62.21  -57.99  48 82.17  287.7    7.6  21830  1.00  916.0  1  2  -75.2   3.2  69.7
    20200516 173549 -62.18  -58.84   3 83.92  292.8    3.2  22277  0.81  623.0  6  2  -66.3  15.0  55.1
    20200517 181856 -61.71  -59.66  66 85.93  287.4    4.9  23693  0.27  264.0 18  2  -58.9  36.0  27.2
    20200518 172330 -62.20  -58.79   6 83.81  282.3    3.4  21959  0.62  866.0  4  2  -69.8   4.6  60.5
    20200519 180635 -61.97  -59.02  24 85.94  286.4    5.0  24459  0.61  648.0 14  2  -73.3  13.9  36.0
    20200520 171110 -62.17  -59.83  48 83.55  296.0    3.4  22802  1.00  516.0  2  2  -66.8  19.9  66.4
    average is
    """
    # month = np.arange(2, 5)
    dtmax = dt.datetime(2022, 5, 18, tzinfo=dt.timezone.utc)
    dtmin = dt.datetime(2022, 5, 9, tzinfo=dt.timezone.utc)
    if max(dtime) < dtmax and min(dtime) >= dtmin:
        # fmt: off
        vals = [279.2, 289.2, 288.3, 259.3, 270.5, 287.7, 292.8, 287.4, 282.3, 286.4, 296.0]
        # fmt: on
        return np.mean(vals)
    else:
        return 257.53


def calc_libradtran(esc_params, dtimes, outfile) -> np.ndarray:
    """
    Get shortwave downward radiation and forcing betwen start and end dates
    from libradtran and save to file
    @param  dtimes: Dates between start and end dates
    @return  libradtran result
    """

    # Location and time-dependent parameters
    lat = esc_params.LATITUDE
    lon = esc_params.LONGITUDE

    ndates = len(dtimes)
    inds = list(range(0, ndates, 100))  # step by 100?
    librad = np.zeros([len(inds), 6])
    lastdayofyear = -1
    count = 0
    # Open the file that will be written to
    with open(outfile, "w") as fid:
        # Write header to file
        content = "Date, swd, o2, o3, o4, o5, o6,\n"
        fid.write(content)

        for i in inds:
            # Current datetime and day of year
            thisdate = dtimes[i]
            dayofyear = thisdate.timetuple().tm_yday
            dtimestr = thisdate.strftime("%Y%m%d %H:%M:%S")

            # Albedo
            if dayofyear != lastdayofyear:
                albedo = get_albedo(dayofyear)
            lastdayofyear = dayofyear

            # Ozone
            ozone = get_ozone(dtimes)

            # SZA
            sunalt = pysolar.solar.get_altitude(lat, lon, thisdate)
            sza0 = 90 - sunalt

            # Run LIBRADTRAN
            if sza0 > 120:
                result = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            else:
                result = run_libradtran(dayofyear, sza0, albedo, ozone)

            # Tack on the libradtran result
            librad[count, :] = result
            count += 1

            # Write to file, since this is so slow
            content = (
                dtimestr
                + f", {result[0]}, {result[1]}, {result[2]}, {result[3]}, {result[4]}, {result[5]},\n"
            )
            fid.write(content)

    return dtimes[inds], librad


def get_libradtran(params, swd_date, clrfile, redo):
    """
    Get the libradtran results. If redo set to true or file does not exist,
    recompute the results using libradtran
    @param swd_date  The dates
    @param clrfile  The file to load or save results to
    @param redo  Whether to redo the calculation if the file exists.
    """
    if not exists(clrfile) or redo:
        libdate, libclear = calc_libradtran(params, swd_date, clrfile)
        swd_clear = libclear[:, 0]
    else:
        # Load the libradtran results, which have columns: Date,swd,o2,o3,o4,o5,o6
        libclear = pd.read_csv(clrfile, delimiter=",\s")
        libdate = pd.to_datetime(libclear["Date"])
        swd_clear = libclear["swd"]

    return libdate, swd_clear
