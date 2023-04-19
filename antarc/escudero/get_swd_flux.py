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
    Notes: average given by:
        (268.5+260.9+240.9+243.3+275.5+256.1)/6
        2021/02/05 - 2021/02/08 OMI
    """
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
