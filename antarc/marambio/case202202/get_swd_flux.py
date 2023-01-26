#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 16:22:16 2019

@author: prowe
"""

# Dependencies
import numpy as np
import pysolar
import matplotlib.pyplot as plt

# My modules
from antarc.load_rad_flux import load_rad_flux
from antarc.run_libradtran import run_libradtran

# Parameter modules
from antarc.escudero.parameters import esc_params, esc202202, swd_params


def get_swd() -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
]:
    """
    Get shortwave downward radiation and forcing betwen start and end dates
    @return  Dates between start and end dates
    @return   Measured shortwave downward radiation
    @return  Simulated clear sky downward radiation
    @return  Diffuse radiation according to pysolar
    @return  Factor for diffuse radiation
    """
    # Location and time-dependent parameters
    lat = esc_params.LATITUDE
    lon = esc_params.LONGITUDE
    pyr_dir = swd_params.STAND_DIR
    pyr_fmt = swd_params.PYR_FILEFORMAT
    pfac = swd_params.SW_FACTOR

    start_date = esc202202.DATE1
    end_date = esc202202.DATE2

    # Get measured swd and dates
    pyr_dir += str(start_date.year) + "/"
    dtimes, rad = load_rad_flux(pyr_dir, pyr_fmt, start_date, end_date)
    dtimes = np.array(dtimes)

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

    plotfigs = False
    if plotfigs:
        # .. Plot it
        plt.figure(1)
        plt.clf()
        plt.plot(dtimes, rad, label="pyranometer")
        plt.plot(dtimes, diffuse, label="pysolar")
        plt.plot(dtimes, diffuse * diffuse_fac, label="clear")

        # .. Cloud forcing
        sw_forcing = diffuse * diffuse_fac - rad
        plt.plot(dtimes, sw_forcing, label="forcing")
        plt.legend()

    # Get clear-sky radiation from libradtran for comparison
    albedo = 0.2
    ozone = 257.0

    first = True
    librad = np.array([])
    ndates = len(dtimes)
    night_result = run_libradtran(1, 90.0, albedo, ozone)
    # This is very slow, so run every xth point
    for i in range(0, ndates, 50):  # ):
        thisdate = dtimes[i]

        dayofyear = thisdate.timetuple().tm_yday

        sza0 = 90 - sunalt[i]
        if sza0 > 90:
            result = np.array([0 * x for x in night_result])
        else:
            result = run_libradtran(dayofyear, sza0, albedo, ozone)

        if first:
            librad = result
            first = False
        else:
            librad = np.vstack([librad, result])

    return dtimes, rad, librad, diffuse * diffuse_fac, diffuse, diffuse_fac
