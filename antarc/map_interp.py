#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 12:31:53 2023

@author: prowe
"""

from scipy import interpolate, ndimage
import numpy as np


def map_interp(amat, date, lat, lon, dates, lats, lons):
    """
    Interpolate to lat/lon/date
    dims of amat must be as follows:
    amat = time, level, lat, lon
    e.g.   (8,    25,   90, 120)
    """

    # Use scipy.interpolate instead of numpy.interp because the former
    # can handle descending data
    flat = interpolate.interp1d(
        lats, np.arange(lats.shape[0]), assume_sorted=False
    )
    flon = interpolate.interp1d(
        lons, np.arange(lons.shape[0]), assume_sorted=False
    )
    ilat = flat(lat)
    ilon = flon(lon)

    ddate = [d - date for d in dates]
    dday = [d.days + d.seconds / 60 / 60 / 24 for d in ddate]
    itime = np.interp(0, dday, np.arange(len(dday)))

    # As arrays
    nlevel = amat.shape[1]
    ilevs = np.arange(nlevel)
    ilats = ilat * np.ones(nlevel)
    ilons = ilon * np.ones(nlevel)
    itimes = itime * np.ones(nlevel)

    dims = [[itimes], [ilevs], [ilats], [ilons]]

    return ndimage.map_coordinates(amat, dims)[0]
