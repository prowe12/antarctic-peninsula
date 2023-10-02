#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:38:16 2019

Copyright 2022 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

Return the columns for a given set of GRAW files

@author: prowe
"""


def graw_get_column_map(cols, units):
    """
    Various possible Headers for GRAW raw radiosonde data from simulations
    """

    # Expected columns of the file
    cols1 = [
        "Time",  # 0
        "P",  # 1
        "T",  # 2
        "Hu",  # 3
        "Ws",  # 4
        "Wd",  # 5
        "Long.",  # 6
        "Lat.",  # 7
        "Alt",  # 8
        "Geopot",  # 9
        "MRI",  # 10
        "RI",  # 11
        "Dewp.",  # 12
        "Virt. Temp",
        "Rs",
        "Elevation",
        "Azimuth",
        "Range",
        "D",
    ]

    # Expected units of the file
    # MRI and RI do not have units
    units1 = [
        "[sec]",  # 0, time
        "[hPa]",  # 1, press
        "[?C]",  # 2, temp
        "[%]",  # 3, rhw
        "[m/s]",  # 4, wspd
        "[?]",  # 5, wdir
        "[?]",  # 6, long
        "[?]",  # 7, lat
        "[m]",  # 8, alt
        "[m]",  # 9, geopot
        "[?C]",  # 10, dewpt
        "[?C]",  # 11
        "[m/s]",  # 12
        "[?]",
        "[?]",
        "[m]",
        "[kg/m3]",
    ]

    # Map the neccessary profile variables to the columns of the file
    # the order is cols, units.  That is, the first index maps to the
    # cols list, the second to the units list. Two sets are needed because
    # units are missing for some columns
    colmap1 = {
        "time": [0, 0],
        "press": [1, 1],
        "alt": [8, 8],
        "temp": [2, 2],
        "rhw": [3, 3],
        "dewpt": [12, 10],
        "wspd": [4, 4],
        "wdir": [5, 5],
    }

    # Expected columns of the file
    cols2 = [
        "Time",  # 0
        "P",  # 1
        "T",  # 2
        "Hu",  # 3
        "Ws",  # 4
        "Wd",  # 5
        "Long.",  # 6
        "Lat.",  # 7
        "Alt",  # 8
        "Geopot",  # 9
        "MRI",  # 10
        "RI",  # 11
        "Dewp.",  # 12
        "Virt. Temp",
        "Rs",
        "Elevation",
        "Azimuth",
        "Range",
        "D",
    ]

    # Expected units of the file
    # MRI and RI do not have units
    units2 = [
        "[sec]",  # 0, time
        "[hPa]",  # 1, press
        "[°C]",  # 2, temp
        "[%]",  # 3, rhw
        "[kn]",  # 4, wspd
        "[°]",  # 5, wdir
        "[°]",  # 6, long
        "[°]",  # 7, lat
        "[m]",  # 8, alt
        "[m]",  # 9, geopot
        "[°C]",  # 10, dewpt
        "[°C]",  # 11
        "[m/s]",  # 12
        "[°]",
        "[°]",
        "[m]",
        "[kg/m3]",
    ]

    # Map the neccessary profile variables to the columns of the file
    # the order is cols, units.  That is, the first index maps to the
    # cols list, the second to the units list. Two sets are needed because
    # units are missing for some columns
    colmap2 = {
        "time": [0, 0],
        "press": [1, 1],
        "alt": [8, 8],
        "temp": [2, 2],
        "rhw": [3, 3],
        "dewpt": [12, 10],
        "wspd": [4, 4],
        "wdir": [5, 5],
    }

    # Expected columns of the file
    cols3 = [
        "Time",  # 0
        "P",  # 1
        "T",  # 2
        "Hu",  # 3
        "Ws",  # 4
        "Wd",  # 5
        "Long.",  # 6
        "Lat.",  # 7
        "Alt",  # 8
        "Geopot",  # 9
        "MRI",  # 10
        "RI",  # 11
        "Dewp.",  # 12
        "Virt. Temp",
        "Rs",
        "D",
    ]

    # Expected units of the file
    # MRI and RI do not have units
    units3 = [
        "[sec]",  # 0
        "[hPa]",  # 1
        "[°C]",  # 2
        "[%]",  # 3
        "[m/s]",  # 4
        "[°]",  # 5
        "[°]",  # 6
        "[°]",  # 7
        "[m]",  # 8
        "[m]",  # 9
        "[°C]",  # 10, dewpt
        "[°C]",  # 11, virt temp (not used)
        "[m/s]",  # 12, Rs (not used)
        "[kg/m3]",  # 13, D (not used)
    ]

    # Map the neccessary profile variables to the columns of the file
    colmap3 = {
        "time": [0, 0],
        "press": [1, 1],
        "alt": [8, 8],
        "temp": [2, 2],
        "rhw": [3, 3],
        "dewpt": [12, 10],
        "wspd": [4, 4],
        "wdir": [5, 5],
    }

    # Expected columns of the file
    cols4 = [
        "Time",  # 0
        "P",  # 1
        "T",  # 2
        "Hu",  # 3
        "Ws",  # 4
        "Wd",  # 5
        "Long.",  # 6
        "Lat.",  # 7
        "Alt",  # 8
        "Geopot",  # 9
        "O3",  # 10
        "I",  # 11
        "Ti",  # 12
        "TO3",  # 13
        "Dewp.",  # 14
        "Virt. Temp",
        "Rs",
    ]

    # Expected units of the file
    units4 = [
        "[sec]",  # 0
        "[hPa]",  # 1
        "[°C]",  # 2
        "[%]",  # 3
        "[kn]",  # 4
        "[°]",  # 5
        "[°]",  # 6
        "[°]",  # 7
        "[m]",  # 8
        "[m]",  # 9
        "[mPa]",  # 10
        "[µA]",  # 11
        "[°C]",  # 12
        "DU",  # 13, TO3
        "[°C]",  # 14, dewpt
        "[°C]",  # 15
        "[m/s]",  # 16
    ]

    # Map the neccessary profile variables to the columns of the file
    colmap4 = {
        "time": [0, 0],
        "press": [1, 1],
        "alt": [8, 8],
        "temp": [2, 2],
        "rhw": [3, 3],
        "dewpt": [14, 14],
        "wspd": [4, 4],
        "wdir": [5, 5],
    }

    cols5 = [
        "Time",  # 0
        "P",  # 1
        "T",  # 2
        "Hu",  # 3
        "Ws",  # 4
        "Wd",  # 5
        "Long.",  # 6
        "Lat.",  # 7
        "Alt",  # 8
        "Geopot",  # 9
        "Dewp.",  # 10
        "Virt. Temp",
        "Rs",
    ]

    units5 = [
        "[sec]",  # 0
        "[hPa]",  # 1
        "[°C]",  # 2
        "[%]",  # 3
        "[kn]",  # 4, wspd
        "[°]",  # 5, wdir
        "[°]",  # 6
        "[°]",  # 7
        "[m]",  # 8
        "[m]",  # 9
        "[°C]",  # 10
        "[°C]",
        "[m/s]",
    ]

    colmap5 = {
        "time": [0, 0],
        "press": [1, 1],
        "alt": [8, 8],
        "temp": [2, 2],
        "rhw": [3, 3],
        "dewpt": [10, 10],
        "wspd": [4, 4],
        "wdir": [5, 5],
    }

    if cols == cols1 and units == units1:
        return colmap1
    if cols == cols2 and units == units2:
        return colmap2
    if cols == cols3 and units == units3:
        return colmap3
    if cols == cols4 and units == units4:
        return colmap4
    if cols == cols5 and units == units5:
        return colmap5
    else:
        raise ValueError("Columns or units not correct")
