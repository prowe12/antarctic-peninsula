#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 09:21:30 2022

@author: prowe
"""

import numpy as np


def get_daily_ave_rad(dtimes, rad):
    if len(dtimes) != len(rad):
        raise ValueError("dates and radiation lists/arrays must have same len")

    radave = []
    radsum = 0
    npts = 0
    day = dtimes[0].day
    dayofmonth = [day]
    for i, dtm in enumerate(dtimes):
        if dtm.day != day:
            # wrap up the previous day
            radave.append(radsum / npts)

            # set up for the new day
            radsum = 0
            npts = 0
            day = dtm.day
            dayofmonth.append(day)

        # Add the radiance. If we have nan, find previous and next index that
        # are not nan and interpolate
        if np.isnan(rad[i]) and i == 0:
            raise ValueError("first value cannot be nan")
        elif np.isnan(rad[i]):
            dayfrac = [
                x.day
                + x.hour / 24
                + x.minute / 60 / 24
                + x.second / 60 / 24 / 24
                for x in dtimes[i - 1 : i + 2]
            ]
            print()
            print("old rad / new rad:")
            print(rad[i - 1], rad[i], rad[i + 1])
            rad[i] = np.interp(dayfrac[1], dayfrac, rad[i - 1 : i + 2])
            print(rad[i - 1], rad[i], rad[i + 1])

        radsum += rad[i]
        npts += 1

    # wrap up the previous day
    radave.append(radsum / npts)

    return dayofmonth, np.array(radave)
