#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 12:29:24 2022

@author: prowe
"""

import datetime as dt
import pytz
import numpy as np

from antarc.plot_sun_elev import PlotSunElevation


class test_plot_sun_elev:
    def domec():
        lat = -(75 + 5 / 60 + 59 / 60 / 60)
        lon = 123.0 + 19 / 60 + 56 / 60 / 60
        alt = 3233
        dtime = dt.datetime(2022, 3, 16, tzinfo=pytz.utc)

        # Create a plot for the day
        plotSunElevation = PlotSunElevation(lat, lon, alt, dtime)

        solarnoon = plotSunElevation.diurnal(True)

        # Create plot for the month
        plotSunElevation.month()

    def escudero():
        # Setup
        lat = -(62 + 12 / 60 + 5 / 60 / 60)
        lon = 58.0 + 57 / 60 + 44 / 60 / 60
        alt = 60
        dtime = dt.datetime(2022, 2, 8, 8, 13, 0, tzinfo=pytz.utc)

        plotSunElevation = PlotSunElevation(lat, lon, alt, dtime)
        elev = plotSunElevation.get_solar_elevation()

        # Assert
        assert np.isclose(elev, 42.76)

        # Create a plot for the day
        # solarnoon = plotSunElevation.diurnal(True)

        # Create plot for the month
        # plotSunElevation.month()


# Tests that still need to be run
# Day starts at nighttime
# No sunrise
# No sunset
