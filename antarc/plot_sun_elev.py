#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 15:37:17 2022

@author: prowe
"""

# Built-in modules
import datetime as dt
from calendar import monthrange
import numpy as np
import matplotlib.pyplot as plt

# My modules
from antarc.sun_position import sun_position


class PlotSunElevation:
    """
    Create plots of the sun elevation
    """

    def __init__(self, lat: float, lon: float, alt: float, dtime: dt.datetime):
        """
        Initialize
        @param lat  Latitude, decimal degrees (S is negative)
        @param lon  Longitude, decimal degrees (W is negative)
        @param alt  Altitude, m
        @param dtime  Datetime of interest
        """
        self.lat = lat
        self.lon = lon
        self.alt = alt
        self.dtime = dtime

    def get_solar_elevation(self):
        """
        Get the solar elevation for a particular date and location
        @param dtimes  The list of datetimes
        @return  The solar elevation
        """
        sun = sun_position(self.dtime, self.lat, self.lon, self.alt)
        return 90 - sun.zenith[0]

    def get_solar_elevations_for_dtimes(self, dtimes: list):
        """
        Get the solar elevation for a particular date and location
        @param dtimes  The list of datetimes
        @return  The solar elevation
        """
        elevation = np.zeros(len(dtimes))
        for i, datey in enumerate(dtimes):
            sun = sun_position(datey, self.lat, self.lon, self.alt)
            elevation[i] = 90 - sun.zenith[0]
        return elevation

    def get_dtimes_for_day(self) -> list:
        """
        Get the datetimes for the day
        @return list of datetimes
        """
        # Remove hours, minutes, seconds
        dtime_daystart = self.dtime.replace(hour=0, minute=0, second=0)

        # Date range: March 16 full day
        timeofday = np.arange(0.0, 24.0, 1 / 60)
        return [dtime_daystart + dt.timedelta(hours=x) for x in timeofday]

    def dtime_to_hour(self, dtimes: list) -> np.array:
        """
        Convert datetime to hour, assuming all same day
        @param dtime  List of datetimes
        @return hour
        """
        return np.array(
            [x.hour + x.minute / 60 + x.second / 60 / 60 for x in dtimes]
        )

    def diurnal(self):
        """
        Plot the solar elevation over the given day
        @return figno
        """
        dtimes = self.get_dtimes_for_day()
        elev = self.get_solar_elevations_for_dtimes(dtimes)
        hour = self.dtime_to_hour(dtimes)

        # Plot the solar elevation for a day
        figno = plt.figure(1)
        plt.subplot(position=[0.12, 0.12, 0.8, 0.8])
        plt.plot(hour, elev, label="solar elevation")
        plt.ylabel("solar elevation ($^o$)")
        plt.plot(hour, np.zeros(len(hour)), "k:")
        plt.xlabel("Hour on " + self.dtime.strftime("%B %-d, %Y"))
        plt.xlim([0, 24])
        plt.ylim(
            [np.floor(np.min(elev) / 5) * 5, np.ceil(np.max(elev) / 5) * 5]
        )

        return figno, hour, elev

    def annotate_diurnal(self, figno, hour, elev):
        """
        Annotate the diurnal figure with sunset, rise, and max
        @return solarnoon
        @return tset
        @return trise
        """

        plt.figure(figno)

        # Get max elevation: solar noon
        maxelev = np.max(elev)
        inoon = np.where(elev == maxelev)[0][0]
        tnoon = hour[inoon]  # .hour + dtimes[inoon].minute/60 \
        # + dtimes[inoon].second/60/60

        # Get sunset and sunrise
        if elev[0] > 0:
            # Starting at daytime, look for onset of nighttime,
            # then daytime again
            iset = np.where(elev < 0)[0][0]
            irise = np.where(elev[iset:] > 0)[0][0] + iset
        elif elev[0] < 0:
            # starting at night, look for onset of day, then night
            irise = np.where(elev > 0)[0][0]
            iset = np.where(elev[irise:] < 0)[0][0] + irise

        tset = hour[iset]  # dtimes[iset].hour + dtimes[iset].minute/60 \
        # + dtimes[iset].second/60/60

        trise = hour[irise]  # dtimes[irise].hour + dtimes[irise].minute/60 \
        # + dtimes[irise].second/60/60

        # Plot the solar elevation for a day
        figno = plt.figure(1)

        tnoon_str = str(round(tnoon, 2))
        tset_str = str(round(tset, 2))
        trise_str = str(round(trise, 2))

        plt.plot(
            tnoon,
            np.interp(tnoon, hour, elev),
            "r*",
            label="solar noon: " + tnoon_str + " UT",
        )
        plt.plot(
            tset,
            np.interp(tset, hour, elev),
            "kv",
            label="sunset:" + tset_str + " UT",
        )
        plt.plot(
            trise,
            np.interp(trise, hour, elev),
            "k^",
            label="sunrise:" + trise_str + " UT",
        )
        plt.legend()

        return tnoon, tset, trise

    def month(self):
        """
        Plot the solar elevation over the month including the given datetime
        @return figno
        """
        # number of days in this month
        ndays = monthrange(self.dtime.year, self.dtime.month)[1]

        # Plot solar elevation for month
        timeofday = np.arange(0.0, 24.0 * ndays, 1 / 15)
        dtime = [
            dt.datetime(2022, 3, 12) + dt.timedelta(hours=x) for x in timeofday
        ]

        # Get elevation
        elev = self.get_solar_elevations_for_dtimes(dtime)

        # Plot the solar elevation for the month
        figno = plt.figure(2)
        plt.subplot(position=[0.12, 0.12, 0.85, 0.85])
        plt.plot(dtime, elev, label="solar elevation")
        plt.ylabel("solar elevation ($^o$)")
        plt.xlim([dtime[0], dtime[-1]])
        plt.ylim(
            [np.floor(np.min(elev) / 5) * 5, np.ceil(np.max(elev) / 5) * 5]
        )

        return figno

    def daterange(self, date1, date2, ax):
        """
        Plot the solar elevation over the given period on the given ax
        @return dtimes
        @return elev
        """
        ndays = date2.day - date1.day
        timeofday = np.arange(0.0, 24.0 * ndays, 1 / 15)
        dtimes = [date1 + dt.timedelta(hours=x) for x in timeofday]

        elev = self.get_solar_elevations_for_dtimes(dtimes)

        # Plot the solar elevation for the month
        ax.plot(dtimes, elev, label="solar elevation")
        ax.plot(dtimes, np.zeros(len(dtimes)), "k:")
        ax.set_ylabel("solar elevation ($^o$)")
        ax.set_xlim([dtimes[0], dtimes[-1]])
        ax.set_ylim(
            [np.floor(np.min(elev) / 5) * 5, np.ceil(np.max(elev) / 5) * 5]
        )

        return dtimes, elev
