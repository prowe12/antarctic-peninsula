#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:25 2022

@author: prowe
"""

import pandas as pd
import numpy as np
import datetime as dt
from netCDF4 import Dataset
from matplotlib import pyplot as plt


class Sonde:
    def __init__(self, fdir, fname, ffmt):
        self.dtime = dt.datetime.strptime(fname, ffmt)
        snd = Dataset(fdir + fname, "r")

        ikeep = np.where(snd["temp"][:].mask == False)[0]
        self.temp = snd["temp"][ikeep].data
        self.lev = snd["lev"][ikeep].data
        self.rh = snd["rh"][ikeep].data
        self.wdir = snd["wdir"][ikeep].data
        self.wspeed = snd["wspeed"][ikeep].data
        self.height = snd["height"][ikeep].data

    # def get_diffs(self, itm, pwrf):
    #     """
    #     Interpolate onto sonde levels and then get differences from sonde
    #     """
    #     slev = np.flipud(self.lev)
    #     plev = np.flipud(pwrf.level)
    #     temp = np.interp(slev, plev, np.flipud(pwrf.temp[:, itm].data))
    #     rh = np.interp(slev, plev, np.flipud(pwrf.rh[:, itm].data))
    #     wind = np.interp(slev, plev, np.flipud(pwrf.wind[:, itm].data))
    #     temp = np.flipud(temp)
    #     rh = np.flipud(rh)
    #     wind = np.flipud(wind)

    #     # Find the tropopause and plot the sonde and the interpolated PWRF
    #     itrp = self.find_tropoause()
    #     self.plot_snd_and_pwrf(temp, rh, wind, itrp)

    #     # Get the differences in the troposphere
    #     tempdiffs = get_diffs(temp, self.temp, itrp)
    #     rhdiffs = get_diffs(rh, self.rh, itrp)
    #     winddiffs = get_diffs(wind, self.wspeed, itrp)

    #     return tempdiffs, rhdiffs, winddiffs

    def get_diffs_pwrf(self, pwrf):
        """
        Interpolate onto sonde levels and then get differences from sonde
        """
        # sonde levels (pressures)
        sndlev = np.flipud(self.lev)

        # PWRF levels (pressures)
        plev = np.flipud(pwrf.level)

        sndtemp = np.interp(plev, sndlev, np.flipud(self.temp))
        sndrh = np.interp(plev, sndlev, np.flipud(self.rh))
        sndwind = np.interp(plev, sndlev, np.flipud(self.wind))
        sndht = np.interp(plev, sndlev, np.flipud(self.height)) / 1000
        sndtemp = np.flipud(sndtemp)
        sndrh = np.flipud(sndrh)
        sndwind = np.flipud(sndwind)
        sndht = np.fliupd(sndht)

        itrp = pwrf.find_tropopause()
        self.plot_pwrf_and_snd(pwrf, sndht, sndtemp, sndrh, sndwind)

        tdiff, rhdiff, wnddiff = self.get_trop_diffs(
            sndtemp, sndrh, sndwind, itrp
        )
        return tdiff, rhdiff, wnddiff

    def plot_pwrf_and_snd(self, pwrf, sndht, temp, rh, wind, itrp):
        """Create a figure as a QC"""
        plt.figure()
        plt.subplot(311)
        plt.plot(temp, sndht, label="PWRF")
        plt.plot(pwrf.temp, sndht, label="sonde")
        plt.plot(pwrf.temp[itrp], sndht[itrp], "*", label="tropopause")
        plt.legend()
        plt.ylabel("Height (km")
        plt.xlabel("Temperature (C)")
        plt.subplot(312)
        plt.plot(rh, sndht, label="PWRF")
        plt.plot(pwrf.rh, sndht, label="sonde")
        plt.plot(pwrf.rh[itrp], sndht[itrp], "*", label="tropopause")
        plt.xlabel("RH (%)")
        plt.subplot(313)
        plt.plot(wind, sndht, label="PWRF")
        plt.plot(pwrf.wspeed, sndht, ".-", label="sonde")
        plt.xlabel("Wind speed (m/s)")

    def get_diffs(self, pwrf):
        """
        Interpolate onto sonde levels and then get differences from sonde
        """
        # sonde levels (pressures)
        slev = np.flipud(self.lev)

        # PWRF levels (pressures)
        plev = np.flipud(pwrf.level)

        # Interpolate PWRF to sonde pressures, up to top of PWRF
        # itop = np.where(slev < min(plev))[0]

        temp = np.interp(slev, plev, np.flipud(pwrf.temp))
        rh = np.interp(slev, plev, np.flipud(pwrf.rh))
        wind = np.interp(slev, plev, np.flipud(pwrf.wind))
        temp = np.flipud(temp)
        rh = np.flipud(rh)
        wind = np.flipud(wind)

        itrp = self.find_tropopause()
        self.plot_snd_and_pwrf(temp, rh, wind, itrp)

        tdiff, rhdiff, wnddiff = self.get_trop_diffs(temp, rh, wind, itrp)
        return tdiff, rhdiff, wnddiff

    def find_tropopause(self):
        # Find the tropopause
        ht_km = self.height / 1000
        itrp = np.where(np.logical_and(ht_km > 9, ht_km < 15))[0]
        iitrop = np.where(np.diff(self.temp[itrp]) > 0)[0][0]
        return itrp[iitrop] + 1

    def plot_snd_and_pwrf(self, temp, rh, wind, itrp):
        ht_km = self.height / 1000

        # Create a figure as a QC
        plt.figure()
        plt.subplot(311)
        plt.plot(temp, ht_km, label="PWRF")
        plt.plot(self.temp, ht_km, label="sonde")
        plt.plot(self.temp[itrp], ht_km[itrp], "*", label="tropopause")
        plt.legend()
        plt.ylabel("Height (km")
        plt.xlabel("Temperature (C)")
        plt.subplot(312)
        plt.plot(rh, ht_km, label="PWRF")
        plt.plot(self.rh, ht_km, label="sonde")
        plt.plot(self.rh[itrp], ht_km[itrp], "*", label="tropopause")
        plt.xlabel("RH (%)")
        plt.subplot(313)
        plt.plot(wind, ht_km, label="PWRF")
        plt.plot(self.wspeed, ht_km, ".-", label="sonde")
        plt.xlabel("Wind speed (m/s)")

    def get_trop_diffs(self, temp, rh, wind, itrp):
        # Get the differences in the troposphere
        tempdiffs = get_diffs(temp, self.temp, itrp)
        rhdiffs = get_diffs(rh, self.rh, itrp)
        winddiffs = get_diffs(wind, self.wspeed, itrp)

        return tempdiffs, rhdiffs, winddiffs


class IgraSonde(Sonde):
    def __init__(self, fdir, fname, ffmt):
        snd = pd.read_csv(fdir + fname, comment="#", delim_whitespace=True)

        self.dtime = dt.datetime.strptime(fname, ffmt)
        self.temp = snd["TEMP"] / 10 - 273.15
        self.lev = snd["PRESS"] / 100
        self.rh = snd["CALCRH"] / 10
        self.wspeed = np.sqrt(
            (snd["UWND"] / 10) ** 2 + (snd["VWND"] / 10) ** 2
        )
        self.height = snd["CALCGPH"]

        i = np.where(snd["UWND"] < 0)[0]
        i = np.intersect1d(i, np.where(snd["VWND"] < 0)[0])
        self.wspeed[i] = np.nan


class DenialSonde(Sonde):
    def __init__(self, fdir, fname, ffmt):
        snd = pd.read_csv(
            fdir + fname,
            comment="#",
            delim_whitespace=True,
            skiprows=[0, 1, 2, 3, 4, 6],
        )

        self.dtime = dt.datetime.strptime(fname, ffmt)
        self.temp = snd["Temp"]
        self.lev = snd["P"]
        self.rh = snd["RH"]
        self.wspeed = snd["Speed"]  # speed in knots
        self.height = snd["HeightMSL"]

        # convert wind speed knots -> m/s
        self.wspeed = self.wspeed * 0.5144444444


class DenialSondePR(Sonde):
    def __init__(self, fdir, fname, ffmt):
        snd = pd.read_csv(
            fdir + fname,
            comment="#",
            delim_whitespace=True,
            skiprows=[0, 1, 2, 4],
        )

        self.dtime = dt.datetime.strptime(fname, ffmt)
        self.temp = snd["Temp"]
        self.lev = snd["P"]
        self.rh = snd["RH"]
        self.wspeed = snd["Speed"]  # speed in knots
        self.height = snd["HeightMSL"]

        # convert wind speed knots -> m/s
        self.wspeed = self.wspeed * 0.5144444444

        # Replace out-of-bound heights with nan
        self.height[self.height < 0] = np.nan


def get_diffs(val1, val2, itrp):
    max_diff = max(val1[:itrp] - val2[:itrp])
    min_diff = min(val1[:itrp] - val2[:itrp])
    mean_diff = np.mean(val1[:itrp] - val2[:itrp])
    rms_diff = np.sqrt(np.mean((val1[:itrp] - val2[:itrp]) ** 2))

    return [max_diff, min_diff, mean_diff, rms_diff]
