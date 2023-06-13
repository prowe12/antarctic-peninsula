#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:25 2022

@author: prowe
"""

import pandas as pd  # Add via poetry
import numpy as np
import datetime as dt


def plot_broadband(fnames_swd, fnames_lwd, ax1, ax2):
    firsttime = True
    for filename in fnames_swd:

        # Read in the csv file for ifile
        df = pd.read_csv(filename, delimiter=",")

        # Get dates, skipping the 0th element, which is the lower part of the header
        dates_str = df["Date"] + " " + df["Time"]

        # Get rad and samples as float
        rad = np.array([float(x) for x in df["Radiation"]])

        # Index to finite numbers
        # nsamples = np.array([float(x) for x in df["#Samples"]])
        # ikeep = np.intersect1d(
        #     np.where(np.isfinite(rad))[0], np.where(np.isfinite(nsamples))[0]
        # )

        # Get the date
        mydates = [
            dt.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in dates_str
        ]

        # Get the hourly averages
        radsum = 0
        naves = 0
        hour = mydates[0].hour
        dateave = []
        radave = []
        for i, mydate in enumerate(mydates):
            if mydate.hour == hour:
                radsum += rad[i]
                naves += 1
            else:
                dt0 = mydates[i - 1]
                thisdate = dt.datetime(
                    dt0.year, dt0.month, dt0.day, dt0.hour, 30, 30
                )
                dateave.append(thisdate)
                radave.append(radsum / naves)
                hour = mydate.hour
                radsum = rad[i]
                naves = 1

        if firsttime:
            ax1.plot(
                dateave, radave, "b", linewidth=8, label="measured, hourly ave"
            )
            ax1.plot(mydates, rad, "c", label="measured")
            ax1.set_ylabel("Downward Shortwave (W/m$^{2}$)")
            firsttime = False
        else:
            ax1.plot(dateave, radave, "b", linewidth=8)
            ax1.plot(mydates, rad, "c")

    # Add on the longwave
    firsttime = True
    for filename in fnames_lwd:

        # Read in the csv file for ifile
        df = pd.read_csv(filename, delimiter=",")

        # Get dates, skipping the 0th element, which is the lower part of the header
        dates_str = df["Date"] + " " + df["Time"]

        # Get rad and samples as float
        rad = np.array([float(x) for x in df["Radiation"]])

        # Index to finite numbers
        # nsamples = np.array([float(x) for x in df["#Samples"]])
        # ikeep = np.intersect1d(
        #     np.where(np.isfinite(rad))[0], np.where(np.isfinite(nsamples))[0]
        # )

        # Get the date
        mydates = [
            dt.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in dates_str
        ]

        # Get the hourly averages
        radsum = 0
        naves = 0
        hour = mydates[0].hour
        dateave = []
        radave = []
        for i, mydate in enumerate(mydates):
            if mydate.hour == hour:
                radsum += rad[i]
                naves += 1
            else:
                dt0 = mydates[i - 1]
                thisdate = dt.datetime(
                    dt0.year, dt0.month, dt0.day, dt0.hour, 30, 30
                )
                dateave.append(thisdate)
                radave.append(radsum / naves)
                hour = mydate.hour
                radsum = rad[i]
                naves = 1

        if firsttime:
            ax2.plot(dateave, radave, "b", linewidth=8)
            ax2.plot(mydates, rad, "c")
            ax2.set_ylabel("Downward Longwave (W/m$^{2}$)")
            firsttime = False

        else:
            ax2.plot(dateave, radave, "b", linewidth=8)
            ax2.plot(mydates, rad, "c")
