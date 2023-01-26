#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:25 2022

@author: prowe
"""

import pandas as pd  # Add via poetry
import numpy as np
import datetime as dt

# import pytz  # Add via poetry
# from matplotlib import pyplot as plt


def plot_broadband(fnames_swd, fnames_lwd, ax1, ax2):
    firsttime = True
    for filename in fnames_swd:
        # filename = dir_swd + fname

        # Read in the csv file for ifile
        df = pd.read_csv(filename, delimiter=",")

        # Get dates, skipping the 0th element, which is the lower part of the header
        dates_str = df["Date"] + " " + df["Time"]

        # Get rad and samples as float
        rad = np.array([float(x) for x in df["Radiation"]])
        nsamples = np.array([float(x) for x in df["#Samples"]])

        # Index to finite numbers
        ikeep = np.intersect1d(
            np.where(np.isfinite(rad))[0], np.where(np.isfinite(nsamples))[0]
        )

        # Get the date and add 3 hours because of UTC offset
        mydates = [
            dt.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in dates_str
        ]

        if firsttime:
            ax1.plot(mydates, rad, "b", label="measured")
            ax1.set_ylabel("Downward Shortwave (W/m$^{2}$)")
            firsttime = False
        else:
            ax1.plot(mydates, rad, "b")

    # Add on the longwave
    firsttime = True
    for filename in fnames_lwd:
        # filename = dir_lwd + fname

        # Read in the csv file for ifile
        df = pd.read_csv(filename, delimiter=",")

        # Get dates, skipping the 0th element, which is the lower part of the header
        dates_str = df["Date"] + " " + df["Time"]

        # Get rad and samples as float
        rad = np.array([float(x) for x in df["Radiation"]])
        nsamples = np.array([float(x) for x in df["#Samples"]])

        # Index to finite numbers
        ikeep = np.intersect1d(
            np.where(np.isfinite(rad))[0], np.where(np.isfinite(nsamples))[0]
        )

        # Get the date and add 3 hours because of UTC offset
        mydates = [
            dt.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in dates_str
        ]

        if firsttime:
            ax2.plot(mydates, rad, "b", label="measured")
            ax2.set_ylabel("Downward Longwave (W/m$^{2}$)")
            firsttime = False
        else:
            ax2.plot(mydates, rad, "b")
