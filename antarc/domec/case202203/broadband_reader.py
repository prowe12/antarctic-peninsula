#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 09:53:04 2022

@author: prowe
"""

import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from get_daily_ave_rad import get_daily_ave_rad


class Rad:
    """ """

    def __init__(self, file: str):
        """ """
        # File format
        datefmt = "%Y-%m-%d %H:%M:%S"
        # cols = [
        #     "date",
        #     "time",
        #     "SWD",
        #     "DIR",
        #     "DIF",
        #     "SWU",
        #     "LWD",
        #     "LWU",
        #     "GLO_spn1",
        #     "DIF_spn1",
        #     "DIR_spn1",
        # ]

        # Load the data
        data = pd.read_csv(file, sep="\s+", header=0)
        ndata = data.shape[0]

        # The datetimes
        dtime = [
            dt.datetime.strptime(
                data["date"][i] + " " + data["time"][i], datefmt
            )
            for i in range(ndata)
        ]

        # Quick QC - replace values < 0 with nan
        data.loc[data["SWD"] < -99, "SWD"] = np.nan
        data.loc[data["LWD"] < -99, "LWD"] = np.nan

        self.dtimearr = np.array(dtime)
        # TODO dayof month repeats!
        self.dayofmonth, self.swdave = self.getaverage(data["SWD"])
        self.dayofmonth, self.swuave = self.getaverage(data["SWU"])
        self.dayofmonth, self.lwdave = self.getaverage(data["LWD"])
        self.dayofmonth, self.lwuave = self.getaverage(data["LWU"])

        self.data = data
        self.dtime = dtime
        self.ndata = ndata

    def clean(self):
        i = np.where(self.data["SWD"] - self.data["SWU"] < -25)[0]
        self.data.loc[i, ("SWU")] = np.nan

    def getaverage(self, rad):
        """
        Get the daily average
        @param dtm  datetimes
        @param rad  Radiation to get the average of, W/m2
        """
        dtimes = self.dtimearr

        alldays = [
            x.day + x.hour / 24 + x.minute / 24 / 60 + x.second / 24 / 60 / 60
            for x in dtimes
        ]
        alldays = np.array(alldays)
        # radorig = copy.deepcopy(rad)

        # Interpolate over nans
        rad = rad.interpolate()

        # inan = np.where(np.isnan(radorig))[0]
        # plt.figure()
        # plt.plot(alldays, rad)
        # plt.plot(alldays[inan], rad[inan], '.')
        # plt.plot(alldays, radorig)

        dayofmonth, radave = get_daily_ave_rad(dtimes, rad)
        return dayofmonth, radave

    def plotqc(self):
        # Sample plot of all fields for QC
        cols = self.data.columns
        plt.figure(1)
        plt.title("All fields")
        plt.clf()
        for i in range(1, 10):
            plt.subplot(3, 3, i)
            plt.plot(self.dtime, self.data[cols[i + 1]][: self.ndata])
            plt.ylabel(cols[i + 1])

    def getnet(self):
        """
        Get the net LW and SW
        """
        self.data["LWnet"] = self.data["LWD"] - self.data["LWU"]
        self.data["SWnet"] = self.data["SWD"] - self.data["SWU"]
