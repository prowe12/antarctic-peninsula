#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:47:08 2022

@author: prowe
"""

import datetime as dt
from netCDF4 import Dataset
import numpy as np

from get_daily_ave_rad import get_daily_ave_rad


class Pwrf:
    def __init__(self, pwrf_file):
        dat = Dataset(pwrf_file)

        # pwrf.variables.keys():
        # dict_keys(['LWD_tosta', 'SWD_tosta', 'LWnet_tosta', 'SWnet_tosta',
        # 'lon_sta', 'lat_sta', 'elev_sta'])
        # print(f'PWRF lat: {dat["lat_sta"][:].data}')
        # print(f'PWRF lon: {dat["lon_sta"][:].data}')
        # print(f'PWRF elev: {dat["elev_sta"][:].data}')

        # Create needed variables
        pwrfhr = range(0, 109, 1)
        self.dtime = [
            dt.datetime(2022, 3, 16) + dt.timedelta(hours=x) for x in pwrfhr
        ]
        self.lwd = dat["LWD_tosta"][:].data
        self.swd = dat["SWD_tosta"][:].data
        self.lwnet = dat["LWnet_tosta"][:].data
        self.swnet = dat["SWnet_tosta"][:].data

        self.lwu = -self.lwnet + self.lwd
        self.swu = -self.swnet + self.swd

    def getaverage_rad(self, field):
        """
        Get the daily average
        @param dtm  datetimes
        @param rad  Radiation to get the average of, W/m2
        """
        # alldays = [
        #     x.day + x.hour / 24 + x.minute / 24 / 60 + x.second / 24 / 60 / 60
        #     for x in self.dtime
        # ]
        # alldays = np.array(alldays)

        day, radave = get_daily_ave_rad(self.dtime, getattr(self, field))
        return day, radave

    def getnet(self):
        return self.swnet + self.lwnet

    def getnet_dailyave(self):
        net = self.getnet()
        day, netave = get_daily_ave_rad(self.dtime, net)
        return day, netave
