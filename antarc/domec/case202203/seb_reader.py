#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:47:08 2022

@author: prowe
"""

import datetime as dt
from netCDF4 import Dataset
import pandas as pd

from get_daily_ave_rad import get_daily_ave_rad


class SebReader:
    def __init__(self, file):
        dat = Dataset(file)

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
        day, radave = get_daily_ave_rad(self.dtime, getattr(self, field))
        return day, radave

    def getnet(self):
        return self.swnet + self.lwnet

    def getnet_dailyave(self):
        net = self.getnet()
        day, netave = get_daily_ave_rad(self.dtime, net)
        return day, netave


class Era5(SebReader):
    def __init__(self, era5_file):
        """
        @param era5_file  ERA5 filename"
        """

        era5 = Dataset(era5_file)
        # era5.variables.keys()
        # dict_keys(['LWD_tosta', 'time',
        #            'SWD_tosta', 'LWnet_tosta',
        #            'SWnet_tosta', 'lon_sta', 'lat_sta', 'elev_sta'])

        # print(f'ERA5 lat: {era5["lat_sta"][:].data}')
        # print(f'ERA5 lon: {era5["lon_sta"][:].data}')
        # print(f'ERA5 elev: {era5["elev_sta"][:].data}')

        # Create needed variables
        hour = range(0, 336, 1)
        hour = [x - 0.5 for x in hour]
        dtm = [dt.datetime(2022, 3, 12) + dt.timedelta(hours=x) for x in hour]

        self.dtime = dtm
        self.swd = era5["SWD_tosta"][:].data
        self.lwd = era5["LWD_tosta"][:].data
        self.swu = -era5["SWnet_tosta"][:].data + era5["SWD_tosta"][:].data
        self.lwu = -era5["LWnet_tosta"][:].data + era5["LWD_tosta"][:].data
        self.swnet = era5["SWnet_tosta"][:].data
        self.lwnet = era5["LWnet_tosta"][:].data

    def clean(self):
        self.swu[self.swu < 4e-05] = 0


class Pwrf(SebReader):
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


class Snowpack(SebReader):
    def __init__(self, spack_file):
        spack = pd.read_csv(spack_file, sep="\s+", comment="#")

        fmt = "%Y-%m-%dT%H:%M:%S"
        dtm = [dt.datetime.strptime(x, fmt) for x in spack["timestamp"]]

        # Get indices to March 12 - 26

        self.dtime = dtm
        self.swd = spack["incoming_short_wave_radiation"]
        self.lwd = spack["incoming_long_wave_radiation"]
        self.swu = spack["reflected_short_wave_radiation"]
        self.lwu = spack["outgoing_long_wave_radiation"]
        # self.swnet = spack["net_short_wave_radiation"]
        # self.lwnet = spack["net_long_wave_radiation"]
        self.swnet = self.swd - self.swu
        self.lwnet = self.lwd - self.lwu
