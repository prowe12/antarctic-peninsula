#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:47:08 2022

@author: prowe
"""

from netCDF4 import Dataset
import datetime as dt

from get_daily_ave_rad import get_daily_ave_rad


class Era5:
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

    def getnet(self):
        return self.swnet + self.lwnet

    def getnet_dailyave(self):
        net = self.getnet()
        day, netave = get_daily_ave_rad(self.dtime, net)
        return day, netave
