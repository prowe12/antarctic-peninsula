#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:25 2022

@author: prowe
"""

import datetime as dt
import numpy as np
from netCDF4 import Dataset


class LoadPwrf:
    def __init__(self, tempfile, rhfile, wndfile, starthour, dhour=1.0):
        self.level = []
        self.temp = []
        self.rh = []
        self.wind = []
        self.dtime = []

        self.level = self.get_level()
        self.temp = self.data_from_file(tempfile)
        self.rh = self.data_from_file(rhfile)
        self.wind = self.data_from_file(wndfile)
        self.dtime = self.get_dtime(tempfile, starthour, dhour)

    def get_level(self):
        level = np.array(
            [
                10,
                20,
                30,
                50,
                70,
                100,
                125,
                150,
                175,
                200,
                225,
                250,
                300,
                350,
                400,
                450,
                500,
                550,
                600,
                650,
                700,
                750,
                775,
                800,
                825,
                850,
                875,
                900,
                925,
                950,
                975,
                1000,
            ]
        )
        return np.flipud(level)

    def data_from_file(self, filename):
        with Dataset(filename, "r") as ncid:
            data = ncid["var"][:].data
        return data

    def get_dtime(self, filename, starttime, dtime):
        ntime = self.get_ntime(filename)
        passed_hours = [x * dtime for x in range(ntime)]
        return [starttime + dt.timedelta(hours=x) for x in passed_hours]

    def get_ntime(self, filename):
        with Dataset(filename, "r") as ncid:
            ntime = ncid.dimensions["time"].size
        return ntime

    def get_timeind(self, snd_dt):
        return np.argmin(np.abs([x - snd_dt for x in self.dtime]))


class Pwrf:
    def __init__(self, allPwrf, snd_dt):
        itime = allPwrf.get_timeind(snd_dt)
        self.level = allPwrf.level
        self.dtime = allPwrf.dtime[itime]
        self.temp = allPwrf.temp[:, itime]
        self.rh = allPwrf.rh[:, itime]
        self.wind = allPwrf.wind[:, itime]

        self.clean()

    def clean(self):
        igood = np.where(self.temp >= -273.15)[0]
        self.level = self.level[igood]
        self.temp = self.temp[igood]
        self.rh = self.rh[igood]
        self.wind = self.wind[igood]

        if np.any(self.rh < 0) or np.any(self.wind > 1e5):
            raise ValueError("Bad values in pwrf data")

    def trim(self, sndlev):
        """Trim the PWRF data outside the sonde bounds"""

        if sndlev[0] > self.level[0]:
            print("first pwrf point missing")

        ikeep = np.where(
            np.logical_and(
                self.level <= sndlev[0], self.level >= sndlev[len(sndlev) - 1]
            )
        )[0]

        # Get the surface point (highest pressure)
        temp0 = np.interp(sndlev[0], self.level[5::-1], self.temp[5::-1])
        rh0 = np.interp(sndlev[0], self.level[5::-1], self.rh[5::-1])
        wind0 = np.interp(sndlev[0], self.level[5::-1], self.wind[5::-1])

        # Get the last point, if it is within the range of pwrf
        lev1 = sndlev[len(sndlev) - 1]
        if lev1 > self.level[-1]:
            temp1 = np.interp(lev1, self.level[-1:-6:-1], self.temp[-1:-6:-1])
            rh1 = np.interp(lev1, self.level[-1:-6:-1], self.rh[-1:-6:-1])
            wind1 = np.interp(lev1, self.level[-1:-6:-1], self.wind[-1:-6:-1])
        else:
            temp1 = np.nan
            rh1 = np.nan
            wind1 = np.nan

        self.level = self.level[ikeep]
        self.temp = self.temp[ikeep]
        self.rh = self.rh[ikeep]
        self.wind = self.wind[ikeep]

        # Add on the surface point, if missing. It is missing if the
        # new PWRF surface pressure is lower than the sonde pressure
        if self.level[0] < sndlev[0]:
            self.level = np.hstack([sndlev[0], self.level])
            self.temp = np.hstack([temp0, self.temp])
            self.rh = np.hstack([rh0, self.rh])
            self.wind = np.hstack([wind0, self.wind])
        elif self.level[0] > sndlev[0]:
            raise ValueError("PWRF press level improperly trimmed?")

        # Add on the last point, if missing
        if lev1 < self.level[-1]:
            self.level = np.hstack([self.level, lev1])
            self.temp = np.hstack([self.temp, temp1])
            self.rh = np.hstack([self.rh, rh1])
            self.wind = np.hstack([self.wind, wind1])
