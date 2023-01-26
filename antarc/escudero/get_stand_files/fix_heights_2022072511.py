#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 14:40:05 2018

@author: prowe

Purpose: The GRAW radiosonde file on 2022/07/25 11 has erroneous heights.
This code fixes them using the pressures at those levels to recompute
the heights.
"""

# Built-in modules
import datetime as dt
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


from hypsometric import hypsometric_for_z

debug = True


class SondeDataDenial:
    """Radiosonde in data denial format"""

    __slots__ = (
        "file",
        "date",
        "timestr",
        "height",
        "temp",
        "press",
        "rh",
        "dewpt",
        "wdir",
        "wspd",
    )

    def __init__(self, sonde_dir, sonde_file, sonde_date):
        """
        Initialize from a sounding file
        @param sonde_dir  Directory with sounding files
        @param sonde_file  Name or radiosonde file
        @param sonde_date  Date of radiosounding
        """

        file = sonde_dir + sonde_file
        sonde = pd.read_csv(file, header=[3, 4], delim_whitespace=True)
        self.timestr = sonde["TimeUTC"].to_numpy().T[0]
        self.press = sonde["P"].to_numpy().T[0]
        self.temp = sonde["Temp"].to_numpy().T[0]
        self.height = sonde["HeightMSL"].to_numpy().T[0]
        self.rh = sonde["RH"].to_numpy().T[0]
        self.dewpt = sonde["Dewp"].to_numpy().T[0]
        self.wdir = sonde["Dir"].to_numpy().T[0]
        self.wspd = sonde["Speed"].to_numpy().T[0]
        self.file = sonde_file

    def fix_bad_heights(self, imissing):
        """
        Estimate missing heights from the pressures etc.
        """
        if np.any(imissing):
            temp = self.temp + 273.15
            for imiss in imissing:
                zbot = self.height[imiss - 1]

                zfill = hypsometric_for_z(
                    self.press[imiss - 1 : imiss + 1],
                    temp[imiss - 1 : imiss + 1],
                    self.rh[imiss - 1 : imiss + 1],
                    zbot,
                )
                self.height[imiss] = zfill[1]
        if debug:
            plt.figure()
            plt.plot(self.temp, self.height)
            plt.xlabel("Temperature")
            plt.title(self.file)

    def save_in_datadenial_format(
        self, out_dir, datestr, timestr, loc, lat, lon, surfalt
    ):
        """
        WRITE THE OUTPUT FILE
        """
        colhdr = "   TimeUTC          P  HeightMSL       Temp         RH "
        colhdr += "      Dewp        Dir     Speed\n"
        colunits = "  hh:mm:ss        hPa          m       DegC          %"
        colunits += "       DegC    Degrees     Knots\n"
        with open(out_dir + self.file, "w+") as fid:

            # .. Write the header
            fid.write(loc + "\n")
            fid.write(
                "latitude {0:10.6f} longitude {1:10.6f} height {2:2.0f}m\n".format(
                    lat, lon, surfalt
                )
            )
            fid.write("Balloon release date and time                   ")
            fid.write(datestr[:4] + "-" + datestr[4:6] + "-" + datestr[6:8])
            fid.write("T" + timestr + "\n")
            fid.write(colhdr)
            fid.write(colunits)

            for i, height in enumerate(self.height):
                timestr = "  " + self.timestr[i]
                p_z_t_rh = "{0:11.1f} {1:10.0f} {2:10.1f} {3:10.0f}".format(
                    self.press[i], height, self.temp[i], self.rh[i]
                )
                dewpt_wdir_spd = "{0:11.1f} {1:10.0f} {2:9.1f}".format(
                    self.dewpt[i], self.wdir[i], self.wspd[i]
                )
                print(timestr + p_z_t_rh + dewpt_wdir_spd, file=fid)


indir = "/Users/prowe/Sync/measurements/Escudero/radiosondes/datadenial/orig/"
outdir = "/Users/prowe/Sync/measurements/Escudero/radiosondes/datadenial/"
sonde_file = "esc_sonde_dd_2022072511.txt"
sonde_date = dt.datetime(2022, 7, 25, 11, 0)

dtime = dt.datetime.strptime(sonde_file, "esc_sonde_dd_%Y%m%d%H.txt")
datestr = dtime.strftime("%Y%m%d")
timestr = "11:21:40"

loc = "esc_sonde_dd_"
lat = -62.201400
lon = -58.962200
srfht = 60

print(f"Working on: {sonde_file}")
sonde = SondeDataDenial(indir, sonde_file, sonde_date)

ibad = np.where(sonde.press >= 812.8)[0]
ibad = np.intersect1d(ibad, np.where(sonde.press <= 957.3)[0])
sonde.fix_bad_heights(ibad)
sonde.save_in_datadenial_format(outdir, datestr, timestr, loc, lat, lon, srfht)
