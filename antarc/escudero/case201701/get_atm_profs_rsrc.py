#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 14:40:05 2018

@author: prowe
"""

# Built-in modules
import calendar
from datetime import datetime, timedelta
from bisect import bisect, bisect_right
import time
import datetime as dt
import pandas as pd
from scipy import interpolate
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, date2num, date2index
import numpy as np
import numpy.typing as npt

# My modules
from antarc.humidRS80 import humidRS80
from antarc.hypsometric import hypsometric
from antarc.getfilenames import get_filenames_dates
from antarc.hypsometric import hypsometric_for_z


# # # # # # #     Get file names     # # # # # # #
class FileInfo:
    """
    Get files within a directory with a given format,
    and get the dates using the format
    """

    __slots__ = ("dir", "files", "dates")

    def __init__(self, direc: str, fmt: str):
        files_n_dates = get_filenames_dates(direc, fmt)
        self.dir = direc
        self.files = [x[0] for x in files_n_dates]
        self.dates = [x[1] for x in files_n_dates]
        # files = os.listdir(directory)
        # files.sort()
        # files = list(filter(lambda x: len(x) == len(sample_filename), files))
        # dates = list(map(lambda f: datetime.strptime(f[i1:i2], fstr), files))


# # # # # # #     Radiosonde     # # # # # # #
class Sonde:
    """
    Class for radiosonde data from 'data denial' format
    """

    __slots__ = ("file", "date", "time", "z", "T", "P", "rh", "h2o")

    def __init__(
        self, sonde_dir: str, sonde_file: str, sonde_date: dt.datetime
    ):
        """Read in a radiosounding saved in data denial format
        @param sonde_dir  The directory with the soundings
        @param sonde_file  The desired file
        @param sonde_date  The corresponding datetime
        """

        file = sonde_dir + sonde_file
        sonde = pd.read_csv(
            file,
            header=[3, 4],
            # names=[
            #     "hh:mm:ss",
            #     "hPa",
            #     "m",
            #     "DegC",
            #     "%",
            #     "DegC",
            #     "Degrees",
            #     "Knots",
            # ],
            delim_whitespace=True,
        )  # ,
        # parse_dates=['hh:mm:ss']) #, date_parser=mydateparser)

        # sonde.columns:
        # MultiIndex([(  'TimeUTC', 'hh:mm:ss'),
        #             (        'P',      'hPa'),
        #             ('HeightMSL',        'm'),
        #             (     'Temp',     'DegC'),
        #             (       'RH',        '%'),
        #             (     'Dewp',     'DegC'),
        #             (      'Dir',  'Degrees'),
        #             (    'Speed',    'Knots')],
        #            )

        colnames = {
            "time": ("TimeUTC", "hh:mm:ss"),
            "press": ("P", "hPa"),
            "height": ("HeightMSL", "m"),
            "temp": ("Temp", "DegC"),
            "rh": ("RH", "%"),
            "dewpt": ("Dewp", "DegC"),
            "winddir": ("Dir", "Degrees"),
            "windspeed": ("Speed", "Knots"),
        }

        # Stupid way of getting the hours:
        hours = np.zeros(len(sonde[colnames["time"]]))
        for i in range(len(sonde[colnames["time"]])):
            hours[i] = (
                float(sonde[colnames["time"]][i][:2])
                + float(sonde[colnames["time"]][i][3:5]) / 60
                + float(sonde[colnames["time"]][i][6:9]) / 60 / 60
            )

        press = sonde[colnames["press"]].values
        temper = sonde[colnames["temp"]].values + 273.15
        alt = sonde[colnames["height"]].values
        rh = sonde[colnames["rh"]].values

        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(temper,alt)

        self.time = hours  # nci.variables['time'][:] + 0
        self.T = temper
        self.P = press
        self.rh = rh
        self.z = alt / 1000

        # Convert the rh to ppmv
        self.h2o = humidRS80(self.rh, self.T, "rhw", "ppmv", self.P)

        # Interpolate everything onto alt grid
        self.file = sonde_file
        yymmdd = datetime(sonde_date.year, sonde_date.month, sonde_date.day)
        self.date = yymmdd + timedelta(self.time[0] / 24)
        # timedelta()/60/60/24)

    def quality_control(self, srf: dict):
        """
        Implement quality control on sonde data
        surf: Dictionary of surface met data: 'temp', 'rhw', 'press'
        """

        # Check if any surface data is missing from the sounding
        if (
            (np.isnan(self.time[0]))
            or (np.isnan(self.z[0]))
            or (np.isnan(self.T[0]))
            or (np.isnan(self.P[0]))
            or (np.isnan(self.h2o[0]))
            or (np.isnan(self.rh[0]))
        ):
            print("add code here")

        datediff = (
            dt.datetime.strptime(srf["date"], "%Y-%m-%d %H:%M:%S") - self.date
        )
        # If the met data differs by more than 4 hours, do not bother using it
        # Otherwise require they be fairly close
        if abs(datediff.total_seconds()) / 3600 >= 4:
            print("  No met data within 4 hours for comparison")

        elif abs(datediff.total_seconds()) / 3600 < 4:  # 0.25:
            # Make sure surface data ~agrees with surface met data
            if srf["alt"] < self.z[0]:
                raise ValueError("met alt is lower than lowest sonde height")
            iht = np.where(self.z < 0.6)[0]  # lowest 600 m
            sndtemp = np.interp(srf["alt"], self.z[iht], self.T[iht])
            sndpress = np.interp(srf["alt"], self.z[iht], self.P[iht])
            sndrhw = np.interp(srf["alt"], self.z[iht], self.rh[iht])

            if abs(datediff.total_seconds()) / 3600 > 0.25:
                print("Big time diff from met; check this")

            # TODO
            # Use the hypsometric equation to interpolate to the sonde height

            if (
                np.abs(sndtemp - srf["temp"]) > 2
                or np.abs(sndrhw - srf["rhw"]) > 25
            ):
                print(
                    "  Surface met temp and rh data differ by "
                    + f"{round(np.abs(sndtemp - srf['temp']),1)} and "
                    + f"{round(np.abs(sndrhw - srf['rhw']),1)}"
                )

            if np.abs(sndpress - srf["press"]) > 5:
                print(
                    "  Surface pressures differ by "
                    + f"{np.round(np.abs(sndpress - srf['press']),2)} mb"
                )

        # Pressures must decrease with altitude/time
        if any(np.diff(self.P) >= 0):
            raise ValueError("One or more pressures increase with altitude")

        # Do not allow any NaNs in pressure
        if any(np.isnan(self.P)):
            # TODO: fix this
            raise ValueError("One or more pressures is nan")

        # Heights must increase with time
        if any(np.diff(self.z) <= 0):
            ibad = np.where(np.diff(self.z) <= 0)[0]
            if len(ibad) == 1:  # or np.all(np.diff(ibad)>1):
                print("  Warning: one altitude decreases with time, fixing")
                self.z[ibad[0] + 1] = np.mean(
                    [self.z[ibad[0]], self.z[ibad[0] + 2]]
                )
            else:
                print("  Multiple altitudes decrease with time; correcting")
                # Set the bad values to nans and fix below
                self.z[ibad + 1] = np.nan

        if any(np.isnan(self.time)):
            print("  One or more time is nan; ignoring")

        if any(np.isnan(self.T)):
            if np.isnan(self.T[0]):
                raise ValueError("Surface temperature is NaN")
            if np.isnan(self.T[-1]):
                raise ValueError("Final temperature is NaN")

            print("  One or more temperature is nan; interpolating")
            ibad = np.where(np.isnan(self.T))[0]
            ilo = ibad[0] - 1
            ihi = ibad[-1] + 1
            self.T[ibad] = np.interp(
                self.z[ibad],
                [self.z[ilo], self.z[ihi]],
                [self.T[ilo], self.T[ihi]],
            )

        if any(np.isnan(self.h2o)) or any(np.isnan(self.rh)):
            print("  H2O has NaNs - need to correct with ERA")

        # Fix any missing heights - MUST USE hypsometric_for_z
        # BECAUSE SOME VALUES DECREASE !!!!
        if np.any(self.z == -999):
            raise ValueError("This should not happen; see code")
        if any(np.isnan(self.z)):
            if np.isnan(self.z[0]):
                raise ValueError("Need surface alt")
            imiss = np.where(np.isnan(self.z))[0]
            if imiss[0] == 0:
                imiss = np.hstack([imiss, imiss[-1] + 1])
            elif imiss[-1] == len(self.z) - 1:
                imiss = np.hstack([imiss[0] - 1, imiss])
            else:
                imiss = np.hstack([imiss[0] - 1, imiss, imiss[-1] + 1])
            alt0 = self.z[imiss[0]] * 1000
            zfill = (
                hypsometric_for_z(
                    self.P[imiss], self.T[imiss], self.rh[imiss], alt0
                )
                / 1000
            )
            self.z[imiss[:-1]] = zfill[:-1]

            # Double check that final value makes sense
            if len(self.z) < imiss[-1] + 11:
                mean_dz = np.diff(self.z[imiss[-1] + 1 : imiss[-1] + 2])
                print("  Only using 1 value to get mean dz")
            else:
                mean_dz = np.mean(
                    np.diff(self.z[imiss[-1] + 1 : imiss[-1] + 10])
                )
            dztop = self.z[imiss[-1]] - zfill[-1]
            if dztop > 0.600 or dztop < -0.600:
                raise ValueError("Bad interpolation with hypsometric eqn")
            elif dztop < 0 or (self.z[imiss[-1]] - zfill[-1] > 2 * mean_dz):
                # Adjust altitudes so max difference is ~mean_dz
                adj = np.linspace(0, dztop, len(zfill))
                zfill += adj
            self.z[imiss] = zfill


# # # # # # #     Station flask CO2     # # # # # # #


class CO2stationData:
    """Load CO2 data measured at the surface at stations"""

    __slots__ = ["date", "co2", "lat", "lon", "alt"]

    def __init__(self, co2_file):
        """
        Get the data from the co2 file
        @param co2_file  Name of the file containing co2 data with fields below
        """
        # site_code year month day hour minute second datetime time_decimal
        # air_sample_container_id value value_unc latitude longitude altitude
        # elevation intake_height method event_number instrument
        # analysis_datetime qcflag

        stn = pd.read_csv(
            co2_file,
            delim_whitespace=True,
            comment="#",
            usecols=[
                "year",
                "month",
                "day",
                "hour",
                "minute",
                "second",
                "value",
                "latitude",
                "longitude",
                "altitude",
            ],
        )
        list(stn.columns.values)

        # Get datetime
        self.date = list(
            map(
                lambda year, month, day, hour, minute, second: dt.datetime(
                    year, month, day, hour, minute, second
                ),
                stn["year"],
                stn["month"],
                stn["day"],
                stn["hour"],
                stn["minute"],
                stn["second"],
            )
        )
        self.co2 = stn["value"]
        self.lat = stn["latitude"]
        self.lon = stn["longitude"]
        self.alt = stn["altitude"]


class CarbonTracker:
    """
    Load global carbon tracker data for a given date, latitude and longitude
    """

    __slots__ = ("co2", "T", "P", "zbnd", "z", "rh")

    def __init__(self, co2Files, date, lat, lon):
        """
        Load co2 data nc file
        """

        # Find closest date
        iaft = bisect(co2Files.dates, date)
        ibef = iaft - 1

        # Sometimes the "before" date is actually slightly after
        # the first date in the file. In this case, decrease ibef and iaft
        with Dataset(co2Files.dir + co2Files.files[ibef], "r") as nci:
            dates = num2date(
                nci.variables["time"][:], nci.variables["time"].units
            )
        if min(dates) > date:
            ibef -= 1
            iaft -= 1

        # Open the netcdf co2 file and get the data
        with Dataset(co2Files.dir + co2Files.files[ibef], "r") as nci:
            """
            # nci['co2'].shape
            # Out[32]: (8,    25,   90, 120)
            #         time, level, lat, lon
            """
            nlev = nci["level"].shape[0]
            nbound = nci["boundary"][:].shape[0]
            lats = nci["latitude"][:]
            lons = nci["longitude"][:]
            dates = num2date(
                nci.variables["time"][:], nci.variables["time"].units
            )

            """
            # Verify it worked
            print('units = %s, values = %s' % (nci.variables['time'].units, 
                                               nci.variables['time'][:]))
            print(nci['time_components'][:])
            print(dates)
            print([date.strftime('%Y-%m-%d %H:%M:%S') for date in dates])
            """

            if max(dates) < date:
                # Glom on the last with the first
                co2mat = np.zeros((2, nlev, 90, 120))
                tmat = np.zeros((2, nlev, 90, 120))
                pmat = np.zeros((2, nbound, 90, 120))
                zmat = np.zeros((2, nbound, 90, 120))
                co2mat[0, :, :, :] = nci["co2"][-1, :, :, :]
                tmat[0, :, :, :] = nci["temperature"][-1, :, :, :]
                pmat[0, :, :, :] = nci["pressure"][-1, :, :, :]
                zmat[0, :, :, :] = nci["gph"][-1, :, :, :]

                with Dataset(
                    co2Files.dir + co2Files.files[iaft], "r"
                ) as nc_aft:
                    # nc_aft = Dataset(co2_dir + co2Files.files[iaft], 'r')
                    dates_aft = num2date(
                        nc_aft.variables["time"][:],
                        nc_aft.variables["time"].units,
                    )
                    dates = [dates[-1], dates_aft[0]]

                    # Glom on the last with the first
                    co2mat[1, :, :, :] = nc_aft["co2"][0, :, :, :]
                    tmat[1, :, :, :] = nc_aft["temperature"][0, :, :, :]
                    pmat[1, :, :, :] = nc_aft["pressure"][0, :, :, :]
                    zmat[1, :, :, :] = nc_aft["gph"][0, :, :, :]

                    co2 = map_interp(co2mat, date, lat, lon, dates, lats, lons)
                    temp = map_interp(tmat, date, lat, lon, dates, lats, lons)
                    press = map_interp(pmat, date, lat, lon, dates, lats, lons)
                    alt = map_interp(zmat, date, lat, lon, dates, lats, lons)

                    del dates, dates_aft, co2mat, tmat, pmat, zmat

            else:
                co2mat = nci["co2"][:]
                tmat = nci["temperature"][:]
                pmat = nci["pressure"][:]
                zmat = nci["gph"][:]

                co2 = map_interp(co2mat, date, lat, lon, dates, lats, lons)
                temp = map_interp(tmat, date, lat, lon, dates, lats, lons)
                press = map_interp(pmat, date, lat, lon, dates, lats, lons)
                alt = map_interp(zmat, date, lat, lon, dates, lats, lons)

                del dates, co2mat, tmat, pmat, zmat

            del nlev, nbound, lats, lons

            self.co2 = co2
            self.T = temp
            self.P = press / 100
            self.zbnd = alt / 1000

        self.z = self.zbnd[:-1] + np.diff(self.zbnd) / 2


# # # # # # #     Profile     # # # # # # #
class Prof:
    """
    Class for the atmosperic profile
    """

    __slots__ = (
        "file",
        "date",
        "time",
        "z",
        "T",
        "P",
        "rh",
        "h2o",
        "co2",
        "o3",
        "f11",
        "f12",
        "f113",
        "hno3",
        "units",
        "model_extra",
    )

    def __init__(self, sonde, layerbnds):
        nlyr = len(layerbnds)
        self.z = np.array(layerbnds)
        self.T = np.nan * np.ones(nlyr)
        self.P = np.nan * np.ones(nlyr)
        self.rh = np.nan * np.ones(nlyr)
        self.h2o = np.nan * np.ones(nlyr)
        self.co2 = np.nan * np.ones(nlyr)
        self.o3 = np.nan * np.ones(nlyr)

        # Guess values (based on fit to a single measurement,
        #    needs improvement)
        self.f113 = 0.00008 * np.ones(nlyr)  # 800-828, 870-930,
        self.f11 = 0.00026 * np.ones(nlyr)  # 830-860, 1060-1107
        self.hno3 = 0.002 * np.ones(nlyr)  # 860-920
        self.f12 = 0.00053 * np.ones(nlyr)  # 867-937, 1080-1177

        # Units for ozone are jchar = C for g/kg
        self.units = dict(
            {
                ("z", "km"),
                ("T", "K"),
                ("P", "mb"),
                ("h2o", "ppmv"),
                ("co2", "ppmv"),
                ("o3", "gm_kg"),
                ("hno3", "ppmv"),
                ("f11", "ppmv"),
                ("f12", "ppmv"),
                ("f113", "ppmv"),
            }
        )

        # Interpolate everything onto z grid
        self.file = sonde.file
        self.date = sonde.date
        self.time = np.interp(layerbnds, sonde.z, sonde.time)

        inds = np.intersect1d(
            np.where(self.z >= sonde.z[0])[0],
            np.where(self.z <= sonde.z[len(sonde.z) - 1])[0],
        )
        self.T[inds] = np.interp(self.z[inds], sonde.z, sonde.T)
        self.P[inds] = np.interp(self.z[inds], sonde.z, sonde.P)
        self.h2o[inds] = np.interp(self.z[inds], sonde.z, sonde.h2o)
        self.rh[inds] = np.interp(self.z[inds], sonde.z, sonde.rh)

    def set(self, attr_name, new):
        """
        Set an attribute to a new value, interpolated onto the altitudes
        and warn user if any extrapolation occurs
        @param attr_name  the name of the attribute
        @param new  A class with fields that include the attribute and z
        """
        new_val = getattr(new, attr_name)
        new_val = np.interp(self.z, new.z, new_val)

        setattr(self, attr_name, new_val)

        if new.z[0] > self.z[0]:
            print(
                "  Extrapolating ",
                attr_name,
                "for heights: ",
                self.z[self.z < new.z[0]],
            )
        if new.z[-1] < self.z[-1]:
            print(
                "  Extrapolating ",
                attr_name,
                "for prof.",
                attr_name,
                "for heights: ",
                self.z[self.z > new.z[-1]],
            )

        """
        igood = np.where(np.isnan(prof.co2)==False)[0]
        should_be_done = False
        for i in range(len(prof.co2)):
            if np.isnan(prof.co2[i]):
                if should_be_done:
                    raise NameError('Should be done with NaNs in CO2')
                prof.co2[i] = prof.co2[igood[0]]
            else:
                should_be_done = True
        """

    def set_n_scale(self, attr_name, new, surf_value):
        """
        Set a new value for an attribute from a profile and then scale the
        profile to the surface value
        @param attr_name  the name of the attribute
        @param new  A class with fields that include the attribute and z
        @param surf_value  The surface value to scale to
        """

        new_val = getattr(new, attr_name)
        new_val = np.interp(self.z, new.z, new_val)
        new_val = surf_value * new_val / new_val[0]

        setattr(self, attr_name, new_val)

    def set_upper(self, attr_name, new):
        """
        Set new values for an attribute from a profile where the old values
        for that attribute are NaN, which should be in the upper atmosphere
        using linear interpolation
        @param attr_name  the name of the attribute
        @param new  A class with fields that include the attribute and z
        """
        # Get the old and new values to work with
        old_val = getattr(self, attr_name)
        new_val = getattr(new, attr_name)
        new_val = np.interp(self.z, new.z, new_val)

        inds = np.where(np.isnan(old_val))
        if len(inds) > 0:
            inds2 = np.intersect1d(
                np.where(self.z >= new.z[0])[0],
                np.where(self.z <= new.z[-1])[0],
            )
            if len(inds2) > 0:
                inds = np.intersect1d(inds, inds2)

            old_val[inds] = new_val[inds]
            setattr(self, attr_name, old_val)

    def set_upper_spline(self, attr_name, new):
        """
        Set new values for an attribute from a profile where the old values
        for that attribute are NaN, which should be in the upper atmosphere
        using spline interpolation
        @param attr_name  the name of the attribute
        @param new  A class with fields that include the attribute and z
        """
        # Get the old and new values to work with
        old_val = getattr(self, attr_name)
        new_val = getattr(new, attr_name)

        # Do the interpolation. Do not use spline because it can
        # distort the pressures
        # tck = interpolate.splrep(new.z, new_val)
        # new_val = interpolate.splev(self.z, tck)
        new_val = np.interp(self.z, new.z, new_val)

        inds = np.intersect1d(
            np.where(self.z >= new.z[0])[0], np.where(self.z <= new.z[-1])[0]
        )
        inds = np.intersect1d(inds, np.where(np.isnan(old_val)))
        old_val[inds] = new_val[inds]
        # setattr(self, old_attr_name, old_val)

        # Spline does not work very well for the last few points
        inds = np.where(self.z > 47)[0]  # new.z[-2]) #
        if len(inds) == 1:
            raise ValueError("len(inds) must be > 1")
        alti = self.z[inds]
        # print('alti[0] =', alti[0], ' in get_atm_profs_KGI_rsrc.py, line 197')
        # print('alti[-1] =', alti[-1], ' in get_atm_profs_KGI_rsrc.py')
        tempi = self.T[inds]
        rhi = 0 * alti
        ptop = hypsometric(alti * 1000, tempi, rhi, self.P[inds[0]])

        old_val[inds] = ptop
        setattr(self, attr_name, old_val)

    def error_check(self, era):
        """
        Check for errors, including NaNs, pressures that do not decrease, and
        heights that do not increase, and make sure the temperature
        differential across layers is not too large
        @param era  ERA data, used to fill in upper-level H2O
        """
        if (
            np.any(np.isnan(self.z))
            or np.any(np.isnan(self.P))
            or np.any(np.isnan(self.co2))
            or np.any(np.isnan(self.h2o))
            or np.any(np.isnan(self.o3))
            or np.any(np.isnan(self.f11))
            or np.any(np.isnan(self.f12))
            or np.any(np.isnan(self.f113))
            or np.any(np.isnan(self.T))
        ):

            if np.any(np.isnan(self.z)):
                raise NameError("One or more z is NaN!")
            if np.any(np.isnan(self.P)):
                raise NameError("One or more P is NaN!")
            if np.any(np.isnan(self.co2)):
                raise NameError("One or more co2 is NaN!")
            if np.any(np.isnan(self.h2o)):
                print("  Warning: one or more h2o is NaN, filling with era")

                # Convert the ERA rh to ppmv
                h2o = humidRS80(era.rh, era.T, "rhw", "ppmv", era.P)
                inan = np.where(np.isnan(self.h2o))[0]

                h2o_z = np.interp(self.z[inan], era.z, h2o)
                self.h2o[inan] = h2o_z

            if np.any(np.isnan(self.o3)):
                raise NameError("One or more o3 is NaN!")
            if np.any(np.isnan(self.f11)):
                raise NameError("One or more f11 is NaN!")
            if np.any(np.isnan(self.f12)):
                raise NameError("One or more f12 is NaN!")
            if np.any(np.isnan(self.f113)):
                raise NameError("One or more f113 is NaN!")
            if np.any(np.isnan(self.T)):
                raise NameError("One or more T is NaN!")

        # Check for heights (z) that do not monotonically increase or
        # pressures (P) that do not monotonically decrease
        # to within 3 significant figures
        zmr = np.round(self.z, 3)
        pmr = np.round(self.P, 3)
        if np.any(np.diff(zmr) <= 0):
            raise ValueError("zs do not increase, or not enough.")

        if (np.any(np.diff(pmr) >= 0)) or (np.any(pmr == 0)):
            raise ValueError(
                "Problem with pressure: not decreasing enough "
                "or 0 detected."
            )

        # Check for changes in temperature greater than 8 in lower 20 km
        if np.any(np.abs(np.diff(self.T[self.z <= 20])) > 8):
            maxdtemp = round(np.max(np.diff(self.T[self.z <= 20])), 1)
            print(
                "  Warning: Temperature differences in lower 20 km are up "
                + "to "
                + str(maxdtemp)
                + " K but should be < 8 K "
                + "(getAtmProfs_KGI_rsrc line 465."
            )

        # Ignore this for now, although it will be a problem if we
        # eventually run DISORT
        # Check for changes in temperature greater than 10 K anywhere
        # if np.any(np.abs(np.diff(self.T)) > 10):
        #     maxdtemp = round(np.max(np.abs(np.diff(self.T))), 1)
        #     print(self.T)

        #     plt.figure()
        #     plt.plot(self.T, self.z, ".-")
        #     print(
        #         "Warning: Temperature differences above 20 km are up to "
        #         + str(maxdtemp)
        #         + " K but should be < 12 K."
        #     )

    def set_surf_temp(self, zmet, tmet):
        """
        Set the surface temperature using surface meteorological data
        @param zmet  Height of surface data
        @param tmet  Temperature of surface data
        """

        # Interpolate using the first 600 m
        iht = np.where(self.z <= 0.6)[0]

        zsurf = self.z[0]
        if zmet < zsurf:
            raise NameError(
                "Met height is below sounding surface height. "
                + "Add code to interpolate to sounding surface."
            )
        elif (zmet == zsurf) or (zmet > self.z[iht[-1]]):
            raise ValueError("Modify code for this possibility")

        # Interpolate the sonde temperature to the height of the met temp
        tsnd_surf = np.interp(zmet, self.z[iht], self.T[iht])

        # If they are within 0.5 degree, all is well, otherwise need to fix
        if abs(tsnd_surf - tmet) > 1:
            # Get new temperatures for iht indices
            ifirst = np.where(self.z >= zmet)[0][0]
            tnew = np.interp(
                self.z[iht],
                [zmet, self.z[iht[-1]]],
                [tmet, self.T[iht[-1]]],
            )

            # Fix the surface-to-met height points by assuming const dT
            dtemp_dz = (tnew[ifirst + 1] - tnew[ifirst]) / (
                self.z[ifirst + 1] - self.z[ifirst]
            )
            ilow = np.arange(ifirst - 1, -1, -1)
            for i in ilow:
                dtemp = -dtemp_dz * (self.z[i + 1] - self.z[i])
                tnew[i] = tnew[i + 1] + dtemp

            # Do not do the following because of possibility of bad data
            # # Merge gradually to old profile
            # falloff = np.linspace(1, 0, len(iht))
            # self.T[iht] = falloff * tnew + (1 - falloff) * self.T[iht]

    """
    def set_upper_P_carbonTracker(self, new, old_attr_name, new_attr_name):
        # Get the old and new values to work with
        old_val = getattr(self, old_attr_name)
        new_val = getattr(new, new_attr_name)
        
        # Do the spline interpolation 
        tck = interpolate.splrep(new.zbnd, new_val)
        new_val = interpolate.splev(self.z, tck)
        
        inds = np.intersect1d(np.where(self.z >= new.zbnd[0])[0],
                              np.where(self.z <= new.zbnd[-1])[0])
        inds = np.intersect1d(inds, np.where(np.isnan(old_val)))
        old_val[inds] = new_val[inds]
        #setattr(self, old_attr_name, old_val)
    """

    def write_to_netcdf_file(self, outputfile):
        """
        Write the results to a netcdf file
        @param outputfile  The full output file name
        """
        with Dataset(outputfile, "w", format="NETCDF4_CLASSIC") as nci:
            # Create Dimensions
            nci.createDimension("level", len(self.z))
            nci.createDimension("time", 1)
            nci.createDimension("const", 1)
            # lat = nci.createDimension('lat', 1)
            # lon = nci.createDimension('lon', 1)

            # Create variables
            nc_time = nci.createVariable("time", np.float64, ("time",))
            nc_z = nci.createVariable("z", np.float32, ("level"))
            nc_t = nci.createVariable("T", np.float32, ("level"))
            nc_p = nci.createVariable("P", np.float32, ("level"))
            nc_rh = nci.createVariable("rh", np.float32, ("level"))
            nc_h2o = nci.createVariable("h2o", np.float32, ("level"))
            nc_co2 = nci.createVariable("co2", np.float32, ("level"))
            nc_o3 = nci.createVariable("o3", np.float32, ("level"))
            nc_f11 = nci.createVariable("f11", np.float32, ("level"))
            nc_f12 = nci.createVariable("f12", np.float32, ("level"))
            nc_f113 = nci.createVariable("f113", np.float32, ("level"))
            nc_hno3 = nci.createVariable("hno3", np.float32, ("level"))
            nc_model_extra = nci.createVariable(
                "model_extra", np.int8, ("const")
            )

            # Global attributes
            nci.filename = self.file
            nci.history = "Created " + time.ctime(time.time())

            # Variable attributes
            nc_time.units = "hours since 0001-01-01 00:00:00"
            nc_time.calendar = "gregorian"
            nc_z.units = self.units["z"]
            nc_t.units = self.units["T"]
            nc_p.units = self.units["P"]
            nc_rh.units = "%"
            nc_h2o.units = self.units["h2o"]
            nc_co2.units = self.units["co2"]
            nc_o3.units = self.units["o3"]
            nc_hno3.units = self.units["hno3"]
            nc_f11.units = self.units["f11"]
            nc_f12.units = self.units["f12"]
            nc_f113.units = self.units["f113"]

            # Assign values
            nc_time[:] = date2num(
                self.date, units=nc_time.units, calendar=nc_time.calendar
            )
            nc_z[:] = self.z
            nc_t[:] = self.T
            nc_p[:] = self.P
            nc_rh[:] = self.rh
            nc_h2o[:] = self.h2o
            nc_co2[:] = self.co2
            nc_o3[:] = self.o3
            nc_f11[:] = self.f11
            nc_f12[:] = self.f12
            nc_f113[:] = self.f113
            nc_hno3[:] = self.hno3
            nc_model_extra[:] = self.model_extra


def get_co2(sonde_date, sonde_lat, co2_drake, co2_palmer):
    """
    Get surface CO2 data from Drake Passage or Palmer station, if available.
    to correspond with sounding. If not available, use best estimate from data
    @param sonde_date  Date of sounding
    @param sonde_lat  Latitude of sounding
    @param co2_drake  Surface CO2 measurements with time in Drake Passage
    @param co2_palmer  Surface CO2 measurements with time at Palmer
    """
    # # # #
    # Uncomment and put text below this line in class
    # co2 = CO2(co2_file_palmer, co2_file_drake, sonde.date, lat, lon)
    #

    def to_tstamp(dtime):
        """
        Convert date to timestamp
        @param dtime  Date
        """
        return calendar.timegm(dtime.timetuple())

    def interp_co2_to_date(sonde_date, sta_date, sta_co2, sta_lat):
        """
        Interpolate co2 to a date
        @param sonde_date  Date of sounding
        @param co2date  Date corresponding to surface co2 measurement
        @param co2  CO2 measurement
        @param co2lat  Latitude of CO2 measurement
        """
        iaft = bisect(sta_date, sonde_date)
        ibef = iaft - 1

        # Sometimes the "before" date is actually slightly after
        # the first date in the file. In this case, decrease ibef and iaft
        if min(sta_date) > sonde_date:
            ibef -= 1
            iaft -= 1

        # Get the Drake Passage CO2 value
        tstp = np.array([to_tstamp(sta_date[ibef]), to_tstamp(sta_date[iaft])])
        co2 = np.interp(to_tstamp(sonde_date), tstp, sta_co2[ibef : iaft + 1])
        lat = np.interp(to_tstamp(sonde_date), tstp, sta_lat[ibef : iaft + 1])

        return co2, lat

    # QC
    if sonde_date < co2_drake.date[0] and sonde_date < co2_palmer.date[0]:
        msg1 = "sonde_date ({sonde_date}) before first date in co2 file "
        msg2 = "for Drake Passage ({co2_drake.date[0]})"
        msg3 = " and Palmer Station ({co2_palmer.date[0]})"
        raise ValueError(msg1 + msg2 + msg3)

    # Find closest date (Drake)
    use_drake = True
    if sonde_date < co2_drake.date[0]:
        use_drake = False
    elif sonde_date > co2_drake.date[-1]:
        # Cheat for now - use the last date
        print(f"Warning: No CO2 info from Drake Passage for {sonde_date}")
        print(f"Using final date: {co2_drake.date[-1]}")
        co2_drake0 = co2_drake.co2.values[-1]
        lat_drake0 = co2_drake.lat.values[-1]
    else:
        co2_drake0, lat_drake0 = interp_co2_to_date(
            sonde_date, co2_drake.date, co2_drake.co2, co2_drake.lat
        )

    use_palmer = True
    if sonde_date < co2_palmer.date[0]:
        use_palmer = False
    elif sonde_date > co2_palmer.date[-1]:
        # Cheat for now - use the last date
        print(f"Warning: No CO2 info from Palmer station for {sonde_date}")
        print(f"Using final date: {co2_palmer.date[-1]}")
        co2_palmer0 = co2_palmer.co2.values[-1]
        lat_palmer0 = co2_palmer.lat.values[-1]
    else:
        co2_palmer0, lat_palmer0 = interp_co2_to_date(
            sonde_date, co2_palmer.date, co2_palmer.co2, co2_palmer.lat
        )

    # Interpolate to desired latitude
    if use_palmer and use_drake:
        fco2 = interpolate.interp1d(
            [lat_drake0, lat_palmer0], [co2_drake0, co2_palmer0]
        )
        return fco2(sonde_lat)

    if use_palmer:
        return co2_palmer0

    if use_drake:
        return co2_drake

    raise ValueError("No co2 data to use")


class Era:
    """
    Get profile information from ERA data
    """

    def __init__(self, fname, date, lat, lon):
        """
        Get the ERA-Interim temperature, RH, and ozone
        @param era_files  ERA files and dates
        @param era_fileformat  ERA file format as prefix + date + postfix
        @param date  Desired date to interpolate to
        @param lat  Desired latitude to interpolate to
        @param lon  Desired longitude to interpolate to
        """

        # Load ERA-Interim nc file
        with Dataset(fname, "r", format="NETCDF4") as nci:
            # Variables, for convenience
            # 'longitude', 'latitude', 'level', 'time', 't', 'r'
            # Dimensions for t, r, o3, z:
            # (time, level, latitude, longitude)
            lats = nci["latitude"][:]
            lons = nci["longitude"][:]
            dates = num2date(
                nci["time"][:], nci["time"].units, nci["time"].calendar
            )

            # Check units
            # Expected variables and units
            expectedvars = {
                "level": "millibars",
                "t": "K",
                "r": "%",
                "z": "m**2 s**-2",
                "o3": "kg kg**-1",
            }
            # "u": "m s**-1",
            # "v": "m s**-1",

            for exptd in expectedvars:
                if exptd not in nci.variables.keys():
                    print(f"For file {fname}")
                    raise NameError(
                        "Expected variable in ERA file missing: " + exptd
                    )
                if nci[exptd].units != expectedvars[exptd]:
                    msg1 = "Expected units of " + expectedvars[exptd] + " for "
                    msg2 = exptd + ", but got " + nci[exptd].units
                    raise ValueError(msg1 + msg2)

            tmat = nci["t"][:]
            o3mat = nci["o3"][:]
            rhmat = nci["r"][:]
            zmat = nci["z"][:] / 9.80665 / 1000
            # umat = nci["u"][:]
            # vmat = nci["v"][:]

            temp = map_interp(tmat, date, lat, lon, dates, lats, lons)
            ozone = map_interp(o3mat, date, lat, lon, dates, lats, lons)
            rhw = map_interp(rhmat, date, lat, lon, dates, lats, lons)
            height = map_interp(zmat, date, lat, lon, dates, lats, lons)
            # uwind = map_interp(umat, date, lat, lon, dates, lats, lons)
            # vwind = map_interp(vmat, date, lat, lon, dates, lats, lons)

            self.P = nci["level"][:]
            self.T = temp
            self.o3 = ozone * 1000  # kg/kg * 1000 g/kg => g/kg
            self.rh = rhw
            self.z = height
            # self.uwind = uwind
            # self.vwind = vwind

    def tack_on_60km(self, hts: npt.NDArray[float], temps: npt.NDArray[float]):
        """
        Era ends just before 60 km, so tack on a few last values
        @param hts  Heights corresponding to temperatures to tack on
        @param temps  Temperatures to tack on
        Notes:
          For pressure, use the hypsometric equation.
          For temperature, this is necessary for because carbonTracker
          does not have the vertical resolution to make the interpolation
          work well.
          For rh, extrapolate
          For ozone, extrapolate
        """
        alti = np.array([self.z[-1], hts[-1]])
        tempi = np.array([self.T[-1], temps[-1]])
        rhi = np.array([self.rh[-1], self.rh[-1]])

        # Get the final pressures
        ptop = hypsometric(alti * 1000, tempi, rhi, self.P[-1])

        self.T = np.hstack([self.T, np.interp(60, alti, tempi)])
        self.rh = np.hstack([self.rh, self.rh[-1]])
        self.P = np.hstack([self.P, ptop[-1]])
        self.o3 = np.hstack([self.o3, np.interp(60, self.z, self.o3)])
        self.z = np.hstack([self.z, 60.0])

        # Do not let the final delta T be greater than +/- 5 K
        # dtemp = self.T[-1] - self.T[-2]
        # if np.abs(dtemp) > 5:
        #     raise NameError("Debug this!")
        #     # Just set it to be 5 K different from previous
        #     self.T[-1] = 5 * np.sign(dtemp) + self.T[-2]


def get_ind_from_ord(dtimes, dtime):
    """
    Get the index to the date using ordinals
    @param dtimes  Datetimes to choose from
    @param dtime Desired date to interpolate to
    @return  index to dtimes
    """
    # ERA files are per day, so check if there is a file on this day
    ords = [dt.datetime.toordinal(x) for x in dtimes]
    ord0 = dt.datetime.toordinal(dtime)
    if not ord0 in ords:
        raise ValueError("ERA file missing for " + str(dtime))

    return ords.index(ord0)


# # # # # # #
def nearest(items, pivot):
    """
    Get closest date without going over
    @param items  Items to find minimum of
    @param pivot
    """
    return min(items, key=lambda x: abs(x - pivot))


# # # # # # #     Get file names     # # # # # # #
def get_files_dates(direc: str, fmt: str):
    """
    Get files within a directory with a given format,
    and get the dates using the format
    @param direc  Desired directory to get file names from
    @param fmt  Format of desired files from direc
    @return List of files
    @return List of datetime objects
    """
    files_n_dates = get_filenames_dates(direc, fmt)
    files = [x[0] for x in files_n_dates]
    dates = [x[1] for x in files_n_dates]
    return files, dates


# # # # # # #     Get file names     # # # # # # #
def get_files(direc: str, fmt: str):
    """
    Get all filenames having the specified length and prefix, and datetimes
    @param direc  Desired directory to get file names from
    @param fmt  Format of desired files from direc
    @return  A tuple of file names and datetimes
    """
    return get_filenames_dates(direc, fmt)


# def get_files(directory, sample_filename, i1, i2, fstr):
# files = os.listdir(directory)
# files.sort()
# # files = filter(lambda x: len(x)==len(sample_filename), files)
# files = [
#     x
#     for x in files
#     if len(x) == len(sample_filename) and (x[:6] == sample_filename[:6])
# ]
# file_n_date = map(lambda f: (f, datetime.strptime(f[i1:i2], fstr)), files)
# return tuple(file_n_date)


def map_interp(amat, date, lat, lon, dates, lats, lons):
    """
    Interpolate to lat/lon/date
    dims of amat must be as follows:
    amat = time, level, lat, lon
    e.g.   (8,    25,   90, 120)
    """

    # Use scipy.interpolate instead of numpy.interp because the former
    # can handle descending data
    flat = interpolate.interp1d(
        lats, np.arange(lats.shape[0]), assume_sorted=False
    )
    flon = interpolate.interp1d(
        lons, np.arange(lons.shape[0]), assume_sorted=False
    )
    ilat = flat(lat)
    ilon = flon(lon)

    ddate = [d - date for d in dates]
    dday = [d.days + d.seconds / 60 / 60 / 24 for d in ddate]
    itime = np.interp(0, dday, np.arange(len(dday)))

    # As arrays
    nlevel = amat.shape[1]
    ilevs = np.arange(nlevel)
    ilats = ilat * np.ones(nlevel)
    ilons = ilon * np.ones(nlevel)
    itimes = itime * np.ones(nlevel)

    dims = [[itimes], [ilevs], [ilats], [ilons]]
    return ndimage.map_coordinates(amat, dims)[0]


def get_surface_met_frei(surfmetFiles, sonde_date):
    """
    Get surface temperatures from file corresponding to date of
    current sounding, from the surface met data at Frei
    @param surfmetFiles  Frei met files directory and file names
    @param sonde_date  Desired date
    """
    # Get the filename that includes the sonde date
    if sonde_date < surfmetFiles.dates[0]:
        ism = 0
    elif sonde_date > surfmetFiles.dates[-1]:
        ism = -1
    else:
        ism = bisect_right(surfmetFiles.dates, sonde_date) - 1
    metfname = surfmetFiles.dir + surfmetFiles.files[ism]

    # Load in the filename
    dfm = pd.read_csv(metfname)

    # Convert the column with the date to datetime
    dtime = pd.to_datetime(dfm["Date"])

    # Get the index to the time within the met file
    # falling just before the time of interest
    itime = max(bisect_right(dtime, sonde_date) - 1, 0)

    # The Frei met temperatures are every 15 minutes, and they are at a
    # different altitude, so just use this index
    surf = {}
    surf["date"] = dfm["Date"][itime]
    surf["temp"] = dfm["temperature (C)"][itime] + 273.15
    surf["rhw"] = dfm["RH (%)"][itime]
    surf["slp"] = dfm["Sea level pressure (hPa)"][itime]
    surf["press"] = dfm["Pressure (hPa)"][itime]

    return surf


def get_surface_temp_frei(surfmetFiles, sonde_date):
    """
    Get surface temperatures from file corresponding to date of
    current sounding, from the surface met data at Frei
    @param surfmetFiles  Frei met files directory and file names
    @param sonde_date  Desired date
    """
    # Get the filename that includes the sonde date
    ism = bisect_right(surfmetFiles.dates, sonde_date) - 1
    metfname = surfmetFiles.dir + surfmetFiles.files[ism]

    # Load in the filename
    dfm = pd.read_csv(metfname)

    # Convert the column with the date to datetime
    dtime = pd.to_datetime(dfm["Date"])

    # Get the index to the time within the met file
    # falling just before the time of interest
    itime = bisect_right(dtime, sonde_date) - 1

    # The Frei met temperatures are every 15 minutes, and they are at a
    # different altitude, so just use this index
    tmet = dfm["temperature (C)"][itime] + 273.15

    return tmet


def get_surface_met_tarp(surfmetFiles, sonde_date):
    """
    Get surface temperatures from file corresponding to date of
    current sounding.
    """

    # TODO: fix
    ism = bisect_right(surfmetFiles.dates, sonde_date) - 1
    metfname = surfmetFiles.dir + surfmetFiles.files[ism]

    raise ValueError("Need to fix this for json file!")
    with Dataset(metfname, "r", format="NETCDF4") as nci:
        itime = date2index(sonde_date, nci.variables["time"], select="after")

        if (itime == 0) or (sonde_date == nci.variables["time"][itime]):
            tmet = nci.variables["temp_mean"][itime] + 273.15
        else:
            surf_time = num2date(
                nci.variables["time"][:], nci.variables["time"].units
            )
            ddate = surf_time[itime - 1 : itime + 1] - sonde_date
            dmin = [d.days * 24 * 60 + d.seconds / 60 for d in ddate]
            wts = np.flipud(np.abs(dmin))
            wts = wts / np.sum(wts)
            tmet = (
                wts[0] * nci.variables["temp_mean"][itime - 1]
                + wts[1] * nci.variables["temp_mean"][itime]
                + 273.15
            )
        # zmet = nci["alt"][:] / 1000

    return tmet


def get_surf_met(surfmetFiles, sonde_date, surfmet_type: str):
    """
    Get surface temperature from met data
    @param surfmetFiles  Frei met files directory and file names
    @param sonde_date  Desired date
    @return surface temperature
    """
    get_surface_met_types = {
        "tarp": get_surface_met_tarp,
        "frei": get_surface_met_frei,
    }
    get_surface_met_fun = get_surface_met_types[surfmet_type]
    surf = get_surface_met_fun(surfmetFiles, sonde_date)

    return surf


# def get_surf_temp(surfmetFiles, sonde_date, surfmet_type: str):
#     """
#     Get surface temperature from met data
#     @param surfmetFiles  Frei met files directory and file names
#     @param sonde_date  Desired date
#     @return surface temperature
#     """
#     get_surface_temp_types = {
#         "tarp": get_surface_temp_tarp,
#         "frei": get_surface_temp_frei,
#     }
#     get_surface_temp = get_surface_temp_types[surfmet_type]
#     tmet = get_surface_temp(surfmetFiles, sonde_date)

#     return tmet


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# class HMS_to_hours:
#     def __init__(self, default=float("nan")):
#         # self.first = None
#         self.default = default

#     def __call__(self, value):
#         value = value.decode("ascii", "ignore").strip()
#         if not value:  # no value
#             return self.default
#         try:  # specify input format here
#             current = datetime.strptime(value, "%H:%M:%S")
#         except ValueError:  # can't parse the value
#             return self.default
#         else:
#             return current
#             # if self.first is not None:
#             #    return (current - self.first).total_seconds()
#             # else:
#             #    self.first = current
#             #    return 0.0
