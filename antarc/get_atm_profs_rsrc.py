#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 14:40:05 2018

@author: prowe
"""

# .. Built-in modules
import scipy.io.netcdf as netcdf
import numpy as np
from netCDF4 import Dataset
from netCDF4 import num2date, date2num, date2index
import scipy.ndimage as ndimage
import os
import sys
from datetime import datetime, timedelta
from bisect import bisect, bisect_right
from scipy import interpolate
import time
import datetime as dt
import xarray as xr

from matplotlib.dates import datestr2num
import pandas as pd


# My modules
from antarc.humidRS80 import humidRS80
from antarc.hypsometric import hypsometric, hypsometric_for_z
from antarc.getfilenames import get_filenames_dates


# # # # # # #     Get file names     # # # # # # #
class FileInfo:
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


class HMS_to_hours:
    def __init__(self, default=float("nan")):
        # self.first = None
        self.default = default

    def __call__(self, value):
        value = value.decode("ascii", "ignore").strip()
        if not value:  # no value
            return self.default
        try:  # specify input format here
            current = datetime.strptime(value, "%H:%M:%S")
        except ValueError:  # can't parse the value
            return self.default
        else:
            return current
            # if self.first is not None:
            #    return (current - self.first).total_seconds()
            # else:
            #    self.first = current
            #    return 0.0


# # # # # # #     Radiosonde     # # # # # # #
class Sonde:
    __slots__ = ("file", "date", "time", "z", "T", "P", "rh", "h2o")

    def __init__(
        self, sonde_dir: str, sonde_file: str, sonde_date: dt.datetime
    ):

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

        self.time = hours  # nc.variables['time'][:] + 0
        self.T = temper  #
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

    def qc(self):
        """Quality control"""
        if any(np.diff(self.P) >= 0):
            raise ValueError("One or more pressures increase with altitude")

        if any(np.diff(self.z) <= 0):
            raise ValueError("One or more altitudes decrease with time")

        if (
            any(np.isnan(self.time))
            or any(np.isnan(self.z))
            or any(np.isnan(self.T))
            or any(np.isnan(self.P))
            or any(np.isnan(self.h2o))
            or any(np.isnan(self.rh))
        ):
            raise ValueError("One or more values is nan")
            # "file", "date", "time", "z", "T", "P", "rh", "h2o")

        # Fix any missing heights - MUST USE hypsometric_for_z
        # BECAUSE SOME VALUES DECREASE !!!!
        imissing = np.where(self.z == -999)[0]
        if np.any(imissing):
            raise ValueError("This should not happen; see code")
            # alt0 = self.z[imissing[0] - 1]
            # zfill = hypsometric_for_z(
            #     self.P[imissing], self.T[imissing], self.rh[imissing], alt0
            # )
            # self.z[imissing] = zfill
            # plt.plot(temper,alt)


# # # # # # #     Station flask CO2     # # # # # # #

# .. Load station co2
class CO2stationData:
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

        # .. Get datetime
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
        with Dataset(co2Files.dir + co2Files.files[ibef], "r") as nc:
            dates = num2date(
                nc.variables["time"][:], nc.variables["time"].units
            )
        if min(dates) > date:
            ibef -= 1
            iaft -= 1

        # .. Open the netcdf co2 file and get the data
        with Dataset(co2Files.dir + co2Files.files[ibef], "r") as nc:
            """
            # nc['co2'].shape
            # Out[32]: (8,    25,   90, 120)
            #         time, level, lat, lon
            """
            Nlevel = nc["level"].shape[0]
            Nbound = nc["boundary"][:].shape[0]
            lats = nc["latitude"][:]
            lons = nc["longitude"][:]
            dates = num2date(
                nc.variables["time"][:], nc.variables["time"].units
            )

            """
            # Verify it worked
            print('units = %s, values = %s' % (nc.variables['time'].units, 
                                               nc.variables['time'][:]))
            print(nc['time_components'][:])
            print(dates)
            print([date.strftime('%Y-%m-%d %H:%M:%S') for date in dates])
            """

            if max(dates) < date:
                # Glom on the last with the first
                co2mat = np.zeros((2, Nlevel, 90, 120))
                Tmat = np.zeros((2, Nlevel, 90, 120))
                Pmat = np.zeros((2, Nbound, 90, 120))
                zmat = np.zeros((2, Nbound, 90, 120))
                co2mat[0, :, :, :] = nc["co2"][-1, :, :, :]
                Tmat[0, :, :, :] = nc["temperature"][-1, :, :, :]
                Pmat[0, :, :, :] = nc["pressure"][-1, :, :, :]
                zmat[0, :, :, :] = nc["gph"][-1, :, :, :]

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
                    Tmat[1, :, :, :] = nc_aft["temperature"][0, :, :, :]
                    Pmat[1, :, :, :] = nc_aft["pressure"][0, :, :, :]
                    zmat[1, :, :, :] = nc_aft["gph"][0, :, :, :]

                    co2 = map_interp(co2mat, date, lat, lon, dates, lats, lons)
                    T = map_interp(Tmat, date, lat, lon, dates, lats, lons)
                    P = map_interp(Pmat, date, lat, lon, dates, lats, lons)
                    z = map_interp(zmat, date, lat, lon, dates, lats, lons)

                    del dates, dates_aft, co2mat, Tmat, Pmat, zmat

            else:
                co2mat = nc["co2"][:]
                Tmat = nc["temperature"][:]
                Pmat = nc["pressure"][:]
                zmat = nc["gph"][:]

                co2 = map_interp(co2mat, date, lat, lon, dates, lats, lons)
                T = map_interp(Tmat, date, lat, lon, dates, lats, lons)
                P = map_interp(Pmat, date, lat, lon, dates, lats, lons)
                z = map_interp(zmat, date, lat, lon, dates, lats, lons)

                del dates, co2mat, Tmat, Pmat, zmat

            del Nlevel, Nbound, lats, lons

            self.co2 = co2
            self.T = T
            self.P = P / 100
            self.zbnd = z / 1000

        self.z = self.zbnd[:-1] + np.diff(self.zbnd) / 2


# # # # # # #     Profile     # # # # # # #
class Prof:
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

        # .. Guess values (based on fit to a single measurement,
        #    needs improvement)
        self.f113 = 0.00008 * np.ones(nlyr)  # 800-828, 870-930,
        self.f11 = 0.00026 * np.ones(nlyr)  # 830-860, 1060-1107
        self.hno3 = 0.002 * np.ones(nlyr)  # 860-920
        self.f12 = 0.00053 * np.ones(nlyr)  # 867-937, 1080-1177

        # .. Units for ozone are jchar = C for g/kg
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
        # .. Get the old and new values to work with
        new_val = getattr(new, attr_name)
        new_val = np.interp(self.z, new.z, new_val)

        setattr(self, attr_name, new_val)

        if new.z[0] > self.z[0]:
            print(
                "Warning: Extrapolating ",
                attr_name,
                "for heights: ",
                self.z[self.z < new.z[0]],
            )
        if new.z[-1] < self.z[-1]:
            print(
                "Warning: Extrapolating ",
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

        new_val = getattr(new, attr_name)
        new_val = np.interp(self.z, new.z, new_val)
        new_val = surf_value * new_val / new_val[0]

        setattr(self, attr_name, new_val)

    def set_upper(self, attr_name, new):
        # .. Get the old and new values to work with
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
        # .. Get the old and new values to work with
        old_val = getattr(self, attr_name)
        new_val = getattr(new, attr_name)

        # .. Do the spline interpolation
        tck = interpolate.splrep(new.z, new_val)
        new_val = interpolate.splev(self.z, tck)

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
        zi = self.z[inds]
        # print('zi[0] =', zi[0], ' in get_atm_profs_KGI_rsrc.py, line 197')
        # print('zi[-1] =', zi[-1], ' in get_atm_profs_KGI_rsrc.py')
        Ti = self.T[inds]
        rhi = 0 * zi
        Ptop = hypsometric(zi * 1000, Ti, rhi, self.P[inds[0]])

        old_val[inds] = Ptop
        setattr(self, attr_name, old_val)

    def error_check(self, era):
        # Check for NaNs
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
                print("Warning: one or more h2o is NaN, filling with era: ")
                print(self.h2o)

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

        # .. Check for heights (z) that do not monotonically increase or
        #    pressures (P) that do not monotonically decrease
        #    to within 3 significant figures
        zmr = np.round(self.z, 3)
        pmr = np.round(self.P, 3)
        if np.any(np.diff(zmr) <= 0):
            raise ValueError("zs do not increase, or not enough.")

        if (np.any(np.diff(pmr) >= 0)) or (np.any(pmr == 0)):
            raise ValueError(
                "Problem with pressure: not decreasing enough "
                "or 0 detected."
            )

        # .. Check for changes in temperature greater than 8 in lower 20 km
        if np.any(np.abs(np.diff(self.T[self.z <= 20])) > 8):
            max_dT = round(np.max(np.diff(self.T[self.z <= 20])), 1)
            print(
                "Warning: Temperature differences in lower 20 km are up "
                + "to "
                + str(max_dT)
                + " K but should be < 8 K "
                + "(getAtmProfs_KGI_rsrc line 465."
            )

        # .. Check for changes in temperature greater than 8 K anywhere
        if np.any(np.abs(np.diff(self.T)) > 12):
            max_dT = round(np.max(np.abs(np.diff(self.T))), 1)
            print(self.T)
            import matplotlib.pyplot as plt

            plt.figure()
            plt.plot(self.T, self.z, ".-")
            raise NameError(
                "Temperature differences above 20 km are up to "
                + str(max_dT)
                + " K but should be < 12 K."
            )

    def set_surf_T(self, surfmetFiles, sonde_date):
        z0 = self.z[0]

        zmet, Tmet = get_surface_T(surfmetFiles, sonde_date)

        if zmet < z0:
            raise NameError(
                "Met height is below sounding surface height. "
                + "Add code to interpolate to sounding surface."
            )
        elif zmet > z0:
            print(
                "Warning: Met height is above sounding surface height "
                + "but I am pretending surface heights are same."
            )
            # raise NameError('Met height is above sounding surface height')
        elif zmet - self.z != 0:
            raise NameError("Bad value for met height")

        # .. Interpolate
        iht = np.where(self.z <= 0.6)[0]
        Tnew = np.interp(
            self.z[iht], [self.z[0], self.z[iht[-1]]], [Tmet, self.T[iht[-1]]]
        )

        f = np.linspace(1, 0, len(iht))
        self.T[iht] = f * Tnew + (1 - f) * self.T[iht]

    """
    def set_upper_P_carbonTracker(self, new, old_attr_name, new_attr_name):
        # .. Get the old and new values to work with
        old_val = getattr(self, old_attr_name)
        new_val = getattr(new, new_attr_name)
        
        # .. Do the spline interpolation 
        tck = interpolate.splrep(new.zbnd, new_val)
        new_val = interpolate.splev(self.z, tck)
        
        inds = np.intersect1d(np.where(self.z >= new.zbnd[0])[0],
                              np.where(self.z <= new.zbnd[-1])[0])
        inds = np.intersect1d(inds, np.where(np.isnan(old_val)))
        old_val[inds] = new_val[inds]
        #setattr(self, old_attr_name, old_val)
    """

    def write_to_netcdf_file(self, outputfile):
        with Dataset(outputfile, "w", format="NETCDF4_CLASSIC") as nc:
            # .. Create Dimensions
            nc.createDimension("level", len(self.z))
            nc.createDimension("time", 1)
            nc.createDimension("const", 1)
            # lat = nc.createDimension('lat', 1)
            # lon = nc.createDimension('lon', 1)

            # .. Create variables
            nc_time = nc.createVariable("time", np.float64, ("time",))
            nc_z = nc.createVariable("z", np.float32, ("level"))
            nc_T = nc.createVariable("T", np.float32, ("level"))
            nc_P = nc.createVariable("P", np.float32, ("level"))
            nc_rh = nc.createVariable("rh", np.float32, ("level"))
            nc_h2o = nc.createVariable("h2o", np.float32, ("level"))
            nc_co2 = nc.createVariable("co2", np.float32, ("level"))
            nc_o3 = nc.createVariable("o3", np.float32, ("level"))
            nc_f11 = nc.createVariable("f11", np.float32, ("level"))
            nc_f12 = nc.createVariable("f12", np.float32, ("level"))
            nc_f113 = nc.createVariable("f113", np.float32, ("level"))
            nc_hno3 = nc.createVariable("hno3", np.float32, ("level"))
            nc_model_extra = nc.createVariable(
                "model_extra", np.int8, ("const")
            )

            # .. Global attributes
            nc.filename = self.file
            nc.history = "Created " + time.ctime(time.time())

            # .. Variable attributes
            nc_time.units = "hours since 0001-01-01 00:00:00"
            nc_time.calendar = "gregorian"
            nc_z.units = self.units["z"]
            nc_T.units = self.units["T"]
            nc_P.units = self.units["P"]
            nc_rh.units = "%"
            nc_h2o.units = self.units["h2o"]
            nc_co2.units = self.units["co2"]
            nc_o3.units = self.units["o3"]
            nc_hno3.units = self.units["hno3"]
            nc_f11.units = self.units["f11"]
            nc_f12.units = self.units["f12"]
            nc_f113.units = self.units["f113"]

            # .. Assign values
            nc_time[:] = date2num(
                self.date, units=nc_time.units, calendar=nc_time.calendar
            )
            nc_z[:] = self.z
            nc_T[:] = self.T
            nc_P[:] = self.P
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
    # # # #
    # Uncomment and put text below this line in class
    # co2 = CO2(co2_file_palmer, co2_file_drake, sonde.date, lat, lon)
    #
    import calendar

    def toTimestamp(d):
        return calendar.timegm(d.timetuple())

    def interp_co2_to_date(sonde_date, station_date, station_co2, station_lat):
        iaft = bisect(station_date, sonde_date)
        ibef = iaft - 1

        # .. Sometimes the "before" date is actually slightly after
        #    the first date in the file. In this case, decrease ibef and iaft
        if min(station_date) > sonde_date:
            ibef -= 1
            iaft -= 1

        # .. Get the Drake Passage CO2 value

        x = np.array(
            [toTimestamp(station_date[ibef]), toTimestamp(station_date[iaft])]
        )
        co2 = np.interp(
            toTimestamp(sonde_date), x, station_co2[ibef : iaft + 1]
        )
        lat = np.interp(
            toTimestamp(sonde_date), x, station_lat[ibef : iaft + 1]
        )

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
    if use_palmer == True and use_drake == True:
        fco2 = interpolate.interp1d(
            [lat_drake0, lat_palmer0], [co2_drake0, co2_palmer0]
        )
        return fco2(sonde_lat)
    elif use_palmer == True:
        return co2_palmer0
    elif use_drake == True:
        return co2_drake
    else:
        raise ValueError("No co2 data to use")


class Era:
    __slots__ = ("P", "T", "o3", "rh", "z")

    def __init__(self, eraFiles, era_fileformat, date, lat, lon):
        """
        Get the ERA-Interim temperature, RH, and ozone
        @param eraFiles  Class with ERA dates, directory, and files
        @param era_fileformat  ERA file format as prefix + date + postfix
        @param date  Desired date to interpolate to
        @param lat  Desired latitude to interpolate to
        @param lon  Desired longitude to interpolate to
        """

        # ERA files are per day, so check if there is a file on this day
        ords = [dt.datetime.toordinal(x) for x in eraFiles.dates]
        ord0 = dt.datetime.toordinal(date)
        if not ord0 in ords:
            raise ValueError("ERA file missing for " + str(date))

        idate = ords.index(ord0)

        # Load ERA-Interim nc file
        fname = eraFiles.files[idate]
        with Dataset(eraFiles.dir + fname, "r", format="NETCDF4") as nc:
            # Variables, for convenience
            # 'longitude', 'latitude', 'level', 'time', 't', 'r'
            # Dimensions for t, r, o3, z:
            # (time, level, latitude, longitude)
            lats = nc["latitude"][:]
            lons = nc["longitude"][:]
            dates = num2date(
                nc["time"][:], nc["time"].units, nc["time"].calendar
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

            for exptd in expectedvars:
                if exptd not in nc.variables.keys():
                    raise NameError(
                        "Expected variable in ERA file missing: " + exptd
                    )
                if nc[exptd].units != expectedvars[exptd]:
                    msg1 = "Expected units of " + expectedvars[exptd] + " for "
                    msg2 = exptd + ", but got " + nc[exptd].units
                    raise ValueError(msg1 + msg2)

            tmat = nc["t"][:]
            o3mat = nc["o3"][:]
            rhmat = nc["r"][:]
            zmat = nc["z"][:] / 9.80665 / 1000

            temp = map_interp(tmat, date, lat, lon, dates, lats, lons)
            ozone = map_interp(o3mat, date, lat, lon, dates, lats, lons)
            rhw = map_interp(rhmat, date, lat, lon, dates, lats, lons)
            height = map_interp(zmat, date, lat, lon, dates, lats, lons)

            self.P = nc["level"][:]
            self.T = temp
            self.o3 = ozone * 1000  # kg/kg * 1000 g/kg => g/kg
            self.rh = rhw
            self.z = height

        # Interpolate to location of Escudero
        # import pandas as pd
        # import pygmt
        # from urllib.request import urlopen

    # def read_data_from_file(self, filename):
    #     """
    #     Read in the data from the ERA5 file
    #     """
    #     # load into memory
    #     with xr.open_dataset(filename) as ds:
    #         points = pd.DataFrame(
    #             data={
    #                 "lon": [LONGITUDE],
    #                 "lat": [LATITUDE],
    #             }
    #         )
    #         print("coordinates dataframe")
    #         print(points)

    #         temp = []
    #         for ilevel in range(ds.dims["level"]):

    #             grid = xr.DataArray(
    #                 data=ds.t.values[0, ilevel, :, :],
    #                 dims=["lat", "lon"],
    #                 coords={
    #                     "lat": ds.latitude.values,
    #                     "lon": ds.longitude.values,
    #                 },
    #             )

    #             track = pygmt.grdtrack(
    #                 points=points, grid=grid, newcolname="sampled_data1"
    #             )
    #             temp.append(track["sampled_data1"].values[0])

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(temp, ds.level)
    # plt.ylim([1000, 0])
    # plt.ylabel("P (hPa)")
    # plt.xlabel("Temperature (K)")

    # Note: direct from URL method:
    # with urlopen(fl.location) as f:
    #    ds = xr.open_dataset(f.read())

    def tack_on_60km(self, carbonTracker):
        """
        Era ends just before 60 km, so tack on a few last values
          For pressure, use the hypsometric equation.
          For temperature, this is necessary for because carbonTracker
          does not have the vertical resolution to make the interpolation
          work well.
          For rh, extrapolate
          For ozone, extrapolate
        """
        zi = np.array([self.z[-1], carbonTracker.z[-1]])
        Ti = np.array([self.T[-1], carbonTracker.T[-1]])
        rhi = np.array([self.rh[-1], self.rh[-1]])
        # print('zi[0] =', zi[0], ' in get_atm_profs_KGI_rsrc.py, line 517')
        # print('zi[-1] =', zi[-1], ' in get_atm_profs_KGI_rsrc.py')
        Ptop = hypsometric(zi * 1000, Ti, rhi, self.P[-1])

        self.T = np.hstack([self.T, np.interp(60, zi, Ti)])
        self.rh = np.hstack([self.rh, self.rh[-1]])
        self.P = np.hstack([self.P, Ptop[-1]])
        self.o3 = np.hstack([self.o3, np.interp(60, self.z, self.o3)])
        self.z = np.hstack([self.z, 60.0])

        # Do not let the final delta T be greater than +/- 5 K
        # dtemp = self.T[-1] - self.T[-2]
        # if np.abs(dtemp) > 5:
        #     raise NameError("Debug this!")
        #     # Just set it to be 5 K different from previous
        #     self.T[-1] = 5 * np.sign(dtemp) + self.T[-2]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # #  Get closest date without going over
def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))


# # # # # # #     Get file names     # # # # # # #
def get_files(direc, fmt):
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

    dd = [d - date for d in dates]
    dday = [d.days + d.seconds / 60 / 60 / 24 for d in dd]
    itime = np.interp(0, dday, np.arange(len(dday)))

    # As arrays
    nlevel = amat.shape[1]
    ilevs = np.arange(nlevel)
    ilats = ilat * np.ones(nlevel)
    ilons = ilon * np.ones(nlevel)
    itimes = itime * np.ones(nlevel)

    dims = [[itimes], [ilevs], [ilats], [ilons]]
    return ndimage.map_coordinates(amat, dims)[0]


def get_surface_T(surfmetFiles, sonde_date):
    """
    Get surface temperatures from file corresponding to date of
    current sounding.
    """

    ism = bisect_right(surfmetFiles.dates, sonde_date) - 1
    metfname = surfmetFiles.dir + surfmetFiles.files[ism]
    with Dataset(metfname, "r", format="NETCDF4") as nc:
        itime = date2index(sonde_date, nc.variables["time"], select="after")

        if (itime == 0) or (sonde_date == nc.variables["time"][itime]):
            Tmet = nc.variables["temp_mean"][itime] + 273.15
        else:
            surf_time = num2date(
                nc.variables["time"][:], nc.variables["time"].units
            )
            dd = surf_time[itime - 1 : itime + 1] - sonde_date
            dmin = [d.days * 24 * 60 + d.seconds / 60 for d in dd]
            wt = np.flipud(np.abs(dmin))
            wt = wt / np.sum(wt)
            Tmet = (
                wt[0] * nc.variables["temp_mean"][itime - 1]
                + wt[1] * nc.variables["temp_mean"][itime]
                + 273.15
            )
        zmet = nc["alt"][:] / 1000

    return zmet, Tmet


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
