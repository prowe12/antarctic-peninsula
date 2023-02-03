#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 13:46:11 2023

@author: prowe

From Zhang et al 2011

"For all four datasets, SBIs are identified from temperature profile data 
using Kahl’s (1990) algorithm, which scans from the surface upward to 
500 hPa to find cases with temperature increasing with altitude. The 
inversion top is defined as the bottom of the first layer in which 
temperature decreases with altitude, but thin (<100 m) noninversion 
layers, with temperature decreasing with height, are ignored if they are 
embedded within a deeper inversion layer. Soundings are considered to be 
unsuitable for analysis if the surface level is missing, there are fewer 
than 10 upper-air data levels from surface to 500 hPa, or the temperature 
difference across the inversion exceeds 40 K (Serreze et al. 1992). This 
definition represents a true inversion layer (from the surface to the 
inversion top), which is different from defining inversion strength as a 
temperature difference between two prespecified levels or heights in the 
boundary layer (Hudson and Brandt 2005; Kay and Gettelman 2009; 
Pavelsky et al. 2011)."

"For each grid point or radiosonde station, we compute frequency of SBI 
occurrence ( f ) and average SBI depth and intensity for the overall dataset 
and separately for each season, for each month, and for 0000 and 1200 UTC
(61 h to account for the 1–2-h sounding duration). Following Kahl (1990) 
and Serreze et al. (1992), three- month seasons are defined as January–March 
(Arctic winter, Antarctic summer), July–September 
(Arctic summer and Antarctic winter), etc.

"In all comparisons of the various datasets, the models and reanalysis 
are sampled at 0000 and 1200 UTC to mimic the radiosondes. Comparisons of 
radiosonde data and the two climate models are for 1990–2007. Comparisons 
for radiosonde data and ERA-Interim, and analysis of radiosonde data alone, 
are for 1990–2009."

"""

from netCDF4 import Dataset, num2date
import numpy as np
import matplotlib.pyplot as plt
import cdsapi


class PassOnProfPlot:
    def __init__(self, plotfig, temp, press, height):
        pass

    def plot(self, label):
        pass

    def addlabels(self):
        pass

    def finalize(self):
        pass

    def addinv(self, iinv, iswitch):
        pass


class ProfPlotter:
    def __init__(self, plotfig, temp, press, height):
        self.plotfig = plotfig
        self.temp = temp
        self.press = press
        self.height = height

    def plot(self, label):
        if plotfig:
            plt.figure(figno)
            plt.subplot(121)
            plt.plot(self.temp, self.height, ".-", label=label)
            plt.subplot(122)
            plt.plot(self.temp, self.press)

    def addinv(self, iinv, iswitch):
        plt.subplot(121)
        plt.plot(self.temp[iinv], self.height[iinv], "go")
        plt.plot(self.temp[iswitch[0]], self.height[iswitch[0]], "k*")

    def finalize(self):
        if plotfig:
            plt.subplot(121)
            # plt.legend()
            plt.ylim([0, 6])
            plt.ylabel("Height (km)")
            plt.subplot(122)
            plt.ylim([500, 1000])
            plt.gca().invert_yaxis()
            plt.xlabel("Temperature (K)")
            plt.ylabel("Pressure (mb)")


def inversion_exists(temp):
    # ID if surface-based temperature inversion is present
    # True if temperature decreases between the first two
    # temperatures *or* the second two temperatures
    # the latter condition is added because the interpolation
    # to the surface temperature appears to be imperfect
    if temp[1] - temp[0] >= 0 or temp[2] - temp[1] >= 0:
        return True
    else:
        return False


def get_params(year):
    # Pressure levels
    press_levs = np.hstack(
        [
            np.arange(1000, 925 - 5, -25),
            np.arange(900, 450, -50),
        ]
    )
    press_levs_ints = list(press_levs)
    press_levs = [str(p) for p in press_levs_ints]

    north = -88.0
    south = -90.0
    west = 121.0  # -59.0
    east = 122.0  # -58.5
    return {
        "format": "netcdf",
        "product_type": "reanalysis",
        "variable": [
            "pressure",
            "Geopotential",
            "Temperature",
        ],
        "pressure_level": press_levs,
        "year": [year],
        "month": ["07"],
        "day": ["01", "02", "03"],
        "time": ["00:00", "12:00"],
        "grid": [0.5, 0.5],
        "area": [north, west, south, east],
    }


def get_era5(eradir):
    # Get the ERA5 data: one dataset per year

    dataset = "reanalysis-era5-pressure-levels"
    c_era = cdsapi.Client()

    # retrieve the path to the file
    year = "2022"
    prefix = "inversion"
    outfile = prefix + year + ".nc"
    params = get_params(year)
    fid = c_era.retrieve(dataset, params)

    # download the file and save it to the output directory
    fid.download(eradir + outfile)


def get_stats(temp, height, plotProf):
    dtemp = np.diff(temp)
    iinv = np.where(dtemp > 0)[0] + 1
    switch = np.diff(np.sign(dtemp))
    iswitch = np.where(switch != 0)[0] + 1

    # However, thin (<100 m) noninversion layers, with temperature decreasing
    # with height, are ignored if they are embedded within a deeper inversion
    # layer. There should be none here, because height differences should
    # all be > 100 m.  Following is a catch, just in case
    if not np.all(np.diff(height)[1:] > 0.1):
        msg = (
            "Some height differences are < 100 m, so a thin noninversion "
            + "that should be ignored (counted as part of the inversion) "
            + "may exist within in the inversion"
        )
        # The following code could be used to identify thin non-inversions
        # that should be ignored, if the above error is ever tripped
        # The inversion top is defined as the bottom of the first layer in which
        # temperature decreases with altitude

        # # Progress through the noninversion depths and include them while
        # # they are less than 100 m
        # currlev = 1  # start with 0 and 1 switch
        # while len(iswitch) > currlev:
        #     noninvdepth = height[iswitch[currlev]] - height[iswitch[currlev - 1]]
        #     if 1000 * noninvdepth > 100:
        #         # We are out of the inversion, so quit
        #         break
        #     currlev += 1
        raise ValueError(msg)

    try:
        invdepth = height[iswitch[0]]
    except:
        print("pause here")
    invstrength = temp[iswitch[0]] - temp[0]
    plotProf.addinv(iinv, iswitch)

    return invdepth, invstrength


def interp_to_zsrf(zin, tin, pin, zsrf):
    tsrf = np.interp(zsrf, zin, tin)
    psrf = np.interp(zsrf, zin, pin)
    iabove = np.where(zin > zsrf)[0]
    znew = np.hstack([zsrf, zin[iabove]])
    tnew = np.hstack([tsrf, tin[iabove]])
    pnew = np.hstack([psrf, pin[iabove]])
    return znew, tnew, pnew


def get_stats_from_ncid(ncid, itime, ilat, ilon, zsrf):
    zfac = 1 / 9.80665 / 1000
    geop, temp, press = interp_to_zsrf(
        ncid["z"][itime, :, ilat, ilon][i500],
        ncid["t"][itime, :, ilat, ilon][i500],
        ncid["level"][i500],
        zsrf,
    )
    # TODO: convert to geometric height!
    # alt=Re*h/(Re−h)
    height = zfac * geop

    plotProf = PlotProf(figno, temp, press, height)

    if inversion_exists(temp):
        plotProf.plot(filename)
        invdepth, invstrength = get_stats(temp, height, plotProf)
        plotProf.finalize()
        return True, invdepth, invstrength, temp[0]

    # Plot and move on to the next time period
    plotProf.plot(filename + ": none")
    plotProf.finalize()
    return False, None, None, temp[0]


def get_surface_heights():
    direc = "/Users/prowe/Sync/measurements/Antarctica/"
    filename = "geo_1279l4_0.1x0.1.grib2_v4_unpack.nc"

    with Dataset(direc + filename, "r") as ncid:
        ilats = np.where(ncid["latitude"][:] <= -88)[0]
        srfht = {
            "lats": ncid["latitude"][ilats],
            "lons": ncid["longitude"][:],
            "z": ncid["z"][0, ilats, :],
        }

    return srfht


def get_surface_height(srfhts, lat, lon):
    ilat = np.where(srfhts["lats"] == lat)[0]
    if len(ilat) != 1:
        raise ValueError("Bad value for lat")
    ilon = np.where(srfhts["lons"] == lon)[0]
    if len(ilon) != 1:
        raise ValueError("Bad value for lat")
    return srfhts["z"][ilat[0], ilon[0]]


plotfig = True
if plotfig:
    PlotProf = ProfPlotter
else:
    PlotProf = PassOnProfPlot


eradir = "/Users/prowe/Sync/measurements/Antarctica/era5/inversions/"
# TODO: Uncomment to get era data
# get_era5(eradir)

# Get the surface geopotential
srfhts = get_surface_heights()

filenames = ["inversion2022.nc"]

# Get info from the first file that doesn't change
filename = filenames[0]
with Dataset(eradir + filename, "r") as ncid:
    levels = ncid["level"][:]
    i500 = np.where(levels >= 500)[0]
    nlats = ncid.dimensions["latitude"].size
    nlons = ncid.dimensions["longitude"].size

    # Data for output
    lats = ncid["latitude"][:].data
    lons = ncid["longitude"][:].data

    # To get the time:
    # num2date(ncid['time'][:],ncid['time'].units)


# Preallocate variables
press = levels[i500]
depthsum = np.zeros([nlats, nlons])
strengthsum = np.zeros([nlats, nlons])
tsurfsum = np.zeros([nlats, nlons])
ninversions = np.zeros([nlats, nlons])
ncases = np.zeros([nlats, nlons])


figno = 2
if plotfig:
    plt.figure(figno)
    plt.clf()

for filename in filenames:
    depths = np.zeros([nlats, nlons])
    with Dataset(eradir + filename, "r") as ncid:
        for ilat in range(nlats):
            for ilon in range(nlons):
                ncases0 = 0
                ninversion = 0
                depth = []
                strength = []
                tsurf = []

                # Get the surface height for this lat and lon
                zsrf = get_surface_height(srfhts, lats[ilat], lons[ilon])

                for itime in range(ncid.dimensions["time"].size):
                    ncases0 += 1
                    isinv, invdepth, invstrength, tsrf = get_stats_from_ncid(
                        ncid, itime, ilat, ilon, zsrf
                    )
                    tsurf.append(tsrf)
                    if isinv:
                        ninversion += 1
                        depth.append(invdepth)
                        strength.append(invstrength)

                # QC
                if (len(depth) != ninversion) or (len(strength) != ninversion):
                    raise ValueError(
                        "Number of inversion and depths/strengths differs!"
                    )

                ncases[ilat, ilon] = ncases0
                ninversions[ilat, ilon] = ninversion
                depthsum[ilat, ilon] = np.sum(depth)
                strengthsum[ilat, ilon] = np.sum(strength)
                tsurfsum[ilat, ilon] = np.sum(tsurf)


frequency = ninversions / ncases
tsurface = tsurfsum / ncases
depth = depthsum / ninversions
strength = strengthsum / ninversions

plt.figure(3)
plt.clf()
plt.subplot(3, 1, 1)
plt.plot(lats, frequency, ".-")
plt.ylabel("Frequency")
plt.subplot(3, 1, 2)
plt.plot(lats, depth, ".-")
plt.ylabel("Inversion depth (m)")
plt.subplot(3, 1, 3)
plt.plot(lats, strength, ".-")
plt.ylabel("Inversion strength (K)")

plt.figure(4)
plt.clf()
plt.plot(lats, tsurface - 273.15, ".")

print(f"The frequency is: {frequency}")
# print(f"The mean depth is: {round(1000*np.mean(depth))} m")
# print(f"The mean strength is: {round(np.mean(strength),2)} K")


# for filename in filenames:
#     ncid = Dataset(eradir + filename, "r")
#     height = zfac * ncid["z"][0, :, 0, 0]
#     temp = ncid["t"][0, :, 0, 0]
#     levels = ncid["level"][:]
#     plt.subplot(121)
#     plt.plot(temp, height, label=filename)
#     plt.subplot(122)
#     plt.plot(temp, levels)
#     ncid.close()
# plt.subplot(121)
# plt.legend()
# plt.ylim([0, 13])
# plt.subplot(122)
# plt.ylim([100, 1002])
# plt.gca().invert_yaxis()
