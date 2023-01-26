#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 14:40:05 2018

@author: prowe

Copyright 2019 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

Purpose: Create atmospheric profiles for a layered model atmosphere for
running LBLRTM and DISORT, for King George Island, Antarctica.

The atmospheric layer boundaries should be:
1) Thin enough that te temperature differential across layers is <= 10 K,
   preferrable <= 5K. From DISORT.doc:
      "The computational layering is usually constrained by the problem,
      in the sense that each computational layer must be reasonably
      homogeneous and not have a temperature variation of more than
      about 5-10 K across it (if thermal sources are important)."
2) Thick enough that the gaseous optical depth is >= 1e-5 if DISORT is run
   in single precision mode (currently the case). Otherwise either the
   optical depth has to be bumped up to 1e-5 to prevent errors or the
   atmosphere above has to be removed from the calculation.
3) The compromise is made by having layers get thicker going up in the
   atmosphere. In CLARRA, if a tropospheric layer has an optical depth
   < 1e-5, it is bumped up to 1e-5, whereas if a layer above the troposphere
   has an opticla depth < 1e-5, the entire atmosphere above is chopped off.
   This is obviously not the best way to do it, better would be to change
   DISORT to double precision.

Output: A netcdf file (optionally also a pickle file and figures) with:
    Dimensions: level
    Variables:
        time in hours since 0001-01-01 00:00:00 (scalar)
        z, T, P, rh, h2o, co2, o3, f11, f12, f113, hno3, model_extra (x level)
        model_extra: model number for unspecified trace gas amounts (scalar)
    Attributes:
        data (e.g. concentration)
        units
        long_name (to be added later)


Variables:
 - The class instance prof (see get_atm_prof_KGI_rsrc), which is
   sorted out and written to the netcdf file. The fields of prof are as
   follows, where for gases, the naming convention is to use the lower
   case and altitude, pressure, and temperature are z, P, and T.
   It is easy to add more gases, just add the name (with concentration)
   and the units.
     file: filename
     date
     time
     z: Altitude (km)
     P: Pressure (mb)
     T: Temperature (K)
     h2o
     co2
     rh: Relative humidity (for reference, h2o is used)
     h2o
     co2
     o3
     f11
     f12
     f113
     hno3
     units
     model_extra
 - For other classes and variables, see get_atm_prof_KGI_rsrc






"""

# .. Built-in modules
import numpy as np
import matplotlib.pyplot as plt
import pickle
import copy
import datetime as dt


# .. My modules
from get_atm_profs_KGI_rsrc import get_files
from get_atm_profs_KGI_rsrc import FileInfo
from get_atm_profs_KGI_rsrc import Sonde
from get_atm_profs_KGI_rsrc import Era
from get_atm_profs_KGI_rsrc import Prof
from get_atm_profs_KGI_rsrc import get_co2
from get_atm_profs_KGI_rsrc import CO2stationData
from get_atm_profs_KGI_rsrc import CarbonTracker

# from get_atm_profs_KGI_rsrc import get_surface_T


# # # # # # # # #    INPUTS   # # # # # # # # # #
# .. Directories
main_dir = "/Users/prowe/Projects/NSF_AERI_cloud_CRI/McMurdo/case_studies/"
yopp_dir = "/Users/prowe/Projects/YOPP/"
co2_dir = yopp_dir + "Bromwich_etal_2020/work/co2/"
sonde_dir = yopp_dir + "data_denial_experiment/radiosondes_for_SRC_"
era_dir = yopp_dir + "Bromwich_etal_2020/work/era_interim/"
out_dir = yopp_dir + "Bromwich_etal_2020/profFiles/"
surfmet_dir = main_dir + "met/"
fig_dir = out_dir + "/figures/"

# .. Filenames
era_filename = "era_interim_KGI_201912.nc"
co2_file_palmer = co2_dir + "co2_psa_surface-flask_1_ccgg_event.txt"
co2_file_drake = co2_dir + "co2_drp_surface-flask_1_ccgg_event.txt"
ct_file = "CT2016.molefrac_glb3x2_2016-12-01.nc"
# surfmet_file = 'awrmetM1.b1.20160208.000000.cdf'

# date indices to files
iera_date1 = 16
iera_date2 = 22

# .. Input parameters
location = "Escudero"  #'Sejong' #

# .. Desired date range
start_date = dt.datetime(2018, 12, 4, 0, 0)
end_date = dt.datetime(2018, 12, 31, 0, 0)

# # # # # # # # # # # # # # # # # # # # # # # #


# # # # # #    START MAIN CODE   # # # # # #

# . Set up directories, files, etc
sonde_dir = sonde_dir + location + "/"
sonde_file = location + "_2018113023.txt"
id1 = len(location) + 1
id2 = len(location) + 1 + 10

# .. Parameters
if location == "Escudero":
    lat = -62.2016
    lon = -58.9657
    alt_surf = 33.0 / 1000  # m
elif location == "Sejong":
    lat = -62.225000
    lon = -58.789000
    alt_surf = 12.0 / 1000  # m
else:
    raise ValueError("Bad location, must be Escudero or Sejong.")


model_extra = 3  # model number for unspecified moledule concentrations

# .. Flags. The created profiles are always saved. These flags set whether to
#    also plot and save figures and pickled results.
plot_figs = 1  # Flag specifiying whether to plot and save figures
pickle_it = 0  # Flag specifying whether to pickle results

# .. Layer boundaries (refine as needed for desired temperature differential)
zm_lo = (
    np.array(
        [
            0.0,
            0.025,
            0.05,
            0.1,
            0.15,
            0.2,
            0.25,
            0.3,
            0.35,
            0.4,
            0.45,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
        ]
    )
    + alt_surf
)
zm_hi = np.array(
    [
        1.1,
        1.25,
        1.4,
        1.6,
        1.8,
        2.0,
        2.2,
        2.4,
        2.6,
        2.8,
        3.0,
        3.1,
        3.15,
        3.42,
        3.65,
        3.925,
        4.2,
        4.475,
        4.75,
        5.0,
        5.5,
        6.0,
        6.5,
        7.0,
        7.5,
        8.0,
        8.5,
        9.0,
        9.5,
        10.0,
        11.0,
        12.0,
        13.0,
        14.0,
        15.0,
        16.0,
        17.0,
        19.0,
        21.0,
        23.0,
        25.0,
        28.0,
        30.5,
        32.5,
        34.4,
        36.4,
        38.2,
        40.0,
        41.9,
        44.0,
        46.0,
        48.5,
        60.0,
    ]
)
z = np.hstack([zm_lo, zm_hi])
#                  1.4,  1.6  , 1.8 , 2.0,   2.275, 2.55, 2.825, 3.1, 3.375,
# 40., 41.5,  42.9, 44.8, 46.9, 49.,  60. ])
# ([19, 21, 23, 25, 27.5, 29.5, 31.5, 33.5, 35, 37,
#   40.,  42.,  44.,  46.,  48.,  50.,  55.,  58., 60.] )


# # # # # # #     START MAIN CODE     # # # # # # #

# .. Get filenames
sonde_file_n_dates = get_files(sonde_dir, sonde_file, id1, id2, "%Y%m%d%H")
ctFiles = FileInfo(co2_dir, ct_file, 23, 33, "%Y-%m-%d")
eraFiles = FileInfo(era_dir, era_filename, iera_date1, iera_date2, "%Y%m")
# surfmetFiles = FileInfo(surfmet_dir, surfmet_file, 12, 20, '%Y%m%d')

# .. Get CO2 concentrations at stations
co2_drake = CO2stationData(co2_file_drake)
co2_palmer = CO2stationData(co2_file_palmer)

# .. Chop dates down to desired range
start_ord = dt.datetime.toordinal(start_date)
end_ord = dt.datetime.toordinal(end_date)
Ndates = len(sonde_file_n_dates)
ords = [dt.datetime.toordinal(sonde_file_n_dates[i][1]) for i in range(Ndates)]
ibeg = ords.index(start_ord)
iend = ords.index(end_ord)

# .. Loop over radiosounding files, set profile, and save as netcdf file
#    For debugging, remove loop and use:
#    sonde_file = sonde_file_n_dates[0][0]
#    sonde_date = sonde_file_n_dates[0][1]
for sonde_file, sonde_date in sonde_file_n_dates[ibeg : iend + 1]:

    print("Working on", sonde_file)
    sonde = Sonde(sonde_dir, sonde_file, sonde_date)  # , z)

    fake_date = dt.datetime(2016, 12, 1, 12, 0)
    ctracker = CarbonTracker(ctFiles, fake_date, lat, lon)
    co2_surf = get_co2(sonde.date, lat, co2_drake, co2_palmer)

    era = Era(eraFiles, sonde.date, lat, lon)
    era.tack_on_60km(ctracker)  # Set values at 60 km, for interping

    prof = Prof(sonde, z)  # P, T, RH from sonde
    prof.set_n_scale("co2", ctracker, co2_surf)  # CO2 from flask measurements
    prof.set("o3", era)  # Ozone from ERA-Interim
    prof.set_upper("T", era)  # Upper T from ERA-Interim
    prof.set_upper("T", ctracker)  # Upmost T from carbonTracker
    prof.set_upper_spline("P", era)  # Upper P from ERA-Interim
    prof.h2o[prof.z >= 11.5] = 4.0  # Upmost H2o = 4 ppm
    if plot_figs:
        profo = copy.deepcopy(prof)  # For plotting figures
    # prof.set_surf_T(surfmetFiles, sonde.date)   # Set surface temperatures
    prof.model_extra = model_extra  # model for other molecs

    prof.error_check(era)  # Check for Nans, z, P

    # .. Save the output as netcdf file (optional pickle file too)
    fname = out_dir + "prof" + prof.date.strftime("%Y%m%d_%H%m")
    if pickle_it:
        pickle.dump(prof, open(fname + ".pickle", "wb"))
    prof.write_to_netcdf_file(fname + ".nc")

    # .. If specified, plot the results
    if plot_figs:
        # zmet, Tmet = get_surface_T(surfmetFiles, sonde.date)

        # .. Make figures showing results
        ylim = [-5, 65]

        plt.figure(0)
        plt.clf()
        plt.subplot(221)
        plt.cla()
        plt.plot(era.T, era.z, ".-", label="ERA")
        plt.plot(ctracker.T, ctracker.z, "o-", label="Carbon Tracker")
        plt.plot(
            sonde.T, sonde.z, linewidth=4, color=[0.6, 0.6, 0.6], label="sonde"
        )
        plt.plot(prof.T, prof.z, "k.-", label="Profile")
        # plt.plot(profo.T, profo.z, 'k.-', label = 'Profile before surf T')
        # plt.plot(Tmet, zmet, 'g*' )
        # plt.legend();
        plt.ylim(ylim)
        plt.ylabel("Temperature (K)")
        plt.title(sonde_date)

        plt.subplot(222)
        plt.cla()
        plt.plot(era.P, era.z, ".-", label="ERA")
        plt.plot(ctracker.P, ctracker.zbnd, "o-", label="Carbon Tracker")
        plt.plot(
            sonde.P, sonde.z, linewidth=4, color=[0.6, 0.6, 0.6], label="sonde"
        )
        plt.plot(prof.P, prof.z, "k", label="Profile")
        plt.legend()
        plt.ylim(ylim)
        plt.ylabel("Pressure (mb)")

        plt.subplot(223)
        plt.cla()
        plt.plot(era.rh, era.z, ".-", label="ERA")
        plt.plot(
            sonde.rh,
            sonde.z,
            linewidth=4,
            color=[0.6, 0.6, 0.6],
            label="sonde",
        )
        plt.plot(prof.rh, prof.z, "k", label="Profile")
        plt.legend()
        plt.ylim(ylim)
        plt.ylabel("Relative Humidity (%)")

        plt.subplot(224)
        plt.cla()
        plt.plot(ctracker.co2, ctracker.z, ".-", label="Carbon Tracker")
        plt.plot(prof.co2, prof.z, "k", label="Profile")
        plt.legend()
        plt.ylim(ylim)
        plt.ylabel("CO2 (ppm)")

        figname = "prof" + prof.date.strftime("%Y%m%d_%H%m") + ".png"
        plt.pause(0.1)
        plt.savefig(fig_dir + figname)

        plt.figure(1)
        plt.clf()
        plt.plot(np.diff(prof.T), prof.z[:-1] + np.diff(prof.z) / 2, ".-")
        plt.title(sonde_date)
        plt.xlabel("$\Delta$Temperature (K)")

        figname = "dT" + prof.date.strftime("%Y%m%d_%H%m") + ".png"
        plt.pause(0.1)
        plt.savefig(fig_dir + figname)

        plt.figure(2)
        plt.clf()
        plt.plot(era.T, era.z, ".-", label="ERA")
        plt.plot(ctracker.T, ctracker.z, "o-", label="Carbon Tracker")
        plt.plot(
            sonde.T, sonde.z, linewidth=4, color=[0.6, 0.6, 0.6], label="sonde"
        )
        # plt.plot(profo.T, profo.z, '-', label = 'Profile before surf T')
        plt.plot(prof.T, prof.z, "k.-", label="Profile")
        # plt.plot(Tmet, zmet, 'g*', label = 'Surface data')
        plt.legend()
        plt.ylim(ylim)
        plt.xlabel("Temperature (K)")
        maxT = np.max(
            [np.max(era.T[era.z < 1.4]), np.max(sonde.T[sonde.z < 1.4])]
        )
        plt.axis(
            [
                np.floor(np.min(prof.T[prof.z < 1.4] - 2)),
                np.ceil(maxT) + 2,
                -0.2,
                1.4,
            ]
        )

        figname = "Tlow" + prof.date.strftime("%Y%m%d_%H%m") + ".png"
        plt.pause(1)
        plt.savefig(fig_dir + figname)


"""
# .. Number of levels and indices to lat, lon, and time
ilat = np.interp(lat, nc['latitude'][:], np.arange(nc['latitude'].shape[0]))
ilon = np.interp(lon, nc['longitude'][:], np.arange(nc['longitude'].shape[0]))
dd = [d-sonde.date for d in dates]
dday = [d.days + d.seconds/60/60/24 for d in dd]
it = np.interp(0, dday, np.arange(Ntimes))

# .. As arrays
ilevs = np.arange(Nlevel)
ibnds = np.arange(Nbound)
ilatL = ilat * np.ones(Nlevel)
ilonL = ilon * np.ones(Nlevel)
itimeL = it * np.ones(Nlevel)
ilatB = ilat * np.ones(Nbound)
ilonB = ilon * np.ones(Nbound)
itimeB = it * np.ones(Nbound)

# .. Interpolate to the lat/lon of McMurdo for all times
co2b = ndimage.map_coordinates(co2mat, [[itimeL],[ilevs],[ilatL],[ilonL]])[0]
Tb = ndimage.map_coordinates(Tmat, [[itimeL],[ilevs],[ilatL],[ilonL]])[0]
Pb = ndimage.map_coordinates(Pmat, [[itimeB],[ibnds],[ilatB],[ilonB]])[0]/100
zb = ndimage.map_coordinates(zmat, [[itimeB],[ibnds],[ilatB],[ilonB]])[0]/1000
"""
