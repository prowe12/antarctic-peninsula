#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 14:40:05 2018

@author: prowe

Copyright 2019-2023 by Penny M. Rowe and NorthWest Research Associates.
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

Output: A netcdf file with:
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

# Dependencies
import datetime as dt
import bisect
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists

# import copy


# My modules
from antarc.escudero.case202202.get_atm_profs_rsrc import (
    get_files,
    FileInfo,
    Sonde,
    Era,
    Prof,
    get_co2,
    CO2stationData,
    CarbonTracker,
)


def get_atm_profs(atm_prof_params, esc_params, esc_case, redo: bool = True):
    """
    Get atmospheric profiles from sonde data, as well as ERA5, carbon tracker
    and surface co2 measurements, using imported parameters
    """
    start_date = esc_case.DATE1
    end_date = esc_case.DATE2

    lat = esc_params.LATITUDE
    lon = esc_params.LONGITUDE
    # alt_surf = esc_params.ALTITUDE
    # location = esc_params.LOCATION
    layerbnds_in = atm_prof_params.LAYERBNDS
    model_extra = atm_prof_params.MODEL_EXTRA
    sonde_dir = atm_prof_params.SONDE_DIR
    sonde_filefmt = atm_prof_params.SONDE_FILEFORMAT
    era_dir = atm_prof_params.ERA_DIR
    era_fileformat = atm_prof_params.ERA_FILEFORMAT
    co2_dir = atm_prof_params.CO2_DIR
    co2_file_drake = atm_prof_params.CO2_FILE_PALMER
    co2_file_palmer = atm_prof_params.CO2_FILE_DRAKE
    ct_filefmt = atm_prof_params.CT_FILEFORMAT
    out_dir = atm_prof_params.OUT_DIR
    fig_dir = atm_prof_params.FIG_DIR

    # Flags. The created profiles are always saved. These flags set whether to
    # also plot and save figures and pickled results.
    plotfigs = True  # Flag specifiying whether to plot and save figures

    # # # # # # #     START MAIN CODE     # # # # # # #

    # Pad out dates by a few days
    start_date = dt.datetime(
        start_date.year,
        start_date.month,
        start_date.day - 2,
        tzinfo=dt.timezone.utc,
    )
    end_date = dt.datetime(
        end_date.year,
        end_date.month,
        end_date.day + 1,
        tzinfo=dt.timezone.utc,
    )

    # Get filenames
    sonde_file_n_dates = get_files(sonde_dir, sonde_filefmt)
    eraFiles = FileInfo(era_dir, era_fileformat)

    ctFiles = FileInfo(co2_dir, ct_filefmt)  # , 23, 33, "%Y-%m-%d")
    # surfmetFiles = FileInfo(surfmet_dir, surfmet_file, 12, 20, '%Y%m%d')

    # Get CO2 concentrations at stations
    co2_drake = CO2stationData(co2_file_drake)
    co2_palmer = CO2stationData(co2_file_palmer)

    # Chop dates down to desired range
    start_ord = dt.datetime.toordinal(start_date)
    end_ord = dt.datetime.toordinal(end_date)
    ndates = len(sonde_file_n_dates)
    ords = [
        dt.datetime.toordinal(sonde_file_n_dates[i][1]) for i in range(ndates)
    ]

    # QC date range of sondes
    if start_ord < ords[0] or end_ord > ords[-1]:
        raise ValueError("sondes missing for the desired date range")

    # Get indices for date range
    ibeg = max(bisect.bisect_right(ords, start_ord) - 1, 0)
    iend = bisect.bisect_left(ords, end_ord)  # ords.index(end_ord)

    # Loop over radiosounding files, set profile, and save as netcdf file
    for sonde_file, sonde_date in sonde_file_n_dates[ibeg : iend + 1]:
        # If the file exists, do not recreate it unless redo set to False
        if not redo and exists(sonde_dir + sonde_file):
            continue

        print("Working on", sonde_file)
        sonde = Sonde(sonde_dir, sonde_file, sonde_date)
        sonde.quality_control()

        # Use the lowest sonde height as the surface height
        layerbnds = np.array(layerbnds_in)
        layerbnds[layerbnds < 1] += sonde.z[0]  # km

        # If the CT data ends a while before the desired date, try one year before
        date_for_ct = sonde_date
        lastdate = ctFiles.dates[-1]
        if date_for_ct < ctFiles.dates[0]:
            raise ValueError("No useable CarbonTracker data found.")

        while date_for_ct > lastdate and (date_for_ct - lastdate).days > 30:
            if date_for_ct < ctFiles.dates[0]:
                raise ValueError("No useable CarbonTracker data found.")
            date_for_ct = dt.datetime(
                date_for_ct.year - 1,
                date_for_ct.month,
                date_for_ct.day,
                date_for_ct.hour,
            )
        # Get the nearest date
        # ddate = [x - date_for_ct for x in ctFiles.dates]
        # idate = np.argmin(np.abs(ddate))

        ctracker = CarbonTracker(ctFiles, date_for_ct, lat, lon)
        co2_surf = get_co2(sonde.date, lat, co2_drake, co2_palmer)

        era = Era(eraFiles, era_fileformat, sonde.date, lat, lon)
        era.tack_on_60km(ctracker)  # Set values at 60 km, for interping

        prof = Prof(sonde, layerbnds)  # P, T, RH from sonde
        prof.set_n_scale("co2", ctracker, co2_surf)  # CO2 from flask measmnts
        prof.set("o3", era)  # Ozone from ERA
        prof.set_upper("T", era)  # Upper T from ERA
        prof.set_upper("T", ctracker)  # Upmost T from carbonTracker
        prof.set_upper_spline("P", era)  # Upper P from ERA-Interim
        prof.h2o[prof.z >= 11.5] = 4.0  # Upmost H2o = 4 ppm
        # if plotfigs:
        #    profo = copy.deepcopy(prof)  # For plotting figures
        # prof.set_surf_T(surfmetFiles, sonde.date)  # Set surface temperatures
        prof.model_extra = model_extra  # model for other molecs

        prof.error_check(era)  # Check for Nans, z, P

        # Save the output as netcdf file
        fname = out_dir + "prof" + prof.date.strftime("%Y%m%d_%H%m")
        prof.write_to_netcdf_file(fname + ".nc")

        # If specified, plot the results
        if plotfigs:
            # zmet, Tmet = get_surface_T(surfmetFiles, sonde.date)

            # .. Make figures showing results
            ylim = [-5, 65]

            plt.figure(0)
            plt.clf()
            plt.subplot(221)
            plt.cla()
            plt.plot(era.T, era.z, ".-", label="ERA")
            # plt.plot(ctracker.T, ctracker.z, "o-", label="Carbon Tracker")
            plt.plot(
                sonde.T,
                sonde.z,
                linewidth=4,
                color=[0.6, 0.6, 0.6],
                label="sonde",
            )
            plt.plot(prof.T, prof.z, "k.-", label="Profile")
            # plt.plot(profo.T, profo.z, 'k.-', label = 'Profile before srf T')
            # plt.plot(Tmet, zmet, 'g*' )
            # plt.legend();
            plt.ylim(ylim)
            plt.ylabel("Temperature (K)")
            plt.title(sonde_date)

            plt.subplot(222)
            plt.cla()
            plt.plot(era.P, era.z, ".-", label="ERA")
            # plt.plot(ctracker.P, ctracker.zbnd, "o-", label="Carbon Tracker")
            plt.plot(
                sonde.P,
                sonde.z,
                linewidth=4,
                color=[0.6, 0.6, 0.6],
                label="sonde",
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
            # plt.plot(ctracker.co2, ctracker.z, ".-", label="Carbon Tracker")
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
            # plt.plot(ctracker.T, ctracker.z, "o-", label="Carbon Tracker")
            plt.plot(
                sonde.T,
                sonde.z,
                linewidth=4,
                color=[0.6, 0.6, 0.6],
                label="sonde",
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
            xmin = np.floor(np.min(prof.T[prof.z < 1.4] - 2))
            plt.axis([xmin, np.ceil(maxT) + 2, -0.2, 1.4])

            figname = "Tlow" + prof.date.strftime("%Y%m%d_%H%m") + ".png"
            plt.pause(1)
            plt.savefig(fig_dir + figname)
