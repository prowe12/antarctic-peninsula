#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:38:16 2019

Copyright 2022 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

Purpose: Convert output in the standard format output by GRAW after the
simulation (UPP RAW) into the format desired for the data denial experiment

@author: prowe
"""

import os
from time import strptime, strftime
from os.path import exists
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import datetime as dt


# My modules
from antarc.escudero.all.get_atm_profs_rsrc import get_files, Sonde


def get_inputfiles(direc, sample_file):
    """
    Return a list of all files of the form of sample_file
    (currently same length and same last 3 characters
    @param direc  The input directory
    @param sample_file  An example file name
    @return A list of filenames (strings)
    """
    all_input_files = np.sort(os.listdir(direc))
    input_files = []
    for file in all_input_files:
        if len(file) == len(sample_file) and file[-3:] == sample_file[-3:]:
            input_files.append(file)
    return input_files


def qcdatelist(datelist):
    """
    Quality control the datelist
    @param datelist: day, month, year, time as hh:mm:ss
    @throws exception if format of any is incorrect
    """
    months = [
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
    ]
    year = int(datelist[2])
    month = datelist[1]
    day = int(datelist[0])
    if year < 2014 or year > 2050:
        raise ValueError("Year falls outside allowable range!")
    if month not in months:
        raise ValueError("Month not in correct format")
    if day < 0 or day > 31:
        raise ValueError("Day of month not witin 1-31")


def getdatetuplefrom_datelist(datelist):
    """
    Return a date tuple from a list with dd, Month, yyyy, hh:mm:ss
    @param list with dd, Month, yyyy, hh:mm:ss
    @return date tuple
    """
    dd_month_yyyy = "-".join(datelist[:3])
    dd_month_yyyy_time = "".join([dd_month_yyyy, " ", datelist[3]])
    return strptime(dd_month_yyyy_time, "%d-%B-%Y %H:%M:%S")


def write_header(
    fid,
    location: str,
    lat: float,
    lon: float,
    height: float,
    datestr: str,
    timestr: str,
):
    """
    Create the header for the output file in 'data denial' format
    @parm noqc  The id to the file to write to
    @parm location  The radiosonde launch location
    @parm lat  The latitude
    @parm lon  The longitude
    @parm height  The altitude
    @parm datestr  String with the date
    @parm timestr  String with the time
    """
    str_latlon = "latitude {0:10.6f} longitude {1:10.6f} height {2:2.0f}m\n"
    fid.write(location + "\n")
    fid.write(str_latlon.format(lat, lon, height))
    fid.write("Balloon release date and time                   ")
    fid.write(datestr[:4] + "-" + datestr[4:6] + "-" + datestr[6:8])
    fid.write("T" + timestr + "\n")
    fid.write("   TimeUTC          P  HeightMSL       Temp         RH       ")
    fid.write("Dewp        Dir     Speed\n")
    fid.write("  hh:mm:ss        hPa          m       DegC          %       ")
    fid.write("DegC    Degrees     Knots\n")


def qc_burst(met: np.ndarray, ialt: int, ipress: int):
    """
    Chop the profile data at the balloon burst height
    @param met  Array of profile data
    @param ialt  Index to column with altitudes
    @param ipress  Index to column with pressures
    @returns  The profile data up to the burst height
    """
    # Get indices to the location in the array of the maximum altitude
    # and the minimum pressure.  Ideally these would both be at the burst.
    imaxalt = np.where(met[:, ialt] == max(met[:, ialt]))[0]
    iminp = np.where(met[:, ipress] == min(met[:, ipress]))[0]

    # Are they length 1?
    if len(imaxalt) == 1 and len(iminp) == 1:
        # If they are length 1, are they the same?
        if imaxalt[0] == iminp[0] and imaxalt[0] == len(met) - 1:
            return met

        if abs(imaxalt - iminp)[0] < 6:
            # If they are length 1, and close, end sounding at lowest
            met = met[: min(imaxalt, iminp)[0] + 1, :]
            return met

        # Otherwise return the highest
        met = met[: max(imaxalt, iminp)[0] + 1, :]
        return met

        raise ValueError("Bad value for met data")

    # If there are more than one found, are they all near each other?
    ilowest = min(imaxalt.min(), iminp.min()) + 1
    ihighest = max(imaxalt.max(), iminp.max(0)) + 1
    if ihighest - ilowest < 16:
        # If they are all near each other, use the lowest height
        # where pressure is the minimum or altitude the maximum
        met = met[:ilowest, :]
        return met

    if np.all(np.diff(met[ilowest:, ipress]) >= 0) and np.all(
        np.diff(met[ilowest:, ialt]) <= 0
    ):
        # If all the pressures are increasing or staying same and
        # all the altitudes are dropping or staying the same
        # starting from the lowest max alt and min pressure
        # Then the balloon is falling and we can crop it
        met = met[:ilowest, :]
        return met

    met = met[: iminp[0] + 1, :]
    return met


def compare_profiles(
    dir1: str, dir2: str, fmt: str, lab1: str, lab2: str, savedir: str
):
    """
    Create figures that compare radiosonde profiles between radiosonde results
    in data denial format in two directories, e.g. where one has some QC and
    The other does not
    @param sampfile  Name of a sample file to match
    @param dir1  The first directory
    @param dir2  The second directory
    @param fmt  The format of the files to be loaded in
    @param lab1  Label for curves corresponding to dir1
    @param lab2  Label for curves corresponding to dir2
    @param savedir  Directory to save created figures
    """

    # Get filenames
    sonde_file_n_dates = get_files(dir1, fmt)

    for sonde_file, sonde_date in sonde_file_n_dates:
        if not exists(dir2 + sonde_file):
            continue
        print("Working on", sonde_file)
        sonde1 = Sonde(dir1, sonde_file, sonde_date)
        sonde2 = Sonde(dir2, sonde_file, sonde_date)

        zmax = min(max(np.nanmax(sonde1.z), np.nanmax(sonde2.z)), 30) + 1
        h2omax = (
            min(max(np.nanmax(sonde1.h2o), np.nanmax(sonde2.h2o)), 10000) + 100
        )
        rhmax = min(max(np.nanmax(sonde1.rh), np.nanmax(sonde2.rh)), 120) + 2
        tmax = min(max(np.nanmax(sonde1.T), np.nanmax(sonde2.T)), 300) + 5
        tmin = max(min(np.nanmin(sonde1.T), np.nanmin(sonde2.T)), 190) - 10

        plt.figure(1)
        plt.clf()
        plt.subplot(221)
        plt.plot(sonde2.P, sonde2.z, ".-", label=lab2)
        plt.plot(sonde1.P, sonde1.z, label=lab1)
        plt.title(sonde_file)
        plt.xlabel("Pressure")
        plt.ylabel("Altitude (km)")
        plt.legend()
        plt.ylim([-0.4, zmax])

        plt.subplot(222)
        plt.plot(sonde2.T, sonde2.z, ".-")
        plt.plot(sonde1.T, sonde1.z)
        plt.xlabel("Temperature")
        plt.ylabel("Altitude (km)")
        plt.ylim([-0.4, zmax])
        plt.xlim([tmin, tmax])

        plt.subplot(223)
        plt.plot(sonde2.h2o, sonde2.z, ".-")
        plt.plot(sonde1.h2o, sonde1.z)
        plt.plot([-100, h2omax], [0, 0], "k:")
        plt.plot([0, 0], [-0.4, zmax], "k:")
        plt.xlabel("H2O")
        plt.xlim([-100, h2omax])
        plt.ylim([-0.4, zmax])

        plt.subplot(224)
        plt.plot(sonde2.rh, sonde2.z, ".-")
        plt.plot(sonde1.rh, sonde1.z)
        plt.xlabel("RH (%)")
        plt.xlim([-2, rhmax])
        plt.ylim([-0.4, zmax])

        figname = sonde_file[:-3] + "png"
        plt.savefig(savedir + figname)


def qc_graw_std(
    input_dir: str,
    output_dir: str,
    sample_file: str,
    location,
    lat,
    lon,
):
    """
    Convert the GRAW-format radiosounding file into the format for the
    YOPP-SH data denial experiements
    @param input_dir  Directory with GRAW-format files
    @param output_dir  Directory with data-denial format files
    @param sample_file  A sample file name
    @param location  The location of the site ("Escudero") for prepending to
                     the output file name
    @param lat  Site latitude
    @param lon  Site longitude
    @param height  The surface altitude


    Inputfile header style (excluding unused lines and removing most white space):

    Launch Date: Saturday, 05 February 2022   Launch Time: 12:05:03 End of Ascent: 00:47:20
    ...
    Profile Data:
    Time    P    T   Hu   Ws  Wd  Long.     Lat.   Alt  Geopot Dewp. Virt. Temp Rs
    [sec] [hPa] [∞C] [%] [kn] [∞] [∞]       [∞]    [m]   [m]   [∞C]    [∞C]  [m/s]
      0   996.8  2.1  94 4.0  304 -58.965700 -62.201600  28  28  1.2  2.8  0.0
    ...
    Tropopauses:
    ...


    Ouput file header style:
    Escudero
    latitude -62.201600 longitude -58.965700 height 33m
    Balloon release date and time                   2022-02-05T12:05:03
       TimeUTC          P  HeightMSL       Temp         RH       Dewp        Dir     Speed
      hh:mm:ss        hPa          m       DegC          %       DegC    Degrees     Knots
      12:05:03      996.8         28        2.1         94        1.2        304       4.0

    Quality Control
    QC is minimal and includes:
        -  Chop off values above the burst height
        - Check if pressures monotonically decreasing or same
          and remove second instance of any same or increasing pressure
        - Check if alts monotonically increasing or same
          and remove decreasing values. Report to log.
        - Replace values that fall outside allowed boundaries with nans

    Not removed by QC:
        nan values are left where values were outside allowed limits
        altitudes may decrease with time / decreasing pressure
        (but pressures should never increases with time)
        This should be done in the method qc of the class Sonde
    """

    # # # # # # # # # # # # #     MAIN CODE    # # # # # # # # # # # # # # # #

    plot_qc_fig1 = False
    plot_qc_fig2 = False

    # Thresholds for quality control
    # verbose = False
    # time_max = 2 * 24 * 3600
    alt_max = 40000
    tempmin = -100
    tempmax = 100
    dewptmin = -273
    dewptmax = 100
    wspd_max = 999
    # metmin = [0, 0, 0, tempmin, 0, dewptmin, 0, 0]
    # metmax = [time_max, 1100, alt_max, tempmax, 150, dewptmax, 360, wspd_max]
    # Exluding time!
    #     press, alt, temp, rhw, dewpt, wdir, wspd]).T
    metmin = [0, 0, tempmin, 0, dewptmin, 0, 0]
    metmax = [1100, alt_max, tempmax, 150, dewptmax, 360, wspd_max]

    # maxheaderlines = 200
    pztrhstr = "{0:11.1f} {1:10.0f} {2:10.1f} {3:10.0f}"
    dewstr = "{0:11.1f} {1:10.0f} {2:9.1f}"

    # Input files
    input_files = get_inputfiles(input_dir, sample_file)

    # Log file
    logfile = "log.txt"
    lid = open(output_dir + logfile, "w", encoding="utf8")
    firstdatestr = "notset"

    for input_file in input_files:
        # input_file = input_files[0]

        # Radiosonde Data
        df = pd.read_csv(input_dir + input_file, na_values="nan")
        df.columns = df.columns.str.replace(" ", "")
        for col in df.columns:
            if col != "Date":
                df[col] = df[col].astype(float)

        temp = df["Temp(C)"].to_numpy().astype(float)
        press = df["P(hPa)"].to_numpy().astype(float)
        alt = df["Alt(m)"].to_numpy().astype(float)
        wdir = df["wdir(deg)"].to_numpy().astype(float)
        wspd = df["wspd(kn)"].to_numpy().astype(float)
        rhw = df["RH(%)"].to_numpy().astype(float)
        dewpt = df["Dewp(C)"].to_numpy().astype(float)

        try:
            dtime = np.array(
                [
                    dt.datetime.strptime(x, "%Y-%m-%d_%H:%M:%S")
                    for x in df["Date"]
                ]
            )

        except:
            print("pause here")

        genmsg = f"{input_file}"

        # output_file = location + "_" + datestr + hourstr + ".txt"
        # start_time_s = launchhour * 3600 + launchmin * 60 + launchsec

        # logmsg = genmsg + f"; Date: {datestr} {timestr}; Output: {output_file}"
        # print(logmsg, file=lid)
        # if verbose:
        #     print(logmsg)

        # nlevs = len(lines[i:])
        # # time_s, press, alt, temp, rhw, dewpt, wdir, wspd
        # met = np.nan * np.ones((nlevs, 8))

        met = np.vstack([dtime, press, alt, temp, rhw, dewpt, wdir, wspd]).T
        # # time_s, press, alt, temp, rhw, dewpt, wdir, wspd
        # met = [
        #     start_time_s + float(columns[colmap["time"][0]]),
        #     float(columns[colmap["press"][0]]),
        #     float(columns[colmap["alt"][0]]),
        #     float(columns[colmap["temp"][0]]),
        #     float(columns[colmap["rhw"][0]]),
        #     float(columns[colmap["dewpt"][0]]),
        #     float(columns[colmap["wdir"][0]]),
        #     wspd,
        # ]

        datestr = dtime[0].strftime("%Y%m%d")

        if firstdatestr == "notset":
            firstdatestr = datestr

        # hourstr = dtime[0].strftime("%H")
        timestr = dtime[0].strftime("%H:%M:%S")
        # timestr = hourstr + ":" + f"{launchmin:02}" + ":" + f"{launchsec:02}"
        # output_file = location + "_" + datestr + hourstr + ".txt"
        output_fig = dtime[0].strftime("esc_sonde_dd_v2_%Y%m%d%H")
        output_file = dtime[0].strftime("esc_sonde_dd_v2_%Y%m%d%H.txt")

        # Remove extra values
        itm = 0
        ipress = 1
        ialt = 2
        itemp = 3

        # # Write an output file with no qc
        # noqcfile = output_dir + "noQC/" + output_file
        # with open(noqcfile, "w+", encoding="utf8") as noqc:
        #     # Write the header
        #     write_header(noqc, location, lat, lon, height, datestr, timestr)

        #     for k in range(j):
        #         # Time as HH:MM:ss
        #         hour = int(np.floor(met[k, itm] / 3600))  # s x min/s x hr/min
        #         mint = int(np.floor(met[k, itm] / 60) - hour * 60)  # minutes
        #         sec = int(met[k, itm] - hour * 3600 - mint * 60)  # seconds
        #         hhmmss = f"  {hour:02d}:{mint:02d}:{sec:02d}"

        #         # p_z_t_rh = pztrhstr.format(press[k], alt[k], temp[k], rhw[k])
        #         # dewpt_wdir_spd = dewstr.format(dewpt[k], wdir[k], wspd[k])
        #         p_z_t_rh = pztrhstr.format(
        #             met[k, 1], met[k, 2], met[k, 3], met[k, 4]
        #         )
        #         dewpt_wdir_spd = dewstr.format(met[k, 5], met[k, 6], met[k, 7])
        #         print(hhmmss + p_z_t_rh + dewpt_wdir_spd, file=noqc)

        # Uncomment to create a plot for QC
        # Get the minimum pressure and the maximum altitude
        # iminp = np.where(met[:, ipress] == min(met[:, ipress]))[0]
        # imaxalt = np.where(met[:, ialt] == max(met[:, ialt]))[0]
        # plt.figure(1)
        # plt.plot(met[:, ipress], met[:, ialt])
        # plt.plot(met[imaxalt, ipress], met[imaxalt, ialt], "*")
        # plt.plot(met[iminp, ipress], met[iminp, ialt], "*")
        # plt.xlabel("Pressure")
        # plt.ylabel("Altitude")

        # print(output_file)

        if plot_qc_fig1:
            # Replace bad values with large number
            for i in range(1, np.shape(met)[1]):
                met[met[:, i] > 100000, i] = 50000

            plt.figure(1)
            plt.clf()
            plt.subplot(121)
            plt.plot(met[:, ipress], met[:, ialt])
            plt.subplot(122)
            plt.plot(met[:, itemp], met[:, ialt])

        # Replace bad values with nan
        for i in range(1, np.shape(met)[1]):
            met[met[:, i] >= 50000, i] = np.nan

        # Chop off values above the burst height
        met = qc_burst(met, ialt, ipress)

        if plot_qc_fig1:
            plt.figure(1)
            plt.subplot(121)
            plt.plot(met[:, ipress], met[:, ialt])
            plt.subplot(122)
            plt.plot(met[:, itemp], met[:, ialt])
            xmin, xmax, ymin, ymax = plt.axis()
            plt.xlim([round(min(met[:, itemp])) - 5, 10])

            plt.savefig(output_dir + output_fig)

        if plot_qc_fig2:
            plt.figure(2)
            plt.clf()
            plt.subplot(121)
            plt.plot(met[:, ipress], met[:, ialt])
            plt.subplot(122)
            plt.plot(met[:, itemp], met[:, ialt])
            xmin, xmax, ymin, ymax = plt.axis()

        # Remove any layers where both pressure and altitude are nan
        # (Note: I have not found any)
        iznan = np.where(np.isnan(met[:, ialt].astype(float)))[0]
        ipnan = np.where(np.isnan(met[:, ipress].astype(float)))[0]
        inan = np.intersect1d(iznan, ipnan)
        if any(inan):
            print("pause here")

        # Remove any levels with nan for pressure or alt
        inan = np.union1d(iznan, ipnan)
        met = np.delete(met, inan, axis=0)

        # Check if pressures monotonically decreasing or same
        if np.any(np.diff(met[:, ipress]) >= 0):
            # Remove second instance of any same or increasing pressure
            press = met[:, ipress]
            while np.any(np.diff(press) >= 0):
                idup = np.where(np.diff(press) >= 0)[0]
                met = np.delete(met, (idup + 1), axis=0)
                press = met[:, ipress]

        # Check if alts monotonically increasing or same
        if np.any(np.diff(met[:, ialt]) <= 0):
            # Remove second instance of any same or decreasing alt
            alt = met[:, ialt]
            while np.any(np.diff(alt) <= 0):
                idup = np.where(np.diff(alt) <= 0)[0]
                met = np.delete(met, (idup + 1), axis=0)
                alt = met[:, ialt]

        # End at burst height (again)
        if np.any(np.diff(met[:, ialt]) <= 0):
            met = qc_burst(met, ialt, ipress)

        if np.any(np.diff(met[:, ialt]) <= 0):
            msg = "; Warning: an altitude decreased with time"
            print(genmsg + msg, file=lid)

        # Check if altitudes monotonically increasing
        if not np.all(np.diff(met[:, ialt]) > 0):
            raise ValueError("nope")

        # # Uncomment to add to QC plot
        # plt.figure(1)
        # plt.plot(met[:, ipress], met[:, ialt])

        # if np.any(met > metmax) or np.any(met < metmin):
        #     for k in range(np.shape(met)[1]):
        #         if np.any(met[:, k] > metmax[k]):
        #             met[met[:, k] > metmax[k], k] = np.nan
        #         if np.any(met[:, k] < metmin[k]):
        #             met[met[:, k] < metmin[k], k] = np.nan

        # No pressure should increase with time
        if np.any(np.diff(met[:, ipress]) >= 0):
            raise ValueError("Pressures should not increase with time!")

        # No times should be nan
        # if np.any(np.isnan(met[:, itm].astype(float))):
        #     raise ValueError("All times must exist")

        # Quality control: Check for values outside thresholds
        if np.any(met[:, 1:] > metmax):
            print("   P   ALT    Temp  RH  dewpt wdir wspd")
            print(metmax)
            print(met[np.where(met[:, 1:] > metmax)[0], 1:])
            raise ValueError("")

        if np.any(met[:, 1:] < metmin):
            raise ValueError("")

        # Bad data - do not save
        if input_file == "esc_sonde_v0_2019100317.txt":
            continue

        if plot_qc_fig2:
            plt.figure(2)
            plt.subplot(121)
            plt.plot(met[:, ipress], met[:, ialt])
            plt.subplot(122)
            plt.plot(met[:, itemp], met[:, ialt])
            xmin, xmax, ymin, ymax = plt.axis()
            plt.xlim([round(min(met[:, itemp])) - 5, 10])

            plt.savefig(output_dir + output_fig)

        # Time as HH:MM:ss
        # hour = int(np.floor(time_s[j] / 3600))  # s x min/s x hr/min
        # mint = int(np.floor(time_s[j] / 60) - hour * 60)  # minutes
        # sec = int(time_s[j] - hour * 3600 - mint * 60)  # seconds
        # hhmmss = "  {:02}:{:02}:{:02}".format(hour, mint, sec)

        # p_z_t_rh = "{0:11.1f} {1:10.0f} {2:10.1f} {3:10.0f}".format(
        #     press[j], alt[j], temp[j], rhw[j]
        # )
        # dewpt_wdir_spd = "{0:11.1f} {1:10.0f} {2:9.1f}".format(
        #     dewpt[j], wdir[j], wspd[j]
        # )
        # print(hhmmss + p_z_t_rh + dewpt_wdir_spd, file=noqc)

        # press_prev = press

        # Create the QCd file
        height = met[0, ialt]
        outfile = output_dir + output_file
        with open(outfile, "w+", encoding="utf8") as fid:
            write_header(fid, location, lat, lon, height, datestr, timestr)
            for k in range(len(met)):
                # Time as HH:MM:ss
                # hour = int(np.floor(met[k, itm] / 3600))  # s x min/s x hr/min
                # mint = int(np.floor(met[k, itm] / 60) - hour * 60)  # minutes
                # sec = int(met[k, itm] - hour * 3600 - mint * 60)  # seconds
                # hhmmss = f"  {hour:02d}:{mint:02d}:{sec:02d}"
                hhmmss = met[k, 0].strftime("  %H:%M:%S")

                p_z_t_rh = pztrhstr.format(
                    met[k, 1], met[k, 2], met[k, 3], met[k, 4]
                )
                dewpt_wdir_spd = dewstr.format(met[k, 5], met[k, 6], met[k, 7])
                print(hhmmss + p_z_t_rh + dewpt_wdir_spd, file=fid)

    lid.close()
    # Rename log file according to date range
    os.rename(
        output_dir + logfile, output_dir + f"log_{firstdatestr}_{datestr}.txt"
    )
