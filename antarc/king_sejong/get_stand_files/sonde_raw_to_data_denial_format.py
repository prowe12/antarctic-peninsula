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
import numpy as np
from os.path import exists
import datetime as dt


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


def get_datefile(direc: str, infile: str) -> str:
    """
    Return the file to check for the date in from profile file
    @param infile  The input file name
    @returns  The date string
    """
    return infile[: infile.find("_UPP")] + "_UpperAirWind Record.txt"


def get_date(direc: str, infile: str):
    """
    Return the datetime from the wind file for a given profile file
    @param infile  The input file name
    @returns  The datetime
    """
    datefile = get_datefile(direc, infile)

    if not exists(direc + datefile):
        raise ValueError(f"File {direc+infile} does not exist")

    datestr = getdatestr_from_file(direc + datefile)
    return dt.datetime.strptime(datestr, "%m/%d/%Y %H:%M:%S")


def getdatestr_from_file(datefile: str):
    """
    Return the datestring from the wind file for a given profile file
    @param infile  The input file name
    @returns  The date string
    """
    date_line = "Release time: "

    with open(datefile, encoding="ISO-8859-1") as fin:
        lines = fin.readlines()

    # Until we get to the row with the release date and time, skip
    for line in lines:
        if date_line in line:
            datestr = line[line.find(date_line) + len(date_line) :]
            datestr = datestr[: datestr.find("\n")]
            return datestr

    raise ValueError(f"Line containing release time, {date_line}, not found")


def graw_raw_to_datadenial(
    input_dir, output_dir, sample_file, location, latitude, longitude, height
):
    """
    Convert the GRAW-format radiosounding file into the format for the
    YOPP-SH data denial experiements
    @param input_dir  Directory with GRAW-format files
    @param output_dir  Directory with data-denial format files
    @param sample_file  A sample file name
    @param location  The location of the site ("Escudero") for prepending to
                     the output file name
    @param latitude  Site latitude
    @param lognitude  Site longitude
    @param height  The surface altitude
    Notes:
    Inputfile header and first line (whitespace between tokens removed for clarity)
    Profile Data:
    Time      P         T           Hu    Ws     Wd   Geopot     Dewp.
    [sec]   [hPa]     [°∆C]        [%]   [m/s]  [°∆]   [m]       [°∆C]
      0     985.3      4.6         90    6.5    335    24        3.1


    Ouput file header style:
    KingSejong
    latitude -62.201600 longitude -58.965700 height 33m
    Balloon release date and time                   2022-02-05T12:05:03
       TimeUTC          P  HeightMSL       Temp         RH       Dewp        Dir     Speed
      hh:mm:ss        hPa          m       DegC          %       DegC    Degrees     Knots
      12:05:03      996.8         28        2.1         94        1.2        304       4.0
    """

    # GRAW-format header values: when the line containing the launch date and
    # time is split, it better match the following:
    # varfordate = "Launch Date"  # Variable indicating date is in the line
    # dateline = [
    #     "Launch",
    #     "Date:",
    #     "Saturday,",
    #     "05",
    #     "February",
    #     "2022",
    #     "Launch",
    #     "Time:",
    #     "12:05:03",
    #     "End",
    #     "of",
    #     "Ascent:",
    #     "00:47:20",
    # ]
    # # These columns of dateline will be required to match those in each file
    # idateline = [0, 1, 6, 7, 9, 10, 11]
    # iday = 3
    # imonth = 4
    # iyear = 5
    # itime = 8

    # For the column labels, the whitespace seems to vary for unknown reasons
    # therefore we will remove the whitespace before comparison
    # And similarly for units, since sometimes the degree sign is rendered as
    # an infinity sign:
    # Time      P         T           Hu    Ws     Wd   Geopot     Dewp.
    # [sec]   [hPa]     [°∆C]        [%]   [m/s]  [°∆]   [m]       [°∆C]
    # cols = [
    #     "Time",
    #     "P",
    #     "T",
    #     "Hu",
    #     "Ws",
    #     "Wd",
    #     "Geopot",
    #     "Dewp.",
    # ]

    cols = {
        "Time": 0,
        "P": 1,
        "T": 2,
        "Hu": 3,
        "Ws": 4,
        "Wd": 5,
        "Geopot": 6,
        "Dewp.": 7,
    }

    units = [
        "[sec]",
        "[hPa]",
        "[°C]",
        "[%]",
        "[m/s]",
        "[°]",
        "[m]",
        "[°C]",
    ]

    # Units that are always printed correctly (so far)
    units_printed_right = ["[sec]", "[hPa]", "[%]", "[m/s]", "[m]"]
    iunit = [i for i in range(len(units)) if units[i] in units_printed_right]

    # Thresholds for quality control
    verbose = False
    alt_max = 30000
    tempmin = -100
    tempmax = 100
    dewptmin = -273
    wspd_max = 999
    maxheaderlines = 200

    # # # # # # # # # # # # #     MAIN CODE    # # # # # # # # # # # # # # # #

    # Input files
    input_files = get_inputfiles(input_dir, sample_file)

    # Log file
    logfile = "log.txt"
    lid = open(output_dir + logfile, "w")
    firstdatestr = "notset"
    datestr = "notset"
    firsttime = True
    for input_file in input_files:

        dtime = get_date(input_dir, input_file)

        start_time_s = dtime.hour * 60 * 60 + dtime.minute * 60 + dtime.second
        datestr = dtime.strftime("%Y%m%d")
        hourstr = dtime.strftime("%H")
        timestr = dtime.strftime("%H:%M:%S")
        # hourstr = f"{dtime.hour:02}"
        # timestr = hourstr + ":" + f"{dtime.minute:02}:" + f"{dtime.second:02}"

        # Read in the file data
        with open(input_dir + input_file, encoding="ISO-8859-1") as fin:
            lines = fin.readlines()

        genmsg = f"{input_file}"
        i = 0

        if firsttime:
            firstdatestr = datestr
            firsttime = False

        # if not dateset:
        #     logmsg = genmsg + ": Date not set, moving to next file"
        #     print(logmsg, file=lid)
        #     if verbose:
        #         print(logmsg)
        #     continue

        # Until we get to the row with Profile Data, we are in the header
        inheader = True
        while inheader:
            if i > maxheaderlines:
                raise ValueError("Exceeded maximum number of header lines")

            line = lines[i]
            if "Profile Data" in line:
                inheader = False
            i += 1

        # Skip two lines, which are the column labels and the units.
        # First, make sure they are as expected
        if lines[i].split() != list(cols.keys()):
            logmsg = (
                genmsg + f": Columns differ: {lines[i].split()}, moving on"
            )
            print(logmsg, file=lid)
            continue
        # Pass the line with "Profile Data"
        i += 1
        if any([lines[i].split()[j] != units[j] for j in iunit]):
            msg = (
                "On file "
                + input_file
                + ":\n"
                + lines[i]
                + "\n"
                + "The units changed; check file."
            )
            raise ValueError(msg)
        # Pass the line with units
        i += 1

        output_file = location + "_" + datestr + hourstr + ".txt"

        logmsg = genmsg + f"; Date: {datestr} {timestr}; Output: {output_file}"
        print(logmsg, file=lid)
        if verbose:
            print(logmsg)

        # WRITE THE OUTPUT FILE
        press_prev = 2000
        with open(output_dir + output_file, "w+") as fid:

            # .. Write the header
            fid.write(location + "\n")
            fid.write(
                "latitude {0:10.6f} longitude {1:10.6f} height {2:2.0f}m\n".format(
                    latitude, longitude, height
                )
            )
            fid.write("Balloon release date and time                   ")
            fid.write(datestr[:4] + "-" + datestr[4:6] + "-" + datestr[6:8])
            fid.write("T" + timestr + "\n")
            fid.write(
                "   TimeUTC          P  HeightMSL       Temp         RH       "
            )
            fid.write("Dewp        Dir     Speed\n")
            fid.write(
                "  hh:mm:ss        hPa          m       DegC          %       "
            )
            fid.write("DegC    Degrees       m/s\n")

            # Go through the rest of the record until we get to 'Trop' or a bad line
            foundbadalt = False
            for iline, line in enumerate(lines[i:]):

                if len(line) < 4:
                    raise ValueError("The line is too short!")
                if line[:4] == "Trop":
                    # Reached end of sounding data, done
                    break

                # Get the columns
                columns = line.split()
                # time_s = start_time_s + float(columns[0])
                # press = float(columns[1])
                # temp = float(columns[2])
                # rhw = float(columns[3])
                # wspd = float(columns[4])
                # wdir = float(columns[5])
                # alt = float(columns[8])
                # dewpt = float(columns[10])

                time_s = start_time_s + float(columns[cols["Time"]])
                press = float(columns[cols["P"]])
                temp = float(columns[cols["T"]])
                rhw = float(columns[cols["Hu"]])
                wspd = float(columns[cols["Ws"]])
                wdir = float(columns[cols["Wd"]])
                alt = float(columns[cols["Geopot"]])
                dewpt = float(columns[cols["Dewp."]])

                if press > press_prev and iline > 2:
                    # The pressure is increasing - balloon dropping - thus done
                    logmsg = (
                        genmsg
                        + f": Pressure increased from {press_prev} to {press} on line {iline}, assuming sounding over"
                    )
                    print(logmsg, file=lid)
                    break

                # Quality control: Check for values outside thresholds
                if press < 0 or press > 1100:
                    raise ValueError(
                        genmsg + "Pressure falls outside allowed threshold"
                    )
                if temp < tempmin or temp > tempmax:
                    # Skip this line
                    print(genmsg + ": Temperature outside bounds, skipping")
                    continue
                if rhw < 0 or rhw > 100:
                    msg = f": RH ({rhw}) outside 0 - 100 at {alt} m, using nan"
                    print(genmsg + msg, file=lid)
                    rhw = np.nan
                if dewpt < dewptmin or dewpt > temp:
                    msg = f": Dew point ({dewpt}) outside -100 to {temp}) at {alt} m, using nan"
                    print(genmsg + msg, file=lid)
                    dewpt = np.nan

                if wdir < 0 or wdir > 360 or wspd < 0 or wspd > wspd_max:
                    raise ValueError(
                        genmsg + "Wind data falls outside allowed thresholds"
                    )

                if alt < 0 or alt > alt_max:
                    if not foundbadalt:
                        # First time finding bad alt, so print message to log
                        # (next time we won't, to avoid numerous messages)
                        logmsg = (
                            genmsg
                            + f": One or more altitude (e.g. {alt} m) outside 0 - {alt_max} m, using nan;"
                        )
                        print(logmsg, file=lid)
                        foundbadalt = True
                    alt = np.nan

                # Time as HH:MM:ss
                hour = int(np.floor(time_s / 3600))  # s x min/s x hr/min
                mint = int(np.floor(time_s / 60) - hour * 60)  # minutes
                sec = int(time_s - hour * 3600 - mint * 60)  # seconds
                hhmmss = "  {:02}:{:02}:{:02}".format(hour, mint, sec)

                p_z_t_rh = "{0:11.1f} {1:10.0f} {2:10.1f} {3:10.0f}".format(
                    press, alt, temp, rhw
                )
                dewpt_wdir_spd = "{0:11.1f} {1:10.0f} {2:9.1f}".format(
                    dewpt, wdir, wspd
                )
                print(hhmmss + p_z_t_rh + dewpt_wdir_spd, file=fid)

                press_prev = press
    lid.close()
    # Rename log file according to date range
    os.rename(
        output_dir + logfile, output_dir + f"log_{firstdatestr}_{datestr}.txt"
    )
