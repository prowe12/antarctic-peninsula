#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 23:01:01 2022

@author: prowe

Get the sounding data and save it to a common format, with no QC
"""

# Modules
import numpy as np
import os
from time import strptime, strftime

# My parameters and modules
import antarc.params
from antarc.escudero.parameters import esc_params
from antarc.escudero.get_stand_files.graw_get_column_map import (
    graw_get_column_map,
)


def getdatetuplefrom_datelist(datelist):
    """
    Return a date tuple from a list with dd, Month, yyyy, hh:mm:ss
    @param list with dd, Month, yyyy, hh:mm:ss
    @return date tuple
    """
    dd_month_yyyy = "-".join(datelist[:3])
    dd_month_yyyy_time = "".join([dd_month_yyyy, " ", datelist[3]])
    return strptime(dd_month_yyyy_time, "%d-%B-%Y %H:%M:%S")


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


def graw_to_std_fmt(
    input_dir: str,
    output_dir: str,
    sample_file: str,
    location,
    lat,
    lon,
    height,
):
    """
    Save GRAW output into a csv file with a standard format
    """

    # Debug params
    verbose = True

    # GRAW-format header values: when the line containing the launch date and
    # time is split, it better match the following:
    varfordate = "Launch Date"  # Variable indicating date is in the line
    dateline = [
        "Launch",
        "Date:",
        "Saturday,",
        "05",
        "February",
        "2022",
        "Launch",
        "Time:",
        "12:05:03",
        "End",
        "of",
        "Ascent:",
        "00:47:20",
    ]
    # These columns of dateline will be required to match those in each file
    idateline = [0, 1, 6, 7, 9, 10, 11]
    iday = 3
    imonth = 4
    iyear = 5
    itime = 8

    # # Thresholds for quality control
    # verbose = False
    # alt_max = 30000
    # tempmin = -100
    # tempmax = 100
    # dewptmin = -273
    # dewptmax = 100
    # wspd_max = 999
    # time_max = 2 * 24 * 3600
    # metmin = [0, 0, 0, tempmin, 0, dewptmin, 0, 0]
    # metmax = [time_max, 1100, alt_max, tempmax, 150, dewptmax, 360, wspd_max]

    maxheaderlines = 200
    pztrhstr = "{0:8.1f}, {1:8.0f}, {2:8.1f}, {3:8.0f},"
    dewstr = "{0:8.1f}, {1:8.0f}, {2:8.1f}"

    # # # # # # # # # # # # # #     MAIN CODE    # # # # # # # # # # # # # # # #

    # Input files
    input_files = get_inputfiles(input_dir, sample_file)

    # Log file
    logfile = "log.txt"
    lid = open(output_dir + logfile, "w", encoding="utf8")
    firstdatestr = "notset"
    datestr = "notset"
    firsttime = True
    for input_file in input_files:

        # Read in the file data
        with open(input_dir + input_file, encoding="ISO-8859-1") as fin:
            lines = fin.readlines()

        # Read the header until we get the date
        genmsg = f"{input_file}"
        dateset = False
        i = 0
        while not dateset:
            if i > maxheaderlines:
                msg = ": Exceeded maximum number of header lines, moving on"
                print(genmsg + msg, file=lid)
                break

            line = lines[i]
            if varfordate in line:
                linelist = line.split()
                if any((linelist[j] != dateline[j] for j in idateline)):
                    msg = f"{genmsg}{line}\nLine with date, indicated by "
                    msg2 = f"{varfordate} has unexpected format; check file."
                    print(genmsg + msg + msg2, file=lid)

                datelist = [
                    linelist[iday],
                    linelist[imonth],
                    linelist[iyear],
                    linelist[itime],
                ]
                qcdatelist(datelist)
                datetuple = getdatetuplefrom_datelist(datelist)
                datestr = strftime("%Y%m%d", datetuple)
                launchhour = datetuple.tm_hour
                launchmin = datetuple.tm_min
                launchsec = datetuple.tm_sec
                dateset = True
            i += 1

        if firsttime:
            firstdatestr = datestr
            firsttime = False

        if not dateset:
            logmsg = genmsg + ": Date not set, moving to next file"
            print(logmsg, file=lid)
            if verbose:
                print(logmsg)
            continue

        # Until we get to the row with Profile Data, we are in the header
        inheader = True
        while inheader:
            if i > maxheaderlines:
                raise ValueError("Exceeded maximum number of header lines")

            line = lines[i]
            if "Profile Data" in line:
                inheader = False
            i += 1

        # Get column labels and units
        file_cols = [x.strip() for x in lines[i].split("\t")]
        i += 1
        file_units = lines[i].split()

        # Make sure columns and units are among those expected
        colmap = graw_get_column_map(file_cols, file_units)

        # Check if wind speed will need to be converted from m/s to knots
        # For this we use index 1 of the column map, which corresponds to
        # units (0 corresponds to column values)
        if file_units[colmap["wspd"][1]] == "[m/s]":
            wspd_mps_to_knots = True
        elif file_units[colmap["wspd"][1]] == "[kn]":
            wspd_mps_to_knots = False
        else:
            print(f"Units are: {file_units[colmap['wspd'][1]]}")
            raise ValueError("Bad value for units!")

        i += 1

        hourstr = f"{launchhour:02}"
        timestr = hourstr + ":" + f"{launchmin:02}" + ":" + f"{launchsec:02}"

        output_file = location + "_" + datestr + hourstr + ".txt"
        start_time_s = launchhour * 3600 + launchmin * 60 + launchsec

        logmsg = genmsg + f"; Date: {datestr} {timestr}; Output: {output_file}"
        print(logmsg, file=lid)
        if verbose:
            print(logmsg)

        nlevs = len(lines[i:])
        # time_s, press, alt, temp, rhw, dewpt, wdir, wspd
        met = np.nan * np.ones((nlevs, 8))

        # Go through  rest of the record until we get to 'Trop' or a bad line
        # foundbadalt = False
        j = 0
        for j, line in enumerate(lines[i:]):

            if len(line) < 4:
                raise ValueError("The line is too short!")
            if line[:4] == "Trop":
                # Reached end of sounding data, done
                break

            # Get the columns
            columns = line.split()

            # Wind speed is sometimes m/s but needs to be knots
            # The boolean wspd_mps_to_knots, set above, determines if
            # the conversion is needed
            # Recall that the index 0 here gets the index to the columns
            # for the values rather than the units (index 1)
            wspd = float(columns[colmap["wspd"][0]])
            if wspd_mps_to_knots:
                wspd *= 1.9438445

            # time_s, press, alt, temp, rhw, dewpt, wdir, wspd
            met[j, :] = [
                start_time_s + float(columns[colmap["time"][0]]),
                float(columns[colmap["press"][0]]),
                float(columns[colmap["alt"][0]]),
                float(columns[colmap["temp"][0]]),
                float(columns[colmap["rhw"][0]]),
                float(columns[colmap["dewpt"][0]]),
                float(columns[colmap["wdir"][0]]),
                wspd,
            ]

        # Remove extra values
        met = met[:j, :]
        # itm = 0
        # ymd = datestr[:4] + "-" + datestr[4:6] + "-" + datestr[6:8]

        # Get datetime
        import datetime as dt

        dtime0 = dt.datetime.strptime(datestr, "%Y%m%d")
        dtime = [dtime0 + dt.timedelta(seconds=x) for x in met[:, 0]]

        # Write an output file with no qc
        noqcfile = output_dir + output_file
        with open(noqcfile, "w+", encoding="utf8") as fid:
            # Write the header
            fid.write("Date,                  P(hPa),   Alt(m),   Temp(C), ")
            fid.write("  RH(%),  Dewp(C), wdir(deg), wspd(kn)\n")

            for k, dtime0 in enumerate(dtime):
                # Time as HH:MM:ss
                # hour = int(np.floor(met[k, itm] / 3600))  # s x min/s x hr/min
                # mint = int(np.floor(met[k, itm] / 60) - hour * 60)  # minutes
                # sec = int(met[k, itm] - hour * 3600 - mint * 60)  # seconds
                # hhmmss = f"_{hour:02d}:{mint:02d}:{sec:02d},"
                dateout = dtime0.strftime("%Y-%m-%d_%H:%M:%S,")

                # p_z_t_rh = pztrhstr.format(press[k], alt[k], temp[k], rhw[k])
                # dewpt_wdir_spd = dewstr.format(dewpt[k], wdir[k], wspd[k])
                p_z_t_rh = pztrhstr.format(
                    met[k, 1], met[k, 2], met[k, 3], met[k, 4]
                )
                dewpt_wdir_spd = dewstr.format(met[k, 5], met[k, 6], met[k, 7])
                print(dateout + p_z_t_rh + dewpt_wdir_spd, file=fid)

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

    lid.close()
    # Rename log file according to date range
    os.rename(
        output_dir + logfile, output_dir + f"log_{firstdatestr}_{datestr}.txt"
    )


# Directories from parameter file
MEAS_DIR = antarc.params.MEAS_DIR

# Set variables from parameters
LAT = esc_params.LATITUDE
LON = esc_params.LONGITUDE
HEIGHT = esc_params.ALTITUDE


# Directories
SONDE_DIR = MEAS_DIR + "Escudero/GRAW_radiosondes/"
ORIG_DIRS = (
    SONDE_DIR + "simulation/eca54_escudero_profiles_simulation2019/",
    SONDE_DIR + "simulation/eca55_escudero_profiles_a/",
    SONDE_DIR + "simulation/eca55_escudero_profiles_b/",
    SONDE_DIR + "simulation/eca56_escudero_profiles/",
    SONDE_DIR + "simulation/eca58_escudero_profiles/",
    SONDE_DIR + "simulation/yopp2022_escudero_profiles/",
    SONDE_DIR + "simulation/eca59_escudero_profiles/",
    SONDE_DIR + "simulation/simulation_reruns_2023/",
)
OUTDIR = MEAS_DIR + "Escudero/radiosondes/graw/"


PREFIX = "esc_sonde_v0"


# Run the code to convert from the GRAW format to the data denial format
# sampfile = sonde_params.SAMPLEFNAME
SAMPFILE = "20220401120020047677_UPP_RAW_89056_2022040112.txt"
for indir in ORIG_DIRS:
    slash = indir[: indir.rfind("/")].rfind("/")
    graw_to_std_fmt(indir, OUTDIR, SAMPFILE, PREFIX, LAT, LON, HEIGHT)


# For first field season (2017/01), use these too
INDIR = SONDE_DIR + "soundings_txt_YYYYMMDD_hhmm/"
SAMPFILE = "20170112_1205.txt"
graw_to_std_fmt(INDIR, OUTDIR, SAMPFILE, PREFIX, LAT, LON, HEIGHT)
