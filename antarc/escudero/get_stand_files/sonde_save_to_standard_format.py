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


# My modules
# from antarc.graw_sonde_cols import (
#     graw_raw_to_datadenial,
#     compare_profiles,

from antarc.graw_sonde_cols import (
    get_cols_esc_eca53,
    get_cols_esc_eca54,
    get_cols_esc_eca55,
    get_cols_esc_eca56,
    get_cols_esc_yopp2022,
)

# get_cols_esc_eca53_old,

# Parameters
import antarc.params


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
    cols,
    units,
    colmap,
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

        # Skip two lines, which are the column labels and the units.
        # First, make sure they are as expected.
        # Require 2+ spaces between words
        # if re.split(r"\s{2,}", lines[i]) != cols:
        file_cols = [x.strip() for x in lines[i].split("\t")]

        if file_cols != cols:
            # file_col = re.split(r"\s{2,}", lines[i])
            logmsg = (
                genmsg
                + f": Columns differ; file has:\n{file_cols} \n\nnot\n\n{cols}"
            )
            # ", moving on"
            raise ValueError(logmsg)
            print(logmsg, file=lid)
            continue
        i += 1
        # if any(lines[i].split()[j] != units[j] for j in iunit):

        file_units = lines[i].split()
        if file_units != units:
            msg = f"On file {input_file}, the units are:\n"
            msg2 = f"\n{file_units}\n\nbut should be\n\n{units}"
            raise ValueError(msg + msg2)

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

        # Go through the rest of the record until we get to 'Trop' or a bad line
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
        itm = 0
        ymd = datestr[:4] + "-" + datestr[4:6] + "-" + datestr[6:8]

        # Write an output file with no qc
        noqcfile = output_dir + output_file
        with open(noqcfile, "w+", encoding="utf8") as fid:
            # Write the header
            fid.write("Date,                  P(hPa),   Alt(m),   Temp(C), ")
            fid.write("  RH(%),  Dewp(C), wdir(deg), wspd(kn)\n")

            for k in range(j):
                # Time as HH:MM:ss
                hour = int(np.floor(met[k, itm] / 3600))  # s x min/s x hr/min
                mint = int(np.floor(met[k, itm] / 60) - hour * 60)  # minutes
                sec = int(met[k, itm] - hour * 3600 - mint * 60)  # seconds
                hhmmss = f"_{hour:02d}:{mint:02d}:{sec:02d},"

                # p_z_t_rh = pztrhstr.format(press[k], alt[k], temp[k], rhw[k])
                # dewpt_wdir_spd = dewstr.format(dewpt[k], wdir[k], wspd[k])
                p_z_t_rh = pztrhstr.format(
                    met[k, 1], met[k, 2], met[k, 3], met[k, 4]
                )
                dewpt_wdir_spd = dewstr.format(met[k, 5], met[k, 6], met[k, 7])
                print(ymd + hhmmss + p_z_t_rh + dewpt_wdir_spd, file=fid)

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


# Parameters
from antarc.escudero.parameters import esc_params, sonde_params
from antarc.escudero.parameters import atm_prof_params


# Directories with radiosonde data and function to parse them
get_col_fun = {
    "/eca54_escudero_profiles_simulation2019/": get_cols_esc_eca54,  # 1
    "/eca55_escudero_profiles_a/": get_cols_esc_eca55,  # 2
    "/eca55_escudero_profiles_b/": get_cols_esc_eca54,  # 3
    "/eca56_escudero_profiles/": get_cols_esc_eca56,  # 4
    "/eca58_escudero_profiles/": get_cols_esc_yopp2022,  # 5
    "/yopp2022_escudero_profiles/": get_cols_esc_yopp2022,  # 6
    "/eca59_escudero_profiles/": get_cols_esc_yopp2022,  # 7
    "/simulation_reruns_2023/": get_cols_esc_yopp2022,
}

# Set variables from parameters
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
height = esc_params.ALTITUDE

outdir = sonde_params.STAND_DIR
no_qc_dir = sonde_params.NO_QC_DIR
sampfile = sonde_params.SAMPLEFNAME
prefix = "esc_sonde_v0"
fmt = atm_prof_params.SONDE_FILEFORMAT


sonde_dir = antarc.params.MEAS_DIR + "Escudero/GRAW_radiosondes/"
outdir = antarc.params.MEAS_DIR + "Escudero/radiosondes/graw/"


# Reading files back in:
# fname = '/Users/prowe/Sync/measurements/Escudero/radiosondes/graw/esc_sonde_dd_v1_2017011212.txt'


# Run the code to convert from the GRAW format to the data denial format
sampfile = sonde_params.SAMPLEFNAME
for indir in sonde_params.ORIG_DIRS[1:]:
    slash = indir[: indir.rfind("/")].rfind("/")
    get_cols = get_col_fun[indir[slash:]]
    cols, units, icol = get_cols()
    graw_to_std_fmt(
        indir, outdir, sampfile, prefix, lat, lon, height, cols, units, icol
    )


# For first field season, simulations often have high temps and winds of zero
# Whereas these values seem to use the kestrel surface readings
# So use these for no
# ECA53: 2017/01/12 - 2017/01/31
indir = sonde_dir + "soundings_txt_YYYYMMDD_hhmm/"
sampfile = "20170112_1205.txt"
get_cols = get_cols_esc_eca53
cols, units, icol = get_cols()
graw_to_std_fmt(
    indir, outdir, sampfile, prefix, lat, lon, height, cols, units, icol
)


# Plot the original sonde profiles vs version with minor QC
# savedir = outdir + "figures/"
# compare_profiles(outdir, no_qc_dir, fmt, "final", "noQC", savedir)


# def compare_profiles(
#     dir1: str, dir2: str, fmt: str, lab1: str, lab2: str, savedir: str
# ):
#     """
#     Create figures that compare radiosonde profiles between radiosonde results
#     in data denial format in two directories, e.g. where one has some QC and
#     The other does not
#     @param sampfile  Name of a sample file to match
#     @param dir1  The first directory
#     @param dir2  The second directory
#     @param fmt  The format of the files to be loaded in
#     @param lab1  Label for curves corresponding to dir1
#     @param lab2  Label for curves corresponding to dir2
#     @param savedir  Directory to save created figures
#     """

#     # Get filenames
#     sonde_file_n_dates = get_files(dir1, fmt)

#     for sonde_file, sonde_date in sonde_file_n_dates:
#         if not exists(dir2 + sonde_file):
#             continue
#         print("Working on", sonde_file)
#         sonde1 = Sonde(dir1, sonde_file, sonde_date)
#         sonde2 = Sonde(dir2, sonde_file, sonde_date)

#         zmax = min(max(np.nanmax(sonde1.z), np.nanmax(sonde2.z)), 30) + 1
#         h2omax = (
#             min(max(np.nanmax(sonde1.h2o), np.nanmax(sonde2.h2o)), 10000) + 100
#         )
#         rhmax = min(max(np.nanmax(sonde1.rh), np.nanmax(sonde2.rh)), 120) + 2
#         tmax = min(max(np.nanmax(sonde1.T), np.nanmax(sonde2.T)), 300) + 5
#         tmin = max(min(np.nanmin(sonde1.T), np.nanmin(sonde2.T)), 190) - 10

#         plt.figure(1)
#         plt.clf()
#         plt.subplot(221)
#         plt.plot(sonde2.P, sonde2.z, ".-", label=lab2)
#         plt.plot(sonde1.P, sonde1.z, label=lab1)
#         plt.title(sonde_file)
#         plt.xlabel("Pressure")
#         plt.ylabel("Altitude (km)")
#         plt.legend()
#         plt.ylim([-0.4, zmax])

#         plt.subplot(222)
#         plt.plot(sonde2.T, sonde2.z, ".-")
#         plt.plot(sonde1.T, sonde1.z)
#         plt.xlabel("Temperature")
#         plt.ylabel("Altitude (km)")
#         plt.ylim([-0.4, zmax])
#         plt.xlim([tmin, tmax])

#         plt.subplot(223)
#         plt.plot(sonde2.h2o, sonde2.z, ".-")
#         plt.plot(sonde1.h2o, sonde1.z)
#         plt.plot([-100, h2omax], [0, 0], "k:")
#         plt.plot([0, 0], [-0.4, zmax], "k:")
#         plt.xlabel("H2O")
#         plt.xlim([-100, h2omax])
#         plt.ylim([-0.4, zmax])

#         plt.subplot(224)
#         plt.plot(sonde2.rh, sonde2.z, ".-")
#         plt.plot(sonde1.rh, sonde1.z)
#         plt.xlabel("RH (%)")
#         plt.xlim([-2, rhmax])
#         plt.ylim([-0.4, zmax])

#         figname = sonde_file[:-3] + "png"
#         plt.savefig(savedir + figname)
