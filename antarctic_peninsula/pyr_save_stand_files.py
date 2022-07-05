#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 11:24:36 2022

@author: prowe
"""

import datetime as dt

from pyr_get_stand_files_from_orig import (
    StandFilesFromOrigV1,
    StandFilesFromOrigV2,
    LwdFilesFromOrigV1,
    LwdFilesFromOrigV2,
)
from getfilenames import getfilenames


def load_origs_save_stands(
    in_dir,
    out_dir,
    year,
    in_file_fmt,
    out_file_fmt,
    StdFrOrig,
    in_pfix,
    date_fmt,
    time_fmt,
):
    """
    Load files in original format and save in standard format, for a year of data
    @param dir_pyr  Main directory (data must be in subdirectory for year)
    @param dir_stand  Directory for standard ouput (will sub into year)
    @param year  Desired year
    @param samplefname  An example file name
    @param StandFromOrig  The version of StandFilesFromOrig to use
    """
    fnamelen = len(dt.datetime.strftime(dt.datetime.now(), in_file_fmt))
    yeardir = f"{year}/"
    dir_in = in_dir + yeardir
    dir_out = out_dir + yeardir
    fnames = getfilenames(dir_in, fnamelen, in_pfix)

    for filename in fnames:
        stdFromOrig = StdFrOrig(dir_in, filename, in_file_fmt)
        if stdFromOrig.isempty():
            continue
        # datestr = stdFromOrig.get_datestr()
        # timestr = stdFromOrig.get_first_timestr()
        # outfile = "esc_swd20" + datestr + "_" + timestr + ".csv"
        outfile = stdFromOrig.get_output_filename(date_fmt, time_fmt, out_file_fmt)
        stdFromOrig.save_dataframe(dir_out + outfile)


def pyr_save_stand_files(paramfile: str, version: str, year: int):
    """
    Save the data from the original format (based on version number)
    to the standard format
    paramfile  File with all the neeeded inputs
    
    Example inputs:

        # Pyranometer
    paramfile = "parameters.swd_orig_params"

    # Pyrgeometer
    paramfile = "parameters.lwd_orig_params"

    """

    # Classes for versions
    CLASSNAMES = {"v1": StandFilesFromOrigV1, "v2": StandFilesFromOrigV2}

    params = __import__(paramfile, fromlist=["nan"])

    # Make sure the desired year exists for the desired version
    if str(year) not in params.YEARS:
        print(f"{year} not found in params.YEARS: {params.YEARS}")

    # Create standardized files from Version 1 original files
    load_origs_save_stands(
        params.ORIG_DIR[version],
        params.STAND_DIR,
        year,
        params.ORIG_FNAME_FORMAT,
        params.STAND_FNAME_FORMAT,
        CLASSNAMES[version],
        params.ORIG_FNAME_PREFIX,
        params.ORIG_DATE_FMT,
        params.ORIG_TIME_FMT,
    )
