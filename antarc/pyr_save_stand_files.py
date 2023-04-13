#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 11:24:36 2022

@author: prowe
"""

import datetime as dt

from antarc.get_pyr_stand_files_from_orig import (
    StandFilesFromOrigV1,
    StandFilesFromOrigV2,
    LwdStandFilesFromOrigV1,
    LwdStandFilesFromOrigV2,
)
from antarc.getfilenames import getfilenames
import importlib


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
    essential_cols,
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
        stdFromOrig = StdFrOrig(dir_in, filename, in_file_fmt, essential_cols)
        if stdFromOrig.get_qc_flag():
            outfile = stdFromOrig.get_output_filename(
                date_fmt, time_fmt, out_file_fmt
            )
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

    classnames = {
        "StandFilesFromOrigV1": StandFilesFromOrigV1,
        "StandFilesFromOrigV2": StandFilesFromOrigV2,
        "LwdStandFilesFromOrigV1": LwdStandFilesFromOrigV1,
        "LwdStandFilesFromOrigV2": LwdStandFilesFromOrigV2,
    }

    params = importlib.import_module(paramfile)

    # Make sure the desired year exists for the desired version
    if year not in params.YEARS[version]:
        import warnings

        warnings.warn(f"{year} not found in {params.YEARS[version]}, skipping")
        return

    # Create standardized files from Version 1 original files
    load_origs_save_stands(
        params.ORIG_DIR[version],
        params.STAND_DIR,
        year,
        params.ORIG_FNAME_FORMAT,
        params.STAND_FNAME_FORMAT,
        classnames[params.CLASSNAMES[version]],
        params.ORIG_FNAME_PREFIX,
        params.ORIG_DATE_FMT,
        params.ORIG_TIME_FMT,
        params.ESSENTIAL_COLS,
    )
