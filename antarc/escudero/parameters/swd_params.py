#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:46:52 2022

@author: prowe
"""

# Directories
STAND_DIR = "/Users/prowe/sync/measurements/Escudero/pyranometer/stand/"
dir_out = "/Users/prowe/sync/measurements/Escudero/pyranometer/misc/"

# Files and formats
PYR_FILEFORMAT = "esc_swd%Y%m%d_%H%M.csv"

SW_FACTOR = [
    -1.0544898029995496e-07,
    1.8441615557710695e-05,
    -0.001282781420930488,
    0.0391867208184695,
    0.4850361690817559,
]

# EXtra stuff (delete)
# Parameters for pyranometer data files (shortwave downwelling)
# ORIG_SAMPLEFNAME = "LOG161207-2200.csv"
ORIG_FNAME_PREFIX = "LOG"
ORIG_FNAME_FORMAT = "LOG%y%m%d-%H%M.csv"
ORIG_DATE_FMT = "%Y-%m-%d"
ORIG_TIME_FMT = "%H:%M:%S"

STAND_FNAME_FORMAT = "esc_swd%Y%m%d_%H%M.csv"

YEARS = {
    "v1": range(2016, 2023),
    "v2": range(2022, 2023),
}

# Classes for versions
CLASSNAMES = {"v1": "StandFilesFromOrigV1", "v2": "StandFilesFromOrigV2"}

ESSENTIAL_COLS = [
    "Date",
    "Time",
    "#Samples",
    "Radiation",
]


# VERSION_START_FILES = {
#     "v1": "LOG161207-2200.csv",
#     "v2": "LOG220207-1906.csv",
# }

# ".data"
# ".data.hdr"
