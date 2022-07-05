#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:46:52 2022

@author: prowe
"""

# Directories
STAND_DIR = "/Users/prowe/sync/measurements/Escudero/pyrgeometer/stand/"
ORIG_DIR_V1 = "/Users/prowe/sync/measurements/Escudero/pyrgeometer/orig_v1/"
ORIG_DIR_V2 = "/Users/prowe/sync/measurements/Escudero/pyrgeometer/orig_v2/"

# Parameters for pyranometer data files (shortwave downwelling)
# ORIG_SAMPLEFNAME = "LOG161207-2200.csv"
ORIG_FNAME_PREFIX = "LOG"
ORIG_FNAME_FORMAT = "LOG%y%m%d-%H%M.csv"
ORIG_DATE_FMT = "%Y-%m-%d"
ORIG_TIME_FMT = "%H:%M:%S"

STAND_FNAME_FORMAT = "esc_lwd%Y%m%d_%H%M.csv"

YEARS_V1 = range(2017, 2023)
YEARS_V2 = range(2022, 2023)


# ".data"
# ".data.hdr"
