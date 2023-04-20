#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:46:52 2022

@author: prowe
"""

from antarc import params

# Directories
SONDE_DIR = params.MEAS_DIR + "KingSejong/graw_radiosondes/"
ORIG_DIRS = (SONDE_DIR + "sounding_data/",)
SAMPLEFNAME = "20220330180019032482_UPP_RAW_89251_2022033018.txt"

STAND_DIR = params.MEAS_DIR + "KingSejong/radiosondes/datadenial/"
PREFIX_FOR_STANDFILE = "ksj_sonde_dd_v0"
SONDE_FILEFORMAT = "ksj_sonde_dd_v0_%Y%m%d%H.txt"
