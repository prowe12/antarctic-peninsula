#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:46:52 2022

@author: prowe
"""

# Parameter module, specific to user
from antarc import params

# Directories specific to user
MEAS_DIR = params.MEAS_DIR

# Directories
SONDE_DIR = MEAS_DIR + "Escudero/GRAW_radiosondes/"
ORIG_DIRS = (
    SONDE_DIR + "simulation/eca53_profiles_simulation_2019run/",
    SONDE_DIR + "simulation/eca54_escudero_profiles_simulation2019/",
    SONDE_DIR + "simulation/eca55_escudero_profiles/",
    SONDE_DIR + "simulation/eca56_escudero_profiles/",
    SONDE_DIR + "simulation/eca58_escudero_profiles/",
    SONDE_DIR + "simulation/yopp2022_escudero_profiles/",
)
SAMPLEFNAME = "20220401120020047677_UPP_RAW_89056_2022040112.txt"

STAND_DIR = MEAS_DIR + "Escudero/radiosondes/datadenial/"
NO_QC_DIR = STAND_DIR + "noQC/"
PREFIX_FOR_STANDFILE = "esc_sonde_dd_v1"
