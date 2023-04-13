#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:20:11 2022

@author: prowe
"""
# # # # # # # # #    INPUTS   # # # # # # # # # #
# Directories
main_dir = "/Users/prowe/Projects/NSF_AERI_cloud_CRI/McMurdo/case_studies/"
yopp_dir = "/Users/prowe/Projects/YOPP/"
CO2_DIR = "/Users/prowe/Sync/measurements/Escudero/co2/"
SONDE_DIR = "/Users/prowe/Sync/measurements/Escudero/radiosondes/datadenial/"
ERA_DIR = "/Users/prowe/Sync/measurements/Escudero/era5/"
OUT_DIR = "/Users/prowe/Sync/measurements/Escudero/profiles/"
FIG_DIR = "/Users/prowe/Sync/measurements/Escudero/profiles//figures/"
surfmet_dir = main_dir + "met/"

# Filenames and formats
SONDE_FILEFORMAT = "esc_sonde_dd_%Y%m%d%H.txt"
# CT_FILEFORMAT = "CT-NRT.v2022-1.molefrac_glb3x2_%Y-%m-%d.nc"
CT_FILEFORMAT = "CT-NRT.v2023-2.molefrac_glb3x2_%Y-%m-%d.nc"
ERA_FILEFORMAT = "era5_esc_%Y%m%d.nc"

CO2_FILE_PALMER = CO2_DIR + "co2_psa_surface-flask_1_ccgg_event.txt"
CO2_FILE_DRAKE = CO2_DIR + "co2_drp_shipboard-flask_1_ccgg_event.txt"
# surfmet_file = 'awrmetM1.b1.20160208.000000.cdf'

MODEL_EXTRA = 3  # model number for unspecified moledule concentrations

# Layer boundaries (refine as needed for desired temperature differential)
LAYERBNDS = [
    0.0,
    0.025,
    0.05,
    0.1,
    0.15,
    0.2,
    0.25,
    0.3,
    0.35,
    0.4,
    0.45,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.1,
    1.25,
    1.4,
    1.6,
    1.8,
    2.0,
    2.2,
    2.4,
    2.6,
    2.8,
    3.0,
    3.1,
    3.15,
    3.42,
    3.65,
    3.925,
    4.2,
    4.475,
    4.75,
    5.0,
    5.5,
    6.0,
    6.5,
    7.0,
    7.5,
    8.0,
    8.5,
    9.0,
    9.5,
    10.0,
    11.0,
    12.0,
    13.0,
    14.0,
    15.0,
    16.0,
    17.0,
    19.0,
    21.0,
    23.0,
    25.0,
    28.0,
    30.5,
    32.5,
    34.4,
    36.4,
    38.2,
    40.0,
    41.9,
    44.0,
    46.0,
    48.5,
    60.0,
]
