#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:20:11 2022

@author: prowe

Sondes:
2017/01/12 - 2017/01/31: Done
2017/11/28 - 2017/12/13: Done
2018/01/04 - 2018/02/20: Done
2018/11/30 - 2018/12/31: Done
2019/01/01 - 2019/02/22: Done
2019/09/02 - 2019/10/03: Done
2020: No soundings
2021: No soundings
2022/01/21 - 2022/08/28: Done
2023/03/20 - 2023/08/10: Carbon Tracker not available yet (2023/08/23)

"""
# # # # # # # # #    INPUTS   # # # # # # # # # #
# Parameter modules - specific to user
from antarc import params

MEAS_DIR = params.MEAS_DIR
PROJ_DIR = params.PROJECT_DIR

# Directories
SONDE_DIR = MEAS_DIR + "Escudero/radiosondes/datadenial/"
ERA_DIR = MEAS_DIR + "Escudero/era5/"
OUT_DIR = MEAS_DIR + "Escudero/profiles/"
FIG_DIR = MEAS_DIR + "Escudero/profiles//figures/"
FREI_MET_DIR = MEAS_DIR + "Escudero/frei_met/by_month_every_15min/"
TARP_MET_DIR = MEAS_DIR + "Escudero/aws/std/"
CO2_DIR = MEAS_DIR + "Escudero/co2/"
CT_DIR = {
    2017: MEAS_DIR + "Escudero/co2/ct_2022/",
    2018: MEAS_DIR + "Escudero/co2/ct_2022/",
    2019: MEAS_DIR + "Escudero/co2/ct_2022/",
    2022: MEAS_DIR + "Escudero/co2/ct_nrt_2023_4/",
    2023: MEAS_DIR + "Escudero/co2/ct_nrt_2023_4/",
}

# Filenames and formats
SONDE_FILEFORMAT = "esc_sonde_dd_v2_%Y%m%d%H.txt"
CT_FILEFORMAT = {
    2017: "CT2022.molefrac_glb3x2_%Y-%m-%d.nc",
    2018: "CT2022.molefrac_glb3x2_%Y-%m-%d.nc",
    2019: "CT2022.molefrac_glb3x2_%Y-%m-%d.nc",
    2022: "CT-NRT.v2023-4.molefrac_glb3x2_%Y-%m-%d.nc",
    2023: "CT-NRT.v2023-4.molefrac_glb3x2_%Y-%m-%d.nc",
}
ERA_FILEFORMAT = "era5_esc_%Y%m%d.nc"
CO2_FILE_PALMER = CO2_DIR + "co2_psa_surface-flask_1_ccgg_event.txt"
CO2_FILE_DRAKE = CO2_DIR + "co2_drp_shipboard-flask_1_ccgg_event.txt"
# FREI_MET_FILE = "esc_met_201901.txt"
# TARP_MET_FILE = "950007_20210215.json"

TARP_MET_FILEFORMAT = "tarp02_met_%Y%m.txt"

# FREI
FREI_MET_FILEFORMAT = "esc_met_%Y%m.txt"
FREI_MET_STATION = 950001
FREI_MET_LAT = -62.19194
FREI_MET_LON = -58.97972
FREI_MET_ALT = 45.0 / 1000  # km

MODEL_EXTRA = 3  # model number for unspecified moledule concentrations

# Layer boundaries (refine as needed for desired temperature differential)
# fmt: off
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
# fmt: on
