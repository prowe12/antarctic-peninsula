#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:46:52 2022

@author: prowe
"""

# Parameter modules
from antarc import params

MEAS_DIR = params.MEAS_DIR

# Directories
SWD_DIR = MEAS_DIR + "Escudero/pyranometer/stand/"
SWD_CLEAR_DIR = MEAS_DIR + "Escudero/swd_clear/"
# STAND_DIR = MEAS_DIR + "Escudero/pyranometer/stand/"
# dir_out = MEAS_DIR + "Escudero/pyranometer/misc/"

# Files and formats
SWD_FILEFORMAT = "esc_swd%Y%m%d_%H%M.csv"
SWD_CLEAR_FILEFORMAT = "flux%Y%m%d_%H%M.txt"
# PYR_FILEFORMAT = "esc_swd%Y%m%d_%H%M.csv"


# SW_FACTOR = [
#     -1.0544898029995496e-07,
#     1.8441615557710695e-05,
#     -0.001282781420930488,
#     0.0391867208184695,
#     0.4850361690817559,
# ]
