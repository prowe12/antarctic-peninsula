#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:46:52 2022

@author: prowe
"""

from antarc import params


# Directories
MEAS_DIR = params.MEAS_DIR
LWD_DIR = MEAS_DIR + "Escudero/pyrgeometer/stand/"
LWD_CLEAR_DIR = MEAS_DIR + "measurements/Escudero/lwd_clear/"

# Files and formats
LWD_FILEFORMAT = "esc_lwd%Y%m%d_%H%M.csv"
LWD_CLEAR_FILEFORMAT_OUTNU = "flux%Y%m%d_%H%M_0_238.txt"
LWD_CLEAR_FILEFORMAT_LOWNU = "flux%Y%m%d_%H%M_238_1500.txt"
LWD_CLEAR_FILEFORMAT_HINU = "flux%Y%m%d_%H%M_1500_2222.txt"
