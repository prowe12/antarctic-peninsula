#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:20:11 2022

@author: prowe
"""
# # # # # # # # #    INPUTS   # # # # # # # # # #
# Parameter modules - specific to user
from antarc import params

MEAS_DIR = params.MEAS_DIR
PROJ_DIR = params.PROJECT_DIR

# Directories
ERA_DIR = MEAS_DIR + "Rothera/era5/"
ERA_FILEFORMAT = "era5_roth_%Y%m%d.nc"
