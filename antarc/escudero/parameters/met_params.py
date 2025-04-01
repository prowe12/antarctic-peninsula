#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:22:11 2022

@author: prowe

"""

from antarc import params


STATION = 950001

# Directories
MEAS_DIR = params.MEAS_DIR
DIREC_BY_HOUR = MEAS_DIR + "Escudero/frei_met/by_hour/"
DIREC_BY_MINUTE = MEAS_DIR + "Escudero/frei_met/by_minute/"
DIR_15MIN = MEAS_DIR + "Escudero/frei_met/by_month_every_15min/"
PFIX = "esc_met_"
EXT = ".csv"

LAT = -62.19194
LON = -58.97972
ALT = 45
