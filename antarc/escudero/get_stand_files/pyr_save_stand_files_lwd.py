#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 12:29:24 2022

@author: prowe
"""

import logging

logging.basicConfig(
    filename="example.log", encoding="utf-8", level=logging.DEBUG
)


from antarc.pyr_save_stand_files import pyr_save_stand_files


# LWD: pyrgeometer
paramfile = "antarc.escudero.parameters.lwd_orig_params"
logging.info("")
logging.info("LWD (pyrgeometer)")
logging.info("Version 1")
pyr_save_stand_files(paramfile, "v1", 2017)
pyr_save_stand_files(paramfile, "v1", 2018)
pyr_save_stand_files(paramfile, "v1", 2019)
pyr_save_stand_files(paramfile, "v1", 2020)
pyr_save_stand_files(paramfile, "v1", 2021)
pyr_save_stand_files(paramfile, "v1", 2022)

logging.info("")
logging.info("Version 2")
pyr_save_stand_files(paramfile, "v2", 2022)
pyr_save_stand_files(paramfile, "v2", 2023)
