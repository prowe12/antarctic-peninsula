#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 23:01:01 2022

@author: prowe
"""

# Dependencies
import datetime as dt
import numpy as np
import os

from antarc.sonde_raw_to_data_denial_format import graw_raw_to_datadenial
from antarc.get_atm_profs_rsrc import Sonde


indir = "tests/testdata/Escudero_graw_upp_raw/"
outdir = "tests/temp/"
testdatadir = "tests/testdata/Escudero_data_denial_format/"
sample_file = "20220205120020060702_UPP_RAW_89056_2022020512.txt"


graw_raw_to_datadenial(
    indir, outdir, sample_file, "Escudero", -62.2014, -58.9622, 60.0
)

sondefiles = [
    "Escudero_2022020700.txt",
    "Escudero_2022020512.txt",
    "Escudero_2022020712.txt",
    "Escudero_2022020800.txt",
    "Escudero_2022020823.txt",
]
sondefilefmt = "Escudero_%Y%m%d%H.txt"

for sondefile in sondefiles:
    sondedate = dt.datetime.strptime(sondefile, sondefilefmt)
    newsonde = Sonde(outdir, sondefile, sondedate)

    correctsonde = Sonde(testdatadir, sondefile, sondedate)

    # fields; P, T, date, file, h2o, rh, time, z
    assert np.allclose(newsonde.P, correctsonde.P)
    assert np.allclose(newsonde.T, correctsonde.T)
    assert np.allclose(newsonde.h2o, correctsonde.h2o)
    assert np.allclose(newsonde.rh, correctsonde.rh)
    assert np.allclose(newsonde.time, correctsonde.time)
    assert np.allclose(newsonde.z, correctsonde.z)
    assert newsonde.file == correctsonde.file

    # Delete the file from the temp dir
    os.remove(outdir + sondefile)
