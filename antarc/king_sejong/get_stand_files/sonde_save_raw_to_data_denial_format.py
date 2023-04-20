#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 23:01:01 2022

@author: prowe
"""

# My modules
from antarc.king_sejong.get_stand_files.sonde_raw_to_data_denial_format import (
    graw_raw_to_datadenial,
)


# Parameters
from antarc.king_sejong.parameters import ksj_params, sonde_params

# Set variables from parameters
lat = ksj_params.LATITUDE
lon = ksj_params.LONGITUDE
height = ksj_params.ALTITUDE

indir = sonde_params.ORIG_DIRS
outdir = sonde_params.STAND_DIR
# no_qc_dir = sonde_params.NO_QC_DIR
samplefile = sonde_params.SAMPLEFNAME
prefix = sonde_params.PREFIX_FOR_STANDFILE
fmt = sonde_params.SONDE_FILEFORMAT


# Run the code to convert from the GRAW format to the data denial format
for indir in sonde_params.ORIG_DIRS:
    graw_raw_to_datadenial(indir, outdir, samplefile, prefix, lat, lon, height)

# # Plot the original sonde profiles vs version with minor QC
# savedir = outdir + "figures/"
# compare_profiles(outdir, no_qc_dir, fmt, "final", "noQC", savedir)
