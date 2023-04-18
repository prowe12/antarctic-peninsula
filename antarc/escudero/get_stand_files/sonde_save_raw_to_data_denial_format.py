#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 23:01:01 2022

@author: prowe
"""

# My modules
from antarc.sonde_raw_to_data_denial_format import (
    graw_raw_to_datadenial,
    compare_profiles,
)

# Parameters
from antarc.escudero.parameters import esc_params, sonde_params
from antarc.escudero.parameters import atm_prof_params

# Set variables from parameters
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
height = esc_params.ALTITUDE

outdir = sonde_params.STAND_DIR
no_qc_dir = sonde_params.NO_QC_DIR
samplefile = sonde_params.SAMPLEFNAME
prefix = sonde_params.PREFIX_FOR_STANDFILE
fmt = atm_prof_params.SONDE_FILEFORMAT

# Run the code to convert from the GRAW format to the data denial format
for indir in sonde_params.ORIG_DIRS:
    graw_raw_to_datadenial(indir, outdir, samplefile, prefix, lat, lon, height)

# Plot the original sonde profiles vs version with minor QC
savedir = outdir + "figures/"
compare_profiles(outdir, no_qc_dir, fmt, "final", "noQC", savedir)
