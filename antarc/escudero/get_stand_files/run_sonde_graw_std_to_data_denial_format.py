#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 23:01:01 2022

@author: prowe
"""

# My modules
from antarc.escudero.get_stand_files.sonde_qc_graw_std import (
    qc_graw_std,
)

# Parameters
from antarc.escudero.parameters import esc_params, sonde_params

# from antarc.escudero.parameters import atm_prof_params


# Set variables from parameters
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
height = esc_params.ALTITUDE

indir = sonde_params.params.MEAS_DIR + "/Escudero/radiosondes/graw/"
outdir = sonde_params.STAND_DIR
prefix = sonde_params.PREFIX_FOR_STANDFILE

SAMPFILE = "esc_sonde_v0_2017011212.txt"

# Run the code to perform QC and convert from standardized GRAW format
# to data denial format
qc_graw_std(indir, outdir, SAMPFILE, prefix, lat, lon)

# Plot the original sonde profiles vs version with minor QC
# savedir = outdir + "figures/"
# compare_profiles(outdir, no_qc_dir, fmt, "final", "noQC", savedir)
