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
    get_cols_esc_eca53_54,
    get_cols_esc_eca55,
    get_cols_esc_eca56,
    get_cols_esc_yopp2022,
)

# Parameters
from antarc.escudero.parameters import esc_params, sonde_params
from antarc.escudero.parameters import atm_prof_params


# TODO: move to sonde_params file?
get_col_fun = {
    "/eca53_profiles_simulation_2019run/": get_cols_esc_eca53_54,  # 0
    "/eca54_escudero_profiles_simulation2019/": get_cols_esc_eca53_54,  # 1
    "/eca55_escudero_profiles_a/": get_cols_esc_eca55,  # 2
    "/eca55_escudero_profiles_b/": get_cols_esc_eca53_54,  # 3
    "/eca56_escudero_profiles/": get_cols_esc_eca56,  # 4
    "/eca58_escudero_profiles/": get_cols_esc_yopp2022,  # 5
    "/yopp2022_escudero_profiles/": get_cols_esc_yopp2022,  # 6
    "/eca59_escudero_profiles/": get_cols_esc_yopp2022,  # 7
}

# Set variables from parameters
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
height = esc_params.ALTITUDE

outdir = sonde_params.STAND_DIR
no_qc_dir = sonde_params.NO_QC_DIR
sampfile = sonde_params.SAMPLEFNAME
prefix = sonde_params.PREFIX_FOR_STANDFILE
fmt = atm_prof_params.SONDE_FILEFORMAT

# Run the code to convert from the GRAW format to the data denial format
for indir in sonde_params.ORIG_DIRS[7:]:
    slash = indir[: indir.rfind("/")].rfind("/")
    get_cols = get_col_fun[indir[slash:]]
    cols, units, icol = get_cols()
    graw_raw_to_datadenial(
        indir, outdir, sampfile, prefix, lat, lon, height, cols, units, icol
    )

# Plot the original sonde profiles vs version with minor QC
# savedir = outdir + "figures/"
# compare_profiles(outdir, no_qc_dir, fmt, "final", "noQC", savedir)
