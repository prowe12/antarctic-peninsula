#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:18:12 2023

@author: prowe
"""

# Built-in modules
import numpy as np

# My modules
from antarc import params


# Directories
MEAS_DIR = params.MEAS_DIR
esc_dir = MEAS_DIR + "Escudero/"
outdir1 = esc_dir + "od/od_20170131_18a/"
outdir2 = esc_dir + "od/od_20170131_18b/"
outdir3 = esc_dir + "od/od_20170131_18c/"


# parameters
plot_results = False
# Within this range the below-cloud transmittance is <= 1e-6
# so we can just use the clear-sky fluxes
bwn1 = 0.0
ewn1 = 360.0
# Within these ranges the below-cloud transmittance is > 1e-6
# so we need to do the DISORT calculation
bwn2 = 360.0
ewn2 = 1484.7
bwn3 = 1484.7
ewn3 = 2222.2  # pyrgeometer cut-off

# Files
# Compute the clear-sky fluxes as functions of wavenumber
# file_tot1 = f"flux_tot_{int(bwn1)}_{int(ewn1)}.txt"
# file_tot2 = f"flux_tot_{int(bwn2)}_{int(ewn2)}.txt"
# file_tot3 = f"flux_tot_{int(bwn3)}_{int(ewn3)}.txt"
file_nu1 = f"flux_nu_{int(bwn1)}_{int(ewn1)}.txt"
file_nu2 = f"flux_nu_{int(bwn2)}_{int(ewn2)}.txt"
file_nu3 = f"flux_nu_{int(bwn3)}_{int(ewn3)}.txt"



# Load the clear-sky results from LBLRTM
nu_lbl1, flux_lbl1 = np.loadtxt(outdir1 + file_nu1)
nu_lbl2, flux_lbl2 = np.loadtxt(outdir2 + file_nu2)
nu_lbl3, flux_lbl3 = np.loadtxt(outdir3 + file_nu3)

totflux_lbl1 = np.trapz(flux_lbl1, nu_lbl1)
totflux_lbl2 = np.trapz(flux_lbl2, nu_lbl2)
totflux_lbl3 = np.trapz(flux_lbl3, nu_lbl3)


# Load the DISORT results
totflux_dis1 = totflux_lbl1
totflux_dis2 = 0
totflux_dis3 = 0
totflux_dis_cld1 = totflux_lbl1
totflux_dis_cld2 = 0
totflux_dis_cld3 = 0

plot_results = True
if plot_results:
    import matplotlib.pyplot as plt
    
    plt.figure()
    plt.plot(nu_lbl1, flux_lbl1, label = 'LBLRTM')
    plt.plot(nu_lbl2, flux_lbl2, label = 'LBLRTM')
    plt.plot(nu_lbl3, flux_lbl3, label = 'LBLRTM')

# wavenumber region 2
for i in range(10):
    ind = i*10000
    nu_dis2, flux_dis2 = np.loadtxt(f"{outdir2}flux_disort_{ind}")
    totflux_dis2 += np.trapz(flux_dis2, nu_dis2)
    nu_dis2cld, flux_dis2cld = np.loadtxt(f"{outdir2}flux_disort_cldy_{ind}")
    totflux_dis_cld2 += np.trapz(flux_dis2cld, nu_dis2cld)
    if plot_results:
        plt.plot(nu_dis2, flux_dis2, label = 'DISORT')
        plt.plot(nu_dis2cld, flux_dis2cld, label = 'DISORT')

# wavenumber region 3
for i in range(8):
    ind = i*10000
    nu_dis3, flux_dis3 = np.loadtxt(f"{outdir3}flux_disort_{ind}")
    totflux_dis3 += np.trapz(flux_dis3, nu_dis3)
    nu_dis3cld, flux_dis3cld = np.loadtxt(f"{outdir3}flux_disort_cldy_{ind}")
    totflux_dis_cld3 += np.trapz(flux_dis3cld, nu_dis3cld)
    if plot_results:
        plt.plot(nu_dis3, flux_dis3)
        plt.plot(nu_dis3cld, flux_dis3cld)


# Compare the LBLRTM and DIOSRT clear-sky and the DISORT cloudy-sky results
print(f"LBLRTM 1: {totflux_lbl1}")
print(f"LBLRTM 2: {totflux_lbl2}")
print(f"DISORT 2: {totflux_dis2}")
print(f"CLOUDY 2: {totflux_dis_cld2}")
print(f"LBLRTM 3: {totflux_lbl3}")
print(f"DISORT 3: {totflux_dis3}")
print(f"CLOUDY 3: {totflux_dis_cld3}")

totflux_lbl = totflux_lbl1 + totflux_lbl2 + totflux_lbl3
totflux_dis = totflux_lbl1 + totflux_dis2 + totflux_dis3
totflux_dis_cld = totflux_lbl1 + totflux_dis_cld2 + totflux_dis_cld3
print(f"LBLRTM total flux: {totflux_lbl}")
print(f"DISORT total flux: {totflux_dis}")
print(f"CLOUDY total flux: {totflux_dis_cld}")
print(f"DISORT - LBLRTM: {totflux_dis - totflux_lbl}")
print(f"CLOUDY - LBLRTM: {totflux_dis_cld - totflux_lbl}")

