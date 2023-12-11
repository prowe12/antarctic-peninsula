#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 12:31:01 2023

@author: prowe
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Directories
maindir = "/Users/prowe/Sync/"
case_dir = maindir + "projects/NSF_AP/case_studies/"
indir = maindir + "measurements/Escudero/ucass/20170131_18/"
outdir = case_dir + "model_performance_KGI_2017_01/figures/"
fname = "2017013118512935_X_layers_Reff_Tot-Cn_Cm_weighted.dat"

# Runtime parameters
save_figs = False

ucass = pd.read_csv(indir + fname, sep=r"\s+")

# Variables in si units (kg and m)
alt = ucass["layer_center[m]"].to_numpy()
reff = ucass["reff[micron]"].to_numpy() / 1e6  # m
c_n = ucass["num_conc[m-3]"].to_numpy()
lwc = ucass["lwc[kg/m3]"].to_numpy()  # kg/m3

print(f"The LWP is {round(np.trapz(lwc[1:], alt[1:])*1000,1)} g/m2")

# Make plot of Reff, Cn, LWC
fsz = 14
fig, axs = plt.subplots(1, 3, layout="constrained")
ax = axs[0]
ax.plot(reff * 1e6, alt)
ax.set_ylabel("Altitude (m)")
ax.set_xlabel(r"Effective radius ($\mu$m)")
ax.text(5.9, 2010, "a)", fontsize=fsz)

ax = axs[1]
ax.plot(c_n / (100**3), alt)
ax.set_xlabel("C$_n$ (cm$^{-3}$)")
ax.set(yticklabels=[])
ax.text(590, 2010, "b)", fontsize=fsz)

ax = axs[2]
ax.plot(lwc * 1000, alt)
ax.set_xlabel("LWC (g/m$^{-3}$)")
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.set(yticklabels=[])
ax.text(0.35, 2010, "c)", fontsize=fsz)

# Just show the cloud layers
axs[0].set_ylim([2000, 2500])
axs[1].set_ylim([2000, 2500])
axs[2].set_ylim([2000, 2500])


if save_figs:
    plt.savefig(outdir + "reff_cn_lwc.png")
    plt.savefig(outdir + "reff_cn_lwc.eps", format="eps")

# Compute the optical depth in the visible region
qext = 2
ro_w = 1000  # kg/m3
tau_z = 3 / 4 * qext * lwc / ro_w / reff  # kg/m3 / (kg/m3) / m => /m
tau = tau_z * 50  # /m * m

alt_stairs = []
tau_stairs = []
reff_stairs = []
alt_stairs.append(0)
for alt0 in alt:
    alt_edge = alt0 + 25
    alt_stairs.append(alt_edge)
    alt_stairs.append(alt_edge)
    tau_stairs.append(tau[alt == alt0][0])
    tau_stairs.append(tau[alt == alt0][0])
    reff_stairs.append(1e6 * reff[alt == alt0][0])
    reff_stairs.append(1e6 * reff[alt == alt0][0])
alt_stairs.pop()


fig, axs = plt.subplots(1, 2, layout="constrained")
axs[0].plot([0, 0], [2000, 2500], "k:")
axs[0].plot(tau_stairs, alt_stairs)
axs[0].plot(tau[alt >= 2000], alt[alt >= 2000], "*")
axs[0].set_xlabel("optical depth")
axs[0].axis([-0.1, 6, 2050, 2450])
axs[0].set_ylabel("Altitude (m)")

axs[1].plot(reff_stairs, alt_stairs)
axs[1].plot(1e6 * reff[alt >= 2000], alt[alt >= 2000], "o")
axs[1].set_xlabel(r"effective radius ($\mu$m)")
axs[1].axis([2, 7, 2050, 2450])
axs[1].yaxis.set_label_position("right")
axs[1].yaxis.tick_right()

if save_figs:
    plt.savefig(outdir + "reff_opticaldepth_stairs.png")
    plt.savefig(outdir + "reff_opticaldepth_stairs.eps", format="eps")
