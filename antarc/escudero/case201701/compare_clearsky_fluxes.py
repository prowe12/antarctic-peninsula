#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:18:12 2023

@author: prowe
"""

# Built-in modules
import datetime as dt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytz


# My modules
from antarc import params
from antarc.escudero.parameters import esc_params
from antarc.escudero.all.analysis_rsrc import get_era5_escudero_byyear
from antarc.lblrtm_utils import Prof
from antarc.blackbody_cloud_flux import get_bb_fluxes
from antarc.escudero.parameters import swd_params
from antarc.load_rad_flux import load_rad_flux


def rnd3(val: float) -> float:
    """
    Return value rounded to 3rd decimal point
    @param val  Value to round
    @returns float  Rounded value
    """
    return round(val, 3)


# Directories
proj_dir = params.PROJECT_DIR + "case_studies/model_performance_KGI_2017_01/"
prof_dir = proj_dir + "ground_measurements/Escudero/profiles/"
oddir1 = proj_dir + "ground_measurements/Escudero/od/od_20170131_18a/"
oddir2 = proj_dir + "ground_measurements/Escudero/od/od_20170131_18b/"
oddir3 = proj_dir + "ground_measurements/Escudero/od/od_20170131_18c/"
ksj_dir = proj_dir + "ground_measurements/KingSejong/"
era_dir = proj_dir + "models/ERA/"
outdir = proj_dir + "figures/"

# Files
ksj_met_file = "KingSejong_Met_Jan31_2017.csv"
era_bbnd_file = "era5_esc_broadband_20170131.nc"
prof_file = prof_dir + "prof20170131_1801.nc"

# Output files
flux_fig = "downwelling_flux.png"
computed_flux_file = "computed_lwd_fluxes.txt"
esc_swd_file = "escudero_swd.txt"
era5_bbnd_file = "era5_broadband.txt"
prof_file_out = "radiosonde.txt"

# Runtime parameters
plot_results = True
save_results = False
# # # # # # # # # # # #


era_fmt = "era5_esc_broadband_*[0-9]*.nc"
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE

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
file_nu1 = f"flux_nu_{int(bwn1)}_{int(ewn1)}.txt"
file_nu2 = f"flux_nu_{int(bwn2)}_{int(ewn2)}.txt"
file_nu3 = f"flux_nu_{int(bwn3)}_{int(ewn3)}.txt"

snd_date = dt.datetime(2017, 1, 31, 18, 20, 0, tzinfo=pytz.utc)

# Load pyranometer data
date1 = dt.datetime(2017, 1, 31, tzinfo=pytz.utc)
date2 = dt.datetime(2017, 2, 1, tzinfo=pytz.utc)
swd_dir = f"{swd_params.SWD_DIR}{date1.strftime('%Y')}/"
swd_fmt = swd_params.SWD_FILEFORMAT
swd_date0, swd = load_rad_flux(swd_dir, swd_fmt, date1, date2)
swd_date = np.array(swd_date0)
swd[swd < 0] = 0


# Load the clear-sky results from LBLRTM
nu_lbl1, flux_lbl1 = np.loadtxt(oddir1 + file_nu1)
nu_lbl2, flux_lbl2 = np.loadtxt(oddir2 + file_nu2)
nu_lbl3, flux_lbl3 = np.loadtxt(oddir3 + file_nu3)

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

if plot_results:
    plt.figure(1)
    plt.plot(nu_lbl1, flux_lbl1, color="cyan", label="Clear sim.")
    plt.plot(nu_lbl2, flux_lbl2, color="cyan")
    plt.plot(nu_lbl3, flux_lbl3, color="cyan")

# wavenumber region 2
for i in range(10):
    ind = i * 10000
    nu_dis2, flux_dis2 = np.loadtxt(f"{oddir2}flux_disort_{ind}")
    totflux_dis2 += np.trapz(flux_dis2, nu_dis2)
    nu_dis2cld, flux_dis2cld = np.loadtxt(f"{oddir2}flux_disort_cldy_{ind}")
    totflux_dis_cld2 += np.trapz(flux_dis2cld, nu_dis2cld)
    if plot_results:
        if i == 1:
            # plt.plot(nu_dis2, flux_dis2, "k", label="DISORT, clear")
            plt.plot(nu_dis2cld, flux_dis2cld, "b", label="Cloudy sim.")
        else:
            # plt.plot(nu_dis2, flux_dis2, "k", label="DISORT, clear")
            plt.plot(nu_dis2cld, flux_dis2cld, "b")

# wavenumber region 3
for i in range(8):
    ind = i * 10000
    nu_dis3, flux_dis3 = np.loadtxt(f"{oddir3}flux_disort_{ind}")
    totflux_dis3 += np.trapz(flux_dis3, nu_dis3)
    nu_dis3cld, flux_dis3cld = np.loadtxt(f"{oddir3}flux_disort_cldy_{ind}")
    totflux_dis_cld3 += np.trapz(flux_dis3cld, nu_dis3cld)
    if plot_results:
        # plt.plot(nu_dis3, flux_dis3, 'k')
        plt.plot(nu_dis3cld, flux_dis3cld, "b")


# Compare the LBLRTM and DIOSRT clear-sky and the DISORT cloudy-sky results
print(f"LBLRTM 1: {totflux_lbl1}")
print(f"LBLRTM 2: {totflux_lbl2}")
print(f"DISORT 2: {totflux_dis2}")
print(f"CLOUDY 2: {totflux_dis_cld2}")
print(f"LBLRTM 3: {totflux_lbl3}")
print(f"DISORT 3: {totflux_dis3}")
print(f"CLOUDY 3: {totflux_dis_cld3}")

totflux_lbl_clr = totflux_lbl1 + totflux_lbl2 + totflux_lbl3
totflux_dis_clr = totflux_lbl1 + totflux_dis2 + totflux_dis3
totflux_dis_cld = totflux_lbl1 + totflux_dis_cld2 + totflux_dis_cld3
print(f"LBLRTM total flux: {totflux_lbl_clr}")
print(f"DISORT total flux: {totflux_dis_clr}")
print(f"CLOUDY total flux: {totflux_dis_cld}")
print(f"DISORT - LBLRTM: {totflux_dis_clr - totflux_lbl_clr}")
print(f"CLOUDY - LBLRTM: {totflux_dis_cld - totflux_lbl_clr}")


# Get the LWD flux from ERA5
date1 = dt.datetime(2017, 1, 31, tzinfo=pytz.utc)
date2 = dt.datetime(2017, 2, 2, tzinfo=pytz.utc)
era = get_era5_escudero_byyear(era_dir, era_fmt, date1, date2, lat, lon)

# Load in the LWD flux from King Sejong
ksj = pd.read_csv(ksj_dir + ksj_met_file, skiprows=1)
ksj_dtime = pd.to_datetime(ksj["Timestamp(LST=UTC-4)"]).to_list()
ksj_dtime = [x + dt.timedelta(hours=4) for x in ksj_dtime]

# Get LWD fluxes for blackbody clouds
prof = Prof(prof_file)
temp_trop = prof.tm[prof.zm < 10]
temp_max = float(max(temp_trop))
temp_21 = float(prof.tm[prof.zm == 2.1])
temp_23 = float(prof.tm[prof.zm == 2.3])

nutot = [0.0, 2222.22]
nuwin = []

bbflux_max_tot, nu_bb, bbflux_max = get_bb_fluxes(nutot[0], nutot[1], temp_max)
bbflux_2km_tot, nu_bb, bbflux_2km = get_bb_fluxes(nutot[0], nutot[1], temp_21)
b277, nu_bb, _ = get_bb_fluxes(nutot[0], nutot[1], 277.0)
b270, nu_bb, _ = get_bb_fluxes(nutot[0], nutot[1], 270.0)
b265, nu_bb, _ = get_bb_fluxes(nutot[0], nutot[1], 265.0)
b260, nu_bb, _ = get_bb_fluxes(nutot[0], nutot[1], 260.0)
b250, nu_bb, _ = get_bb_fluxes(nutot[0], nutot[1], 250.0)
b255, nu_bb, _ = get_bb_fluxes(nutot[0], nutot[1], 255.0)


xrtick = dt.datetime(2017, 2, 1, 0, 20)
date_range = [dt.datetime(2017, 1, 31), dt.datetime(2017, 2, 1)]

# The King Sejong data is not correct. Perhaps it is something like:
# sigma = 5.67e-8
# tpyrg = 327
# plt.plot(ksj_dtime, sigma*tpyrg**4-ksj[" Rlw"], label="Meas: KSJ")

plt.figure(2)
plt.plot(ksj_dtime, ksj[" Rlw"], label="Meas: KSJ")
plt.plot(era["date"], era["lwd"], label="ERA5")
# plt.plot(snd_date, bbflux_max_tot, "r^", label="BB(Tmax)")
plt.plot(snd_date, totflux_dis_cld, "o", label="Cloudy, sim.")
# plt.plot(snd_date, bbflux_2km_tot, "bs", label="BB(T(2.1 km))")
plt.plot(snd_date, totflux_lbl_clr, "cd", label="Clear, sim.")
plt.plot(date_range, [b277, b277], "k:")
plt.plot(date_range, [b270, b270], "k:")
plt.plot(date_range, [b265, b265], "k:")
plt.plot(date_range, [b260, b260], "k:")
plt.plot(date_range, [b255, b255], "k:")
plt.plot(date_range, [b250, b250], "k:")
plt.text(xrtick, 321, "B$_{277}$")
plt.text(xrtick, 296, "B$_{270}$")
plt.text(xrtick, 277, "B$_{265}$")
plt.text(xrtick, 255, "B$_{260}$")
plt.text(xrtick, 238, "B$_{255}$")
plt.text(xrtick, 219, "B$_{250}$")
plt.legend()
plt.xlim(date_range)
plt.ylabel("Flux (W/m$^2$)")

# Add bb fluxes to wavenumber-dep plot
if plot_results:
    plt.figure(1)
    plt.plot(nu_bb, bbflux_max, label=f"B({round(temp_max)} K)")
    plt.plot(nu_bb, bbflux_2km, label=f"B({round(temp_21)} K)")
    # plt.plot(nu_bb, bbflux_280, label="BB(280 K)")
    plt.legend()
    plt.xlabel("Wavenumber (cm$^{-1}$)")
    plt.ylabel("LWD Flux [W/(m$^2$ cm$^{-1}$)]")
    plt.xlim([0, 2223])

# Plot for SangJong
date_range = [ksj_dtime[0], ksj_dtime[-1]]
bbflux_290_tot, nu_bb, _ = get_bb_fluxes(nutot[0], nutot[1], 290.0)
bbflux_296_tot, nu_bb, _ = get_bb_fluxes(nutot[0], nutot[1], 296.0)
bbflux_280_tot, nu_bb, bbflux_280 = get_bb_fluxes(nutot[0], nutot[1], 280.0)

[b278, b278] = [bbflux_280_tot, bbflux_280_tot]
[b290, b290] = [bbflux_290_tot, bbflux_290_tot]
[b296, b296] = [bbflux_296_tot, bbflux_296_tot]

xrtick = dt.datetime(2017, 2, 2, 5)
plt.figure(3)
plt.plot(ksj_dtime, ksj[" Rlw"], label="Meas: KSJ")
plt.plot(era["date"], era["lwd"], label="ERA5")
plt.plot(snd_date, totflux_dis_cld, "o", label="Cloudy sim.")
plt.plot(snd_date, totflux_lbl_clr, "cd", label="Clear, sim.")
plt.plot(date_range, [b270, b270], "k:")
plt.plot(date_range, [b278, b278], "k:")
plt.plot(date_range, [b290, b290], "k:")
plt.plot(date_range, [b296, b296], "k:")
plt.text(xrtick, 296, "B$_{270}$")
plt.text(xrtick, 360, "B$_{283}$")
plt.text(xrtick, 397, "B$_{290}$")
plt.text(xrtick, 430, "B$_{296}$")
plt.legend()
plt.xlim(date_range)
plt.ylabel("Flux (W/m$^2$)")


# Save the results
if save_results:
    # Save the Escudero SWD data as CSV
    dfm = pd.DataFrame({"date": swd_date, "swd": swd})
    dfm.to_csv(outdir + esc_swd_file, index=False)

    # Save ERA5 data as CSV
    dfm = pd.DataFrame(era)
    dfm.to_csv(outdir + era5_bbnd_file, index=False)

    # Save computed clear and cloudy SWD flux as two-line CSV
    with open(outdir + computed_flux_file, "w", encoding="utf-8") as fid:
        # Header
        fid.write("date,")
        fid.write("flux_lwd_lblrtm_clear,")
        fid.write("flux_lwd_disort_clear,")
        fid.write("flux_lwd_disort_cloudy,")
        fid.write("planck250,planck255,planck260,planck265,")
        fid.write("planck270,planck277\n")
        # Values
        fid.write(f"{snd_date.strftime('%Y%m%d %H:%M:%S')},")
        fid.write(f"{round(totflux_lbl_clr, 3)},")
        fid.write(f"{round(totflux_dis_clr, 3)},")
        fid.write(f"{round(totflux_dis_cld, 3)},")
        fid.write(f"{rnd3(b250)},{rnd3(b255)},{rnd3(b260)},{rnd3(b265)},")
        fid.write(f"{rnd3(b270)},{rnd3(b277)}")

    # Save profile as text file
    dfm = pd.DataFrame(
        {
            "press": prof.pm,
            "temp": prof.tm,
            "alt": prof.zm,
            "h2o": prof.traceGases[:, 0],
            "co2": prof.traceGases[:, 1],
            "ozone": prof.traceGases[:, 2],
        }
    )
    dfm.to_csv(outdir + prof_file_out, index=False)
