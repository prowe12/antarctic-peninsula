#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 13:35:07 2023

@author: prowe

Purpose: Get gaseous optical depths for radiosonde launches of 2017/01/31



"""

# Dependencies
import datetime as dt
import numpy as np
import pytz
import time


# Parameter modules specific to the user: the user must create params.py
# in the directory antarc with the variables MEAS_DIR and PROJ_DIR set
from antarc import params

# Parameter modules
from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import radtran_params

# Our modules
from antarc import lblrtm_utils
from antarc.ssp import get_ssp
from antarc.load_ssp_nu import load_ssp_nu
from antarc.get_surface_albedo_from_file import get_surface_albedo_from_file
from antarc.get_solar_beam import get_solar_beam
from antarc.sun_position import sun_position
from antarc.radtran_inputs import Inputs
from antarc.run_radtran import run_lblrtm_hires_od, run_lblrtm_hires_rad
from antarc.lblrtm_utils import read_od_lblrtm
from antarc.escudero.case201701.run_disort import run_disort


def set_sun_angle_vars(lat:float, lon:float, alt:float, dtime):
    sun = sun_position(dtime, lat, lon, alt)
    sza = sun.zenith
    return np.cos(np.deg2rad(sza[0]))
    
def get_clear_flux(
    prof, 
    outdir: str,
    lblrtm_dir_base: str,
    proc: int,
    bwn: float,
    ewn: float,
    file_flux_tot: str,
    file_flux_nu: str,
):
    """
    Get clear sky fluxes for input beginning and ending wavenumber
    @param outdir  Directory to save results
    @param params  Various necessary parameters (see within)
    @param dates
    @param proc  Process number (LBLRTM directory number)
    @param redo  Whether to redo the calculation if it exists
    @param bwn  Beginning wavenumber
    @param ewn  Ending wavenumber
    """

    lblrtm_dir = lblrtm_dir_base + "run" + str(proc) + "/"
    run_cmd = lblrtm_dir + "lblrtm_dbl"        

    zangles = [0.0, 20.0, 40.0, 60.0, 80.0]
    for i, zangle in enumerate(zangles):
        nu_lbl, rad0 = run_lblrtm_hires_rad(
            prof, lblrtm_dir, run_cmd, bwn, ewn, zangle
        )
        if i == 0:
            rad = np.zeros((len(nu_lbl), len(zangles)))
            nu_prime = nu_lbl
        elif (
              len(nu_lbl) != len(nu_prime) 
              or not np.isclose(nu_lbl[0], nu_prime[0])
              or not np.isclose(nu_lbl[-1], nu_prime[-1])
             ):
            rad0 = np.interp(nu_prime, nu_lbl, rad0)

        rad[:, i] = rad0
        
    nu_lbl = nu_prime

    # Fit to a quartic
    cos_zang = np.cos(np.deg2rad(zangles))
    nlbl = len(nu_lbl)
    flux_down_lbl = np.zeros(nlbl)
    for inu in range(nlbl):
        p = np.polyfit(cos_zang, rad[inu, :], 4)
        flux_down_lbl[inu] = (
            1 / 6 * p[0]
            + 1 / 5 * p[1]
            + 1 / 4 * p[2]
            + 1 / 3 * p[3]
            + 1 / 2 * p[4]
        )
    flux_down_lbl = 2 * np.pi * flux_down_lbl
    # Compute total flux (integrated over wavenumber)
    flux_down_tot_lbl = np.trapz(1e-3 * flux_down_lbl, nu_lbl)

    # Save the wavenumber-dependent partial fluxes
    np.savetxt(outdir + file_flux_nu, (nu_lbl, 1e-3*flux_down_lbl))
    
    # Save the total flux
    with open(outdir + file_flux_tot, "w", encoding="utf-8") as fid:
        fid.write(str(flux_down_tot_lbl))
            
            
            
def get_clear_od_wn(
    prof,
    lblrtm_dir_base: str,
    outdir: str,
    proc: int,
    bwn: float,
    ewn: float,
):
    """
    Get clear sky fluxes for input beginning and ending wavenumber
    @param outdir  Directory to save results
    @param proc  Process number (LBLRTM directory number)
    @param bwn  Beginning wavenumber
    @param ewn  Ending wavenumber
    """
    # # # # # #   MAIN CODE    # # # # # # # #
    lblrtm_dir = lblrtm_dir_base + "run" + str(proc) + "/"
    run_cmd = lblrtm_dir + "lblrtm_dbl"        
    zangle = 0
    nu_lbl, rad = run_lblrtm_hires_od(
        prof, lblrtm_dir, run_cmd, bwn, ewn, zangle, outdir,
    )

def qc_inputs(total_od_vis_liq, dtau_gas, iobs, temper, nssp, nnu):
    # Make sure dims are all ok
    if not type(reff_liq) == float:
        raise ValueError('reff_liq should be float')
    if not type(total_od_vis_liq) == np.float64:
        raise ValueError('total_od_vis_liq should be float')
    if len(cldlyr.shape) != 1 or len(iliq_layer.shape) != 1:
        raise ValueError('cldlyr and iliq_layer should be vectors')
    if len(cldlyr) < len(iliq_layer):
        raise ValueError('num cldlyrs must be >= num liq_layers')
    ncld = len(cldlyr)
    if itemp_liq.shape != (ncld, nssp) or od_vis_liq.shape != (ncld, nssp):
        raise ValueError('itemp_liq and od_vis_liq should be ncld x nssp')
    
    if dtau_gas.shape != (nlayers, nnu):
        raise ValueError('Bad shape for dtau_gas')
    if iobs < 0 or iobs > nlayers:
        raise ValueError('iobs outside model atmosphere')
    if len(temper.shape) != 1 or len(temper) != nlayers + 1:
        raise ValueError('Bad shape for temper')

def get_flux_down(outdir, outfile_pfx):
    wnum0, taulay = read_od_lblrtm(outdir + 'ODdeflt_001', 'd')    
    tau = np.zeros([nlayers, len(wnum0)])
    tau[0, :] = taulay
    
    for i in range(2, nlayers + 1):
        fname = 'ODdeflt_0' + str(i).zfill(2)
        wnum, taulay = read_od_lblrtm(outdir + fname, 'd')
        if len(wnum0) == len(wnum):
            tau[i-1, :] = taulay
        else:
            taulay = np.interp(wnum0, wnum, taulay)
            tau[i-1, :] = taulay
            
    # Flip temp, dtau_gas so they go from TOA down
    dtau_gas = np.flipud(tau)  # nlyr x nnu
    temper = np.flipud(prof.tm[:nlayers+1])
    iobs = np.shape(dtau_gas)[0]
    wnum = wnum0
        
    inds_list = []
    for i in range(int(np.ceil(len(wnum)/10000))):
        inds_list.append([i*10000, (i+1)*10000])
            
    for ind in inds_list:
        wnums = wnum[ind[0]:ind[1]]
        
        ssp_liq = get_ssp(pmom_files, wnums)
        
        nstr = 16
        delv_nu = 1.
        
        surf_albedo = get_surface_albedo_from_file(surf_albedo_file)
        albedo_nu = np.interp(wnums, surf_albedo[:,0], surf_albedo[:,1])
        fbeam_nu = get_solar_beam(wnums, beam_file)
        
        cos_in_beam_angle = set_sun_angle_vars(lat, lon, prof.zm[0], dtime)
        cos_out_angles = [np.cos(np.deg2rad(180 - 0))]
        
                
        # No ice in cloud
        reff_ice = 20.
        total_od_vis_ice = np.array([])
        itemp_ice = np.array([[]])
        od_vis_ice = np.array([[]])
        iice_layer = np.array([[]])
        
        total_od_vis_liq = np.sum(od_vis_liq)


        qc_inputs(total_od_vis_liq, dtau_gas, iobs, temper, len(pmom_files), 
                  len(wnum))
        
        rad, izm, flux_down, flux_up, warn = run_disort(
            wnums, dtau_gas[:, ind[0]:ind[1]], cldlyr,
            reff_liq, total_od_vis_liq,
            reff_ice, total_od_vis_ice,
            itemp_liq, od_vis_liq, itemp_ice, od_vis_ice,
            iliq_layer, iice_layer,
            ssp_liq, ssp_liq, temper, nstr, cos_in_beam_angle, cos_out_angles,
            albedo_nu, fbeam_nu, delv_nu, iobs)
        
        
        # Save the wavenumber-dependent partial fluxes
        np.savetxt(f"{outdir}{outfile_pfx}{ind[0]}", (wnums, flux_down))
    
        elapsed = time.time() - timey
        print(f"Elapsed time: {round(elapsed)}s")
        

timey = time.time()

# Directories
MEAS_DIR = params.MEAS_DIR
esc_dir = MEAS_DIR + "Escudero/"
outdir1 = esc_dir + "od/od_20170131_18a/"
outdir2 = esc_dir + "od/od_20170131_18b/"
outdir3 = esc_dir + "od/od_20170131_18c/"

# Files
prof_file = esc_dir + "profiles_cloud/prof20170131_1801.nc"
surf_albedo_file = esc_dir + "surface_emissivity/soil_emissivity_concord1.txt"
beam_file = esc_dir + "kurucz.dat"
pmomfile_liq_253 = esc_dir + "ssp/pmom_water_T253_RFN_S331.nc" # 0
pmomfile_liq_263 = esc_dir + "ssp/pmom_water_T263_RFN_S331.nc" # 1
# pmomfile_liq_273 = esc_dir + "ssp/pmom_water_T273_RFN_S331.nc" # 2

dtime = dt.datetime(2017, 1, 31, 18, 1, 0)

# Params needed from param files
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
lbldir = params.LBLRTM_DIR

find_wnum_cutoffs = False


# Within this range the below-cloud transmittance is <= 1e-6
# so we can just use the clear-sky fluxes
bwn1 = 0.0
ewn1 = 360.0
# Within this range the below-cloud transmittance is > 1e-6
# so we need to do the DISORT calculation
bwn2 = 360.0
ewn2 = 1484.7
# Need to test this
bwn3 = 1484.7
ewn3 = 2222.2  # pyrgeometer cut-off

# Get the atmospheric profile
prof = lblrtm_utils.Prof(prof_file)

# List of pmom files that will be needed
pmom_files = [pmomfile_liq_253, pmomfile_liq_263]

# Compute the optical depths and save them to the respective directories
# get_clear_od_wn(prof, lbldir, outdir1, 1, bwn1, ewn1)
# get_clear_od_wn(prof, lbldir, outdir1, 2, bwn2, ewn2)
# get_clear_od_wn(prof, lbldir, outdir1, 3, bwn3, ewn3)

# Compute the clear-sky fluxes as functions of wavenumber
file_tot1 = f"flux_tot_{int(bwn1)}_{int(ewn1)}.txt"
file_tot2 = f"flux_tot_{int(bwn2)}_{int(ewn2)}.txt"
file_tot3 = f"flux_tot_{int(bwn3)}_{int(ewn3)}.txt"
file_nu1 = f"flux_nu_{int(bwn1)}_{int(ewn1)}.txt"
file_nu2 = f"flux_nu_{int(bwn2)}_{int(ewn2)}.txt"
file_nu3 = f"flux_nu_{int(bwn3)}_{int(ewn3)}.txt"

# get_clear_flux(prof, outdir1, lbldir, 1, bwn1, ewn1, file_tot1, file_nu1)
# get_clear_flux(prof, outdir2, lbldir, 2, bwn2, ewn2, file_tot2, file_nu2)
# get_clear_flux(prof, outdir3, lbldir, 3, bwn3, ewn3, file_tot3, file_nu3)

# Compute the clear-sky downwelling fluxes
nlayers = len(prof.zm) - 1

# Find wavenumber cut-offs
if find_wnum_cutoffs:
    wnum0, taulay = read_od_lblrtm(outdir + 'ODdeflt_001', 'd')    
    tau = np.zeros([nlayers, len(wnum0)])
    tau[0, :] = taulay
    tautot = taulay
    
    # optical depth where transmittance < .001
    ihi_od_set = set()
    hi_od = -np.log(.00001)
    itop = (nlayers * np.ones(len(wnum0))).astype(int)
    for i in range(2, nlayers + 1):
        fname = 'ODdeflt_0' + str(i).zfill(2)
        wnum, taulay = read_od_lblrtm(outdir + fname, 'd')
        if len(wnum0) == len(wnum):
            tau[i-1, :] = taulay
            tautot += taulay
        else:
            taulay = np.interp(wnum0, wnum, taulay)
            tau[i-1, :] = taulay
            tautot += taulay
        # Find where the transmittmance is less than .001, or
        # where the total optical depth is greater than ~7
        ihi_od_set = set(np.where(tautot > hi_od)[0]) - ihi_od_set
        ihi_od = list(ihi_od_set)
        itop[ihi_od] = np.min([itop[ihi_od], np.ones(len(ihi_od))*(i-1)], axis=0)
            
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(wnum0, np.exp(-tautot))
    ido = np.where(np.exp(-tautot) >= 1e-6)[0]


# Get clear-sky fluxes
# Liquid part of cloud
cldlyr = np.array([44])
reff_liq = 10.
itemp_liq = np.array([[0]])
od_vis_liq = np.array([[1e-4]])  #  (layers x ssps)
iliq_layer = np.array([0])        

# get_flux_down(outdir2, "flux_disort_")
# get_flux_down(outdir3, "flux_disort_")


# Get cloudy-sky fluxes

# The cloud data
# alt bounds  optical depth  effective radius
cloud_data = np.array([
    [ 2.10,  2.15,  0.47963,    3.86],
    [ 2.15,  2.20,  0.2693,     4.28],
    [ 2.20,  2.25,  1.0172,     4.09],
    [ 2.25,  2.30,  2.412,      5.86],
    [ 2.30,  2.35,  5.046956,   5.72],
    [ 2.35,  2.40,  2.81574,    6.28],
])
# [ 2050  2100  8.8715e-5  2.65  omit
# [ 2400  2450  0.0044285  5.52  omit

# Get mean temperature for each cloud layer
itemp_liq = np.zeros([len(cloud_data), 2]).astype(int)
itemp_liq[:,1] = 1
od_vis_liq = np.zeros([len(cloud_data), 2])
# meanalt = []
# meantemp = []
for i, dat in enumerate(cloud_data):
    ily1 = np.where(prof.zm >= dat[0])[0]
    ily2 = np.where(prof.zm <= dat[1])[0]
    ily = np.intersect1d(ily1, ily2)
    meantemp0 = np.mean(prof.tm[ily])
    # meanalt.append(np.mean(prof.zm[ily]))
    # meantemp.append(meantemp0)
    if meantemp0 < 253 or meantemp0 > 263:
        raise ValueError('Not expecting this case: modify code')
    wt0 = (263 - meantemp0) / (263 - 253)
    wt1 = (meantemp0 - 253) / (263 - 253)
    if wt0 + wt1 != 1 or wt0 < 0 or wt1 < 0:
        raise ValueError('Weights for ssp temps are wrong')
    od_vis_liq[i,:] = [wt0 * dat[2], wt1 * dat[2]]

# Liquid part of cloud
alt_from_top = np.flipud(prof.zm)
icldtop = np.where(alt_from_top == 2.4)[0][0]
icldbot = np.where(alt_from_top == 2.15)[0][0]
cldlyr = np.arange(icldtop, icldbot+1)
reff_liq = float(np.mean(cloud_data[:,3]))  # TODO: Use layer-dependent value?
iliq_layer = np.arange(0, 6).astype(int)    

# get_flux_down(outdir2, "flux_disort_cldy_")
get_flux_down(outdir3, "flux_disort_cldy_")
    