"""
Created on Tue Oct 29 13:04:15 2019

@author: prowe

Purpose: Run DISORT using inputs specified in a namelist file

Instructions:
    1) You will need to install f90nml:
       $ pip install f90nml
    2) Install rundisort_py
    3) Change the path below as needed for your setup.
    4) Make sure you have the file "sample.nml" in your path.
    5) Results should agree with those from sampleRun.

Copyright 2019 by Penny Rowe and NorthWest Research Associates
"""


# Set the path as needed for your directory structure)
# import sys

# sys.path.append("/Users/prowe/Git_repos/rundisort_py/installation/")
# sys.path.append("/Users/prowe/Git_repos/run_radtran/sample_run/")
# sys.path.append("/Users/prowe/Git_repos/clarra_py/")
# sys.path.append("/Users/prowe/Git_repos/general_py/")  # sun_position.py


# Dependencies
import numpy as np
import os.path
import glob
import datetime as dt
import matplotlib.pyplot as plt

# Our modules
from antarc import lblrtm_utils
from antarc.radtran_inputs import Inputs
from antarc.humidRS80 import humidRS80
from antarc.plancknu import plancknu


def run_lblrtm_hires_rad(prof, dir_lblrtm, run_lblrtm, v1, v2, angle):
    """
    @param prof, dir_lblrtm, run_lblrtm, v1, v2, angl
    """
    import subprocess

    #  .. Clean out the LBLRTM directory
    for filen in glob.glob(dir_lblrtm + "ODdeflt_*"):
        os.remove(filen)
    for filen in glob.glob(dir_lblrtm + "TAPE12"):
        os.remove(filen)
    for filen in glob.glob(dir_lblrtm + "TAPE5"):
        os.remove(filen)

    # try:
    # .. Write the TAPE5 file using the new python code!
    #    prof better be the name of the profile file, including directory
    lblrtm_utils.write_tape5_prof(prof, v1, v2, angle, "radiance", dir_lblrtm)

    # ..  Run LBLRTM (slow!)
    p = subprocess.Popen([run_lblrtm], cwd=dir_lblrtm)
    p.wait()

    # .. Load in the perfect resolution file
    nu_mono, rad_mono, _ = lblrtm_utils.read_TAPE12_lblrtm(
        dir_lblrtm + "TAPE12", prec="d"
    )

    rad_mono = 1e7 * rad_mono

    return nu_mono, rad_mono


# # # # # # # #      INPUTS      # # # # # # # # # # # # #

# Parameter modules, specific to user
from antarc import params

# Parameter modules
from antarc.escudero.parameters import radtran_params


# Directories
out_dir = params.MEAS_DIR + "Escudero/lwd/"
proj_dir = params.PROJECT_DIR

# Flags
do_run_clear_sky_sims = True

# Inputs
viewAngle = 0  # do we need this?
Xi = [1e-4, 0, 15, 30]  # tau, fice, Reliq, Reice
iTempLiq = 0
wTempLiq = 1
iTempIce = 0
wTempIce = 1
iLiqLayer = [0]
iIceLayer = [0]
bwn = 238.09524
ewn = 1500.0  # 2222.22      #
get_upper_atm = False
get_Bctc_best = False
N = 16384

# Desired date range
start_date = dt.datetime(2022, 2, 5)
end_date = dt.datetime(2022, 2, 8)

# Process number (LBLRTM directory number)
proc = 1
# 1: 2018/12/4: 238 - 1500 (done)
# # # # # # # # # # # # # # # # # # # # # #


# # # # # #   MAIN CODE    # # # # # # # #
savefigs = False
lblrtm_dir_base = radtran_params.lblrtm_dir

ip = Inputs(radtran_params)
ip.prepInputs()

prof_datenums = ip.prof_datenums
prof_files = ip.prof_files
prof_dir = ip.prof_dir
run_lblrtm = ip.run_lblrtm
iprofdate = ip.iprofdate
iprofhour = ip.iprofhour

ords = np.array([x.toordinal() for x in prof_datenums])
iprofs = np.intersect1d(
    np.where(ords >= start_date.toordinal()),
    np.where(ords <= end_date.toordinal()),
)

color = ["blue", "green", "orange", "red", "purple"]
plt.figure(num=1)
plt.clf()
maxtemp = np.zeros(len(iprofs))
# Load in the profile and get the date
for iprof in iprofs:
    thisdatetime = prof_datenums[iprof]
    print(thisdatetime)
    prof_file = prof_files[iprof]
    prof = lblrtm_utils.Prof(prof_dir + prof_file)

    # Convert the ppmv to rh
    h2o = prof.traceGases[:, 0]
    rh = humidRS80(h2o, prof.tm, "ppmv", "rhw", prof.pm)
    datestr = thisdatetime.strftime("%m-%d %H")
    c = color[iprof]
    plt.subplot(121)
    plt.plot(prof.tm, prof.zm, color=c, label=datestr)
    plt.subplot(122)
    plt.plot(rh, prof.zm, color=c, label=thisdatetime)
    maxtemp[iprof] = max(prof.tm[prof.zm < 14])

plt.figure(1)
plt.subplot(121)
plt.plot([273, 273], [0, 20], ":")
plt.legend()
plt.ylim([0, 15])
plt.xlabel("Temperature (K)")
plt.ylabel("Altitude (km)")

plt.subplot(122)
plt.ylim([0, 15])
plt.xlabel("RH (%)")

out_dir = proj_dir + "NSF_AP/case_studies/Feb_2022/figures/"
prof_figname = out_dir + "prof.png"
if savefigs:
    plt.savefig(prof_figname)

nu = np.arange(238.09524, 2222.22, 0.1)
temps = np.arange(220, max(maxtemp), 10)

# Fit to a quartic
zangles = [0.0, 20.0, 40.0, 60.0, 80.0]
mu = np.cos(np.deg2rad(zangles))
flux_down_tot = np.zeros(len(temps))
for i, temp in enumerate(temps):

    bt = plancknu(nu, temp)
    rad = np.vstack([bt, bt, bt, bt, bt])
    rad = rad.T

    flux_down_opaque = np.zeros(len(nu))
    for inu in range(len(nu)):
        p = np.polyfit(mu, rad[inu, :], 4)
        flux_down_opaque[inu] = (
            1 / 6 * p[0]
            + 1 / 5 * p[1]
            + 1 / 4 * p[2]
            + 1 / 3 * p[3]
            + 1 / 2 * p[4]
        )
    flux_down_opaque = 2 * np.pi * flux_down_opaque
    flux_down_tot[i] = np.trapz(flux_down_opaque, nu)

plt.figure(num=2)
plt.plot(temps, flux_down_tot, ".-")
plt.xlabel("Temperature (K)")
plt.ylabel("Max Cloud LWD Flux (W/m$^{2}$)")

figname = out_dir + "max_cloud_lwd.png"
if savefigs:
    plt.savefig(figname)


#     # Set some variables
#     nlyr_trop = len(prof.tm) - 1
#     nlyr_toa = nlyr_trop
#     cos_angle = np.cos(viewAngle)
#     layerCloud = np.array([nlyr_toa - 1])  # surface (clear sky anyway), list

#     lblrtm_dir = lblrtm_dir_base + str(proc) + "/"

#     for i, zangle in enumerate(zangles):
#         nu_lbl, rad0 = run_lblrtm_hires_rad(
#             prof, lblrtm_dir, lblrtm_dir + run_lblrtm, bwn, ewn, zangle
#         )
#         if i == 0:
#             rad = np.zeros((len(nu_lbl), len(zangles)))

#         rad[:, i] = rad0


#     plt.figure(1)
#     plt.cla()
#     plt.plot(nu_lbl, 1e-3 * flux_down_lbl)

#     # Compare total flux (integrated over wavenumber)
#     flux_down_tot_lbl = np.trapz(1e-3 * flux_down_lbl, nu_lbl)
#     print("flux from lblrtm:", flux_down_tot_lbl)

#     # .. Save these results
#     outfile = (
#         out_dir
#         + "flux"
#         + prof_file[iprofdate[0][0] : iprofhour[0][1]]
#         + "_"
#         + str(int(bwn))
#         + "_"
#         + str(int(ewn))
#         + ".txt"
#     )
#     f = open(outfile, "w")
#     f.write(str(flux_down_tot_lbl))
#     f.close()

#     """
#     # .. Verify that ozone is ok by looking at low-res spectrum
#     from reduce_resolution import reduce_resolution
#     nu_lo, rad_lo = reduce_resolution(nu_lbl, rad0, .5, 16384, 2**22)
#     rad_lo = np.real(rad_lo)
#     # .. The ozone looks a little low, but leave for now ...
#     """

#     """


#     # .. Run LBLRTM to get the optical depths
#     lblrtm_utils.do_lblrtm_run(lblrtm_dir, lblrtm_dir + ip.run_lblrtm,
#                                prof, ip.v1, ip.v2, lblrtm_dir)


#     # .. Get low-resolution optical depths for DISORT (slow)
#     nu, rads, tc, _, \
#     _, _, _, _, _ = radtran_lblrtm_od.radtran_lblrtm_od( \
#                 lblrtm_dir, nlyr_trop, nlyr_toa, prof.tm, 0, \
#                 ip.dnu, bwn, ewn, N, prof.zm, cos_angle, 'd', \
#                 get_upper_atm, get_Bctc_best)

#     # # # # # #

#     # .. Get the layer optical depths corresponding to the
#     #    transmittances and get the radiances above
#     Nnu = len(nu)
#     odtot_Lm1 = np.zeros(Nnu)
#     dtau_gas = np.zeros((Nnu, nlyr_toa))
#     for ilayer in range(nlyr_toa):
#         # .. Layer effective-resolution optical depth
#         trans0 = copy.deepcopy(tc[:,ilayer])
#         trans0[trans0>=1] = 1
#         trans0[trans0<=0] = 1e-40
#         odtot = -np.log(trans0)
#         dtau_gas[:,ilayer] = odtot - odtot_Lm1

#         dtau_gas[dtau_gas[:,ilayer]<0,ilayer] = 1e-8

#         odtot_Lm1 = odtot

#     # # # # # #
#     # .. TOA to surface temperature for DISORT
#     TEMPER = np.flipud(prof.tm)

#     delv = 0*nu + 1

#     # .. Get pmoms, even though this is clear sky
#     sspLiq = get_ssp([ip.pmomfiles_liq], nu)            # SSP liquid
#     sspIce = sspLiq                                     # SSP ice - use liquid


#     # .. Solar spectrum
#     FBEAM = getsolarbeam_IR(nu, ip.solarSourceFunFile)

#     # .. Solar angles
#     sun = sun_position(thisdatetime, ip.lat, ip.lon, ip.alt)
#     sza = sun.zenith

#     # .. UMU: cosine of scene viewing angle, defined from TOA
#     #         -1=>looking up (downwelling), 1=>looking down (upwelling)
#     viewAngle_TOA = 180 - viewAngle

#     # Assign values
#     saa = sun.azimuth
#     UMU0 = np.cos(np.deg2rad(sza[0]))
#     UMU = [np.cos(np.deg2rad(viewAngle_TOA))]

#     # .. Surface albedo
#     nu_surf_albedo, surf_albedo = \
#         get_surface_albedo_from_file(ip.surfEmissDataFile)

#     albedo = np.interp(nu, nu_surf_albedo, surf_albedo)
#     if any([any(albedo < 0), any(albedo > 1), \
#        any(np.isnan(albedo))]):
#            raise NameError('bad albedo')


#     cld_tau_liq_tot = Xi[0] * (1-Xi[1])
#     cld_tau_ice_tot = Xi[0] * Xi[1]
#     weight = 1. * np.ones((1,1))
#     Rc, izm, flux_down, flux_up = run_disort_flux(\
#       nu, np.flipud(dtau_gas.T),
#       layerCloud,
#       Xi[2], cld_tau_liq_tot,
#       Xi[3], cld_tau_ice_tot,
#       iTempLiq, weight * cld_tau_liq_tot,
#       iTempIce, weight * cld_tau_ice_tot,
#       iLiqLayer, iIceLayer, sspLiq, sspIce,
#       TEMPER, ip.NSTR,
#       UMU0, UMU,
#       albedo, FBEAM,
#       delv, len(TEMPER))


#     plt.figure(1); plt.cla()
#     plt.plot(nu_lbl, 1e-3*flux_down_lbl, label='LBLRTM down flux')
#     plt.plot(nu, flux_down, label='DISORT down flux')
#     plt.plot(nu, flux_up, label='DISORT up flux')
#     plt.legend()

#     # Compare total flux (integrated over wavenumber)
#     flux_down_tot_lbl = np.trapz(1e-3*flux_down_lbl, nu_lbl)
#     flux_down_tot_dis = np.trapz(flux_down, nu)
#     print('flux from lblrtm:', flux_down_tot_lbl)
#     print('flux from DISORT:', flux_down_tot_dis)
# """
