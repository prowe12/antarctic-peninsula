#  ... This file contains inputs from the namelist file
#      "disortinput.nml" in order to perform a sample run
#      of DISORT using Python.
#
#      The results should be comparable to the results
#      from running DISORT using matlab and python

# You may need to set the path as follows (change the beginning
# as needed for your directory structure):
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

Copyright 2019-2023 by Penny Rowe and NorthWest Research Associates
"""

# Dependencies
import numpy as np
import os.path
import glob
import datetime as dt
import matplotlib.pyplot as plt

# Our modules
from antarc import lblrtm_utils
from antarc.radtran_inputs import Inputs


def run_lblrtm_hires_rad(prof, dir_lblrtm, run_cmd, v1, v2, angle):
    """
    @param prof
    @param dir_lblrtm
    @param run_cmd
    @param v1  first
    @param v2
    @param angl
    """
    import subprocess

    #  .. Clean out the LBLRTM directory
    for filen in glob.glob(dir_lblrtm + "ODdeflt_*"):
        os.remove(filen)
    for filen in glob.glob(dir_lblrtm + "TAPE12"):
        os.remove(filen)
    for filen in glob.glob(dir_lblrtm + "TAPE5"):
        os.remove(filen)

    # Write the TAPE5 file
    lblrtm_utils.write_tape5_prof(prof, v1, v2, angle, "radiance", dir_lblrtm)

    # Run LBLRTM (slow!)
    p = subprocess.Popen([run_cmd], cwd=dir_lblrtm)
    p.wait()

    # Load in the perfect-resolution file
    fname = dir_lblrtm + "TAPE12"
    nu_mono, rad_mono, _ = lblrtm_utils.read_TAPE12_lblrtm(fname, prec="d")

    # except:
    #     #    For unknown reasons, for one file, the following occurs:
    #     #       v2 = 2222.000 works
    #     #       v2 = 2222.200 works
    #     #       v2 = 2222.210 works
    #     #       v2 = 2222.214 works
    #     #       v2 = 2222.216 works
    #     #       v2 = 2222.219 does not work
    #     #       v2 = 2222.220 does not work
    #     #       v2 = 2222.221 does not work
    #     #       v2 = 2222.222 does not work
    #     #       v2 = 2222.223 works!
    #     #       v2 = 2222.300 works
    #     #       v1 = 2000, v2 = 2223 ???
    #     # .. Write the TAPE5 file using the new python code!
    #     #    prof better be the name of the profile file, including directory
    #     if v2 < 2222:
    #         raise NameError("Source of problem unknown.")
    #     print("Ending wavenumber:", str(v2), "mysteriously does not work.")
    #     print("Switching to 2222.223 and retrying.")
    #     lblrtm_utils.write_tape5_prof(
    #         prof, v1, 2222.223, angle, "radiance", dir_lblrtm
    #     )

    #     # ..  Run LBLRTM (slow!)
    #     p = subprocess.Popen([run_cmd], cwd=dir_lblrtm)
    #     p.wait()

    #     # .. Load in the perfect resolution file
    #     nu_mono, rad_mono, _ = lblrtm_utils.read_TAPE12_lblrtm(
    #         dir_lblrtm + "TAPE12", prec="d"
    #     )

    rad_mono = 1e7 * rad_mono

    return nu_mono, rad_mono


def get_clear_fluxes(outdir: str, params, dates, proc: int, redo: bool = True):
    """
    Get clear sky fluxes for two sets of wavenumbers
    @param outdir  Directory to save results
    @param params  Various necessary parameters (see within)
    @param dates
    @param proc  Process number (LBLRTM directory number)
    @param redo  Whether to redo the calculation if it exists
    """
    bwn1 = 238.09524
    ewn1 = 1500.0
    bwn2 = 1500.00
    ewn2 = 2222.22

    get_clear_fluxes_wn(outdir, params, dates, proc, redo, bwn1, ewn1)
    get_clear_fluxes_wn(outdir, params, dates, proc, redo, bwn2, ewn2)


def get_clear_fluxes_wn(
    outdir: str,
    params,
    dates,
    proc: int,
    redo: bool,
    bwn: float,
    ewn: float,
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

    # # # # # #   MAIN CODE    # # # # # # # #
    plot_figs = True
    verbose = True
    lblrtm_dir_base = params.lblrtm_dir

    start_date, end_date = dates

    ip = Inputs(params)
    ip.prepInputs()

    prof_datenums = ip.prof_datenums
    prof_files = ip.prof_files
    prof_dir = ip.prof_dir
    run_cmd = ip.run_lblrtm
    iprofdate = ip.iprofdate
    iprofhour = ip.iprofhour

    ords = np.array([x.toordinal() for x in prof_datenums])
    iprofs = np.intersect1d(
        np.where(ords >= start_date.toordinal()),
        np.where(ords <= end_date.toordinal()),
    )

    # Load in the profile and get the date
    lblrtm_dir = lblrtm_dir_base + str(proc) + "/"
    zangles = [0.0, 20.0, 40.0, 60.0, 80.0]
    for iprof in iprofs:
        # thisdatetime = prof_datenums[iprof]
        prof_file = prof_files[iprof]
        prof = lblrtm_utils.Prof(prof_dir + prof_file)

        # Output file for this profile
        pfile = prof_file[iprofdate[0][0] : iprofhour[0][1]]
        outfile = f"{outdir}flux{pfile}_{int(bwn)}_{int(ewn)}.txt"

        # If we are not redoing this case, check if it exists
        # and if so, skip it and move on
        if not redo:
            if os.path.exists(outfile):
                continue

        for i, zangle in enumerate(zangles):
            nu_lbl, rad0 = run_lblrtm_hires_rad(
                prof, lblrtm_dir, lblrtm_dir + run_cmd, bwn, ewn, zangle
            )
            if i == 0:
                rad = np.zeros((len(nu_lbl), len(zangles)))

            rad[:, i] = rad0

        # Fit to a quartic
        mu = np.cos(np.deg2rad(zangles))

        nlbl = len(nu_lbl)
        flux_down_lbl = np.zeros(nlbl)
        for inu in range(nlbl):
            p = np.polyfit(mu, rad[inu, :], 4)
            flux_down_lbl[inu] = (
                1 / 6 * p[0]
                + 1 / 5 * p[1]
                + 1 / 4 * p[2]
                + 1 / 3 * p[3]
                + 1 / 2 * p[4]
            )
        # For comparison: np.trapz(mu,mu*rad[0,:])
        flux_down_lbl = 2 * np.pi * flux_down_lbl

        if plot_figs:
            plt.figure(1)
            plt.cla()
            plt.plot(nu_lbl, 1e-3 * flux_down_lbl)

        # Compare total flux (integrated over wavenumber)
        flux_down_tot_lbl = np.trapz(1e-3 * flux_down_lbl, nu_lbl)
        if verbose:
            print("flux from lblrtm:", flux_down_tot_lbl)

        # Save these results
        with open(outfile, "w", encoding="utf-8") as fid:
            fid.write(str(flux_down_tot_lbl))


if __name__ == "__main__":
    # Parameter modules
    from antarc.escudero.parameters import radtran_params

    # Directories
    out_dir = "/Users/prowe/sync/measurements/Escudero/lwd_clear/"

    # Desired date range
    start_date = dt.datetime(2022, 2, 5)
    end_date = dt.datetime(2022, 2, 8)

    proc = 1

    get_clear_fluxes(out_dir, radtran_params, (start_date, end_date), 1, False)
