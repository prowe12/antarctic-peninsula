#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:48:55 2022

@author: prowe

Purpose: run libradtran to compute clear sky downwelling shortwave
"""

import os
import numpy as np
import subprocess
import io


def run_libradtran(doy, sza, albedo, ozone):
    """
    @param doy  day of year
    @param sza  Solar zenight angle
    @param albedo  Surface albedo
    @param ozone  Total column ozone amount ()
    """
    libradtran_dir = "/Users/prowe/Documents/libRadtran-2.0.4/"
    cwd = os.getcwd()
    os.chdir(os.path.join(libradtran_dir, "bin"))

    inputstr = "                         # Loc of atmospheric profile file.\n"
    inputstr += "atmosphere_file ../data/atmmod/afglss.dat\n"
    inputstr += "                         # Loc of extraterrestrial spectrum\n"
    inputstr += "source solar ../data/solar_flux/gg6.txt\n"
    inputstr += f"mol_modify O3 {ozone} DU    # Set ozone column\n"
    inputstr += f"day_of_year {doy}          # Correct for Earth-Sun dist\n"
    inputstr += f"albedo {albedo}               # Surface albedo\n"
    inputstr += f"sza {sza}                 # Solar zenith angle\n"
    inputstr += "rte_solver disort        # Radiative transfer eqn solver\n"
    inputstr += "number_of_streams  6     # Number of streams\n"
    inputstr += "wavelength 280.0 2500.0  # Wavelength range [nm]\n"
    inputstr += "output_process integrate        # integrate over wavelength\n"
    inputstr += "#quiet\n"
    inputstr += "aerosol_angstrom 0 0\n"
    inputstr += "aerosol_season                                1  # summer\n"
    inputstr += "aerosol_modify ssa scale 0.99\n"
    inputstr += "altitude 0\n"
    inputstr += "aerosol_default\n"
    inputstr += "output_user  edir  eglo edn eup   enet   esum\n"

    process = subprocess.run(
        [os.getcwd() + "/uvspec"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        input=inputstr,
        encoding="ascii",
    )
    os.chdir(cwd)
    result = np.genfromtxt(io.StringIO(process.stdout))

    # Units are mW/m2 after integrating over nm, so divide by 1000 for W/m2
    # dir  eglo edn eup   enet   esum
    # edir: direct irradiance
    # eglo: global irradiance
    # edn: diffuse downward irradiance
    # eup: diffuse upward irradiance
    # enet: global - upward irradiance
    # esum: global + upward irradiance
    return result / 1000


# Example:
# dayofyear = 170
# sza = 32.0
# albedo = 0.2
# ozone = 257.0

# result = run_libradtran(dayofyear, sza, albedo, ozone)
# print(result)
