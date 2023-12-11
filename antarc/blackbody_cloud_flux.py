#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 10:05:02 2023

@author: prowe

Get the flux for a blackbody cloud with a given temperature

Copyright 2023 by Penny Rowe and NorthWest Research Associates
"""

# Dependencies
import numpy as np

# Our modules
from antarc.plancknu import plancknu


def get_bb_fluxes(bwn: float, ewn: float, temp: float):
    """
    Get fluxes between wavenumbers for a blackbody with the specified
    temperature
    @param bwn  Beginning wavenumber
    @param ewn  Ending wavenumber
    @param temp  Desired temperature
    @returns  The total flux
    @returns  A numpy array of wavenumbers
    @returns  The wavenumber-dependent flux
    """

    zangles = [0.0, 20.0, 40.0, 60.0, 80.0]
    wnum = np.arange(bwn, ewn, 0.1)
    nwnum = len(wnum)
    rad = np.zeros([nwnum, len(zangles)])

    for i in range(len(zangles)):
        rad[:, i] = plancknu(wnum, temp)

    # Fit to a quartic
    mu = np.cos(np.deg2rad(zangles))

    flux_down_wnum = np.zeros(nwnum)
    for inu in range(nwnum):
        p = np.polyfit(mu, rad[inu, :], 4)
        flux_down_wnum[inu] = (
            1 / 6 * p[0]
            + 1 / 5 * p[1]
            + 1 / 4 * p[2]
            + 1 / 3 * p[3]
            + 1 / 2 * p[4]
        )
    # Longwave downard (LWD) flux
    flux_lwd = 2 * np.pi * flux_down_wnum
    flux_tot = np.trapz(flux_lwd, wnum)

    return flux_tot, wnum, flux_lwd
