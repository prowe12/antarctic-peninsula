#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 09:55:46 2022

@author: prowe
"""

# Built in modules
from scipy.interpolate import interp1d

# CLARRA Modules
from antarc.load_ssp_nu import load_ssp_nu


class SSP:
    """Class for holding single-scattering parameters and functions for
    doing interpolations"""

    def __init__(self, pmomfile, wnum):
        """Get single-scattering parameters from file
        @param pmomfile  Legendre moment file with single scattering params
        @param wnum  Wavenumbers
        """
        # Load single-scattering parameters
        self.npmom, self.pmom, self.reff, self.w0, self.qext = load_ssp_nu(
            pmomfile, wnum
        )

        # Create functions for doing interpolations
        self.fun_npmom = interp1d(self.reff, self.npmom.T)
        self.fun_qext = interp1d(self.reff, self.qext.T)
        self.fun_w0 = interp1d(self.reff, self.w0.T)


def get_ssp(pmomfiles, wnum):
    """
    Get single scattering parameters for wavenumber from pmom file
    @param pmomfiles  Legendre moment files with single scattering params
    @param wnum  Wavenumbers
    """
    nssp = len(pmomfiles)
    ssp = {}
    for i in range(nssp):
        ssp[i] = SSP(pmomfiles[i], wnum)

    return ssp
