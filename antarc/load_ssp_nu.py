#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  2 14:40:13 2022

@author: prowe
"""

from scipy.io import netcdf_file
from scipy.interpolate import interp1d
import numpy as np
import bisect


def load_ssp_nu(datafile, nu):
    """
    Return a tuple with the single-scattering parameter data
    @param datafile The filename of the file with the ssp data
    @param nu The wavenumbers of interest
    @return n_pmom  The number of legendre moments
    @return pmom  The legendre moments
    @return reff  The effective radii
    @return w0  The single-scattering albedo
    @return qext  The extinction efficiency Qext

    Default mmap is True, meaning variables are mapped to the netcdf file.
    My understanding is that using mmap = True is an advantage when
    "a file is to be opened and only a small portion of its data is to be
    read and/or written. In this scenario, MMAP will cause only the
    accessed data to be retrieved from disk." Without MMAP the whole file
    will be read into memory. "Thus, MMAP will provide some performance
    improvement in this case."
    However, when mmap = True the data must be copied or the netcdf file
    will not be closed.
    Since the ssp files are very large, we will keep mmap = True
    """
    with netcdf_file(datafile, "r") as nc:

        # Make sure chosen wavenumbers are within range of existing
        if nu[0] < nc.variables["wnum_list"][0]:
            raise ValueError(
                "ssp files do not include lowest selected wavenumbers"
            )

        if nu[-1] > nc.variables["wnum_list"][-1]:
            raise ValueError(
                "ssp files do not include highest selected wavenumbers"
            )

        # Pad out inu to low and high end by one
        inu1 = bisect.bisect(nc.variables["wnum_list"][:], nu[0])
        inu2 = bisect.bisect(nc.variables["wnum_list"][:], nu[-1])
        inu1 = max(inu1 - 1, 0)
        inu2 = min(inu2 + 1, len(nc.variables["wnum_list"][:]))
        inu = np.arange(inu1, inu2)

        n_pmomarray = nc.variables["Npmomarray"][:, inu].astype("int32")
        w0_mesh = nc.variables["w0_mesh"][:, inu].astype("float64")
        qext_mesh = nc.variables["qext_mesh"][:, inu].astype("float64")
        reff_list = nc.variables["reff_list"][:].astype("float64")
        reff = reff_list[:, 0].copy()
        wnum_vec = nc.variables["wnum_list"][inu].astype("float64")[:, 0]

        # Set up empty output arrays
        n_nu = nu.size
        n_reff = reff.size

        # Interpolate qext, w0, get an interpolated number of moments!
        f_q = interp1d(wnum_vec, qext_mesh, axis=1)
        f_w = interp1d(wnum_vec, w0_mesh, axis=1)
        f_npmom = interp1d(wnum_vec, n_pmomarray, axis=1)

        qext = f_q(nu)
        w0 = f_w(nu)
        n_pmom_dec = f_npmom(nu)

        # Use floor so we never interpolate between a moment and 0.
        n_pmom = np.floor(n_pmom_dec).astype(int)
        n_pmom_max = np.max(n_pmom)
        pmomarray = nc.variables["pmomarray"][:, inu, :n_pmom_max]
        pmomarray = pmomarray.astype("float64")

        # Loop over all the moments to do the same
        # This is the slowest code block
        pmom = np.zeros((n_reff, n_nu, n_pmom_max))
        for j in range(n_pmom_max):
            f_pmom = interp1d(wnum_vec, pmomarray[:, :, j])
            pmom[:, :, j] = f_pmom(nu)

        nc.close()

    return (n_pmom, pmom, reff, w0, qext)
