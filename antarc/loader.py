#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 10:57:24 2020

@author: prowe

This file contains the loaders for ancillary data used by CLARRA
(except for radiances and clouds, which have loaders in separate files)
Thus if the input file formats differ, they must be modified to conform to
the loaders, or, alternatively, the respective loader must be modified.

When creating files to be loaded, be sure the fields and units are as
specified below!  If not given, or unsure, check the sample files in the
sampleRun directory.
"""

# Built-in modules
import logging
import scipy.io as spio
from netCDF4 import Dataset, num2date, date2num  # pylint: disable=E0611

# from scipy.io.netcdf import NetCDFFile as DS
from scipy.interpolate import interp2d
import numpy as np
import cftime


log = logging.getLogger(__name__)


def load_ssp_nu(datafile, wnum):
    """
    Return a tuple with the single-scattering parameter data
    @param datafile The filename of the file with the ssp data
    @param wnum The wavenumbers of interest
    @return n_pmom  The number of legendre moments
    @return pmom  The legendre moments
    @return reff  The effective radii
    @return w0  The single-scattering albedo
    @return qext  The extinction efficiency Qext

    File requirements:
    It is suggested to use the files created for this project. Contact the
    authors for access to all available single-scattering parameter files.

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

    Warning: File format and data must be correct. No error-checking is done here!
    """
    dnu = wnum[1] - wnum[0]
    # with DS(datafile, 'r') as nc:
    with Dataset(datafile, "r", format="NETCDF4_CLASSIC") as nc:

        inu = np.where(
            np.logical_and(
                nc.variables["wnum_list"][:] > wnum[0] - 2 * dnu,
                nc.variables["wnum_list"][:] < wnum[-1] + 2 * dnu,
            )
        )[0]
        n_pmomarray = nc.variables["Npmomarray"][:, inu].astype("int32")
        w0_mesh = nc.variables["w0_mesh"][:, inu].astype("float64")
        qext_mesh = nc.variables["qext_mesh"][:, inu].astype("float64")
        reff_list = nc.variables["reff_list"][:].astype("float64")
        reff = reff_list[:, 0].copy()
        wnum_vec = nc.variables["wnum_list"][inu].astype("float64")

        # Set up empty output arrays
        n_nu = wnum.size
        n_reff = reff.size
        qext = np.zeros((n_reff, n_nu))
        w0 = np.zeros((n_reff, n_nu))
        n_pmom_dec = np.zeros((n_reff, n_nu))

        # Interpolate qext, w0, get an interpolated number of moments!
        f_q = interp2d(reff, wnum_vec, qext_mesh.T)
        f_w = interp2d(reff, wnum_vec, w0_mesh.T)
        f_npmom = interp2d(reff, wnum_vec, n_pmomarray.T)

        for i in range(n_reff):
            qext[i, :] = f_q(reff[i], wnum)[:, 0]
            w0[i, :] = f_w(reff[i], wnum)[:, 0]
            n_pmom_dec[i, :] = f_npmom(reff[i], wnum)[:, 0]

        # Use floor so we never interpolate between a moment and 0.
        n_pmom = np.floor(n_pmom_dec).astype(int)
        n_pmom_max = np.max(n_pmom)
        pmomarray = nc.variables["pmomarray"][:, inu, :n_pmom_max]
        pmomarray = pmomarray.astype("float64")

        # Loop over all the moments to do the same
        # This is the slowest code block
        pmom = np.zeros((n_reff, n_nu, n_pmom_max))
        for j in range(n_pmom_max):
            f_pmom = interp2d(reff, wnum_vec, pmomarray[:, :, j].T)
            pmom[:, :, j] = [f_pmom(x, wnum)[:, 0] for x in reff]

        nc.close()

    return (n_pmom, pmom, reff, w0, qext)


def load_profiles(prof_file):
    """
    Load in the profile file
    @param prof_file The name of the profile file
    @return height  Altitudes or heights (km)
    @return pressure  Pressures corresponding heights (mb)
    @return temperature  Temperatures corresponding heights (K)
    @return h2o  Water vapor concentration corresponding to heights (ppmv)

    File requirements:
    fname must be the pathname of a netcdf file in format NETCDF4
    The fields must include the following, with the units as given:
        time: with units specified in netcdf format
        temp_mean: surface temperatures at times, in Celsius
        atmos_pressure: surface atmospheric presure, in mb
        temp_10m: (Optional): surface temperature at 10 m, in Celsius

    Warning: There is no error-checking here!  Run quality control on all
    profiles before using this code.
    """

    # .. Load in the profile
    with Dataset(prof_file, "r", format="NETCDF4_CLASSIC") as nc:
        height = np.double(nc["z"][:].data)
        pressure = np.double(nc["P"][:].data)
        temperature = np.double(nc["T"][:].data)
        h2o = np.double(nc["h2o"][:].data)

    return height, pressure, temperature, h2o


def load_od_gas(odfile):
    """
    Load in radiative transfer parameters for gases, at effective resolution
    @param odfile  The filename
    @return date  The date
    @return view_angle  The viewing angle
    @return wnum  The wavenumber (cm-1)
    @return rads  The radiances up to a desired max height (mW/(m2 sr-1 cm-1) or RU)
    @return tsc  The transmittances from surface to the desired max height
    @return rad_above  The radiance above the desired max height (RU)
    @return Bctc  The Planck function times the transmittance (RU * [trans])
    @return dt_dtau  The layer transmittance times the layer optical depth

    File requirements: Files should be created by clear_sky_sims_runner.py
    or with code that gives the equivalent file structure and units

    Warning: File format and data must be correct. No error-checking is done here!
    """
    odinfo = spio.loadmat(odfile)
    date = odinfo["date"]
    wnum = odinfo["nu"][0]
    rads = odinfo["rads"]
    tsc = odinfo["tsc"]
    view_angle = odinfo["view_angle"]
    bctc = odinfo["Bc_tsc"]
    dt_dtau = odinfo["dt_dtau"]
    rad_above = odinfo["rad_above"][0]

    return date, view_angle, wnum, rads, tsc, rad_above, bctc, dt_dtau


def load_surfmet(fname):
    """
    Load surface meteorology data
    @param fname  The filename for the surface met date
    @return surf_dnum  The datenumber as specified by t_units
    @return surf_t  The temperature (K)
    @return surf_p  The pressure (mb)
    @return surf_t_10m  The 10-meter temperature if present (K), else nans

    File requirements:
    fname must be the pathname of a netcdf file in format NETCDF4
    The fields must include the following, with the units as given:
        time: with units specified in netcdf format
        temp_mean: surface temperatures at times, in Celsius
        atmos_pressure: surface atmospheric presure, in mb
        temp_10m: (Optional): surface temperature at 10 m, in Celsius

    Warning:
        File format and data must be correct. No error-checking is done here!
        The times should be cftime.DatetimeGregorian, and then the following
        will give a value equivalent to the timestamp of a timezone-aware
        datetime in UTC. Be sure to unit test this by adding an example
        from your code to test_loadsurfmet.py
    """
    with Dataset(fname, "r", format="NETCDF4") as nc:
        # Note: num2date is slow
        dtimes = num2date(nc.variables["time"][:], nc.variables["time"].units)
        if not isinstance(dtimes[0], cftime.DatetimeGregorian):
            raise TypeError(
                "load_surfmet currently only handles cftime.DatetimeGregorian"
            )
        surf_time = dtimes.data

        # Convert to same number as timestamp would
        t_units = "seconds since 1970-01-01 00:00:00 0:00"
        surf_dnum = [date2num(x, t_units) for x in surf_time]

        surf_t = nc.variables["temp_mean"][:] + 273.15
        surf_p = nc.variables["atmos_pressure"][:]
        if "temp_10m" in nc.variables:
            surf_t_10m = nc.variables["temp_10m"][:] + 273.15
        else:
            surf_t_10m = np.nan + np.zeros(np.shape(surf_t))

    return surf_dnum, surf_t, surf_p, surf_t_10m


def get_surface_albedo_from_file(surface_emissivity_file):
    """
    Return the surface albedo data (wavenumber, albedo)
    @param surfEmissDataFile  The filename for the surface emissivity data
    @return surf_albedo_data  Wavenumber and surface albedo

    File requirements:
    surface_emissivity_file must be the pathname of a text file
    The file may have comments at the top, with % at the beginning of each line
    columns:
        - The first column (index 0) is ignored.
        - The second column (index 1) is the wavenumber in cm-1
        - The third column (index 2) is the emissivity
    Warning: File format and data must be correct. No error-checking is done here!
    """
    albedo_data = np.loadtxt(surface_emissivity_file, comments="%")
    surf_albedo_data = np.array([albedo_data[:, 1], 1 - albedo_data[:, 2]]).T
    return surf_albedo_data


def getsolarbeam(wnum, solarbeam_file):
    """
    Return the solar beam in the indicated file for the given wavenumbers
    @param wnum  The wavenumbers
    @param solarbeam_file  The solar beam file
    @return The solar beam at the indicated wavenumbers

    File requirements:
    solarbeam_file must be the pathname of a text file
    The file must have the following first two columns
        - The first column (index 0) is wavenumber.
        - The second column (index 1) is the beam data in the required units
    Warning: File format and data must be correct. No error-checking is done here!
    """
    beamdata = np.loadtxt(solarbeam_file)
    beam = np.interp(wnum, beamdata[:, 0], beamdata[:, 1]) / 1000
    return beam


# def get_files_from_dir(dirname, samplefname):
#     """
#     Get all the files of type sampleFnameStr from directory  dirstr
#     """
#     fnames = []
#     dirstuff = sorted(os.listdir(dirname))
#     for fname in dirstuff:
#         if len(fname) == len(samplefname) and fname[:3] == samplefname[:3] \
#           and fname[-3:] == samplefname[-3:]:
#             fnames.append(fname)
#     #dum = [x for x in fnames if (len(x) == len(samplefname) \
#     #                              and x[:3] == samplefname[:3] \
#     #                              and x[-3:] == samplefname[-3:])]
#     return fnames


# def get_surface_albedo_IR(wnum = None, surface_emiss_file = None):
#     """
#     Return the surface albedo
#     @param wnum  The wavenumbers of interest
#     @param surface_emiss_file  The filename for the surface emissivity data
#     @return surface_albedo  Wavenumber and surface albedo
#     """
#     emiss_data = np.loadtxt(surface_emiss_file)
#     emissivity = np.interp(wnum, emiss_data[:,1], emiss_data[:,2])
#     emissivity[emissivity > 1.] = 1.
#     emissivity[emissivity < 0.] = 0.
#     surface_albedo = 1. - emissivity
#     return surface_albedo

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
