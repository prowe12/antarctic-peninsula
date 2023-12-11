#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec 1 2020

@author: prowe
"""

# Built-in modules
import numpy as np


class Blackbody:
    """
    Class for Planck blackbody function
    """

    # function f = plancknu(nu_icm,T);
    #
    # spectral Planck function as function of wavenumbers (cm-1)
    #
    # [h]    = J*s
    # [c]    = m/s
    # [cbar] = cm/s
    # [k]    = J*K-1
    #
    #    Note: LBLRTM uses ...
    #    Constants from NIST 01/11/2002
    #
    #      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
    #     *     CLIGHT / 2.99792458E+10 /,
    #
    #
    # h    = 6.62606896e-34				# J s;  CODATA 2006
    # c    = 2.99792458e8				# m/s;  NIST
    # k    = 1.3806504e-23				# J K-1; CODATA 2006
    # cbar = 2.99792458e10			    # cm/s
    # h_2_cbar4_over_c2 = 1.1910427584934559e-08
    # h_cbar_over_k = 1.4387751601679204
    #
    # h    = 6.62606957e-34				# J s;  CODATA 2010
    # k    = 1.3806488e-23				# J K-1; CODATA 2010
    # h_2_cbar4_over_c2 = 1.1910428681415877e-08
    # h_cbar_over_k = 1.4387769599838156
    #
    # top = h * cbar**3 * 2 * nu**3
    # bottom = c**2 *  ( np.exp(h*cbar*nu/(k*T))-1 )
    # f = cbar * top/bottom
    #
    # h_2_cbar4_over_c2 = h * cbar**3 * 2 *cbar / c**2
    # h_cbar_over_k = h * cbar / k

    def __init__(self, nu_icm: np.ndarray, temp: float):
        """
        Get parts of blackbody function that only depend on wavenumber
        @param nu_icm  Wavenumber (inverse cm): numpy array or float
        @param temp  Temperature (K): float
        """

        # QC: nu_icm must be a numpy array
        if not isinstance(nu_icm, np.ndarray):
            raise TypeError("Blackbody: nu_icm must be numpy array")

        # QC: temperature must be float
        if not isinstance(temp, float):
            raise TypeError("Blackbody: temp must be float")

        h_2_cbar4_over_c2 = 1.1910428681415877e-08
        h_cbar_over_k = 1.4387769599838156

        self.temp = temp
        self.top = h_2_cbar4_over_c2 * nu_icm**3
        self.exp_arg = h_cbar_over_k * nu_icm

        # Catch where wavenumber is 0 and set to small number
        if np.any(nu_icm == 0):
            ized = np.where(nu_icm == 0)[0]
            self.top[ized] = 0.0
            self.exp_arg[ized] = 0.1  # Does not matter b/c numerator is 0

        self.btemp = self.top / (np.exp(self.exp_arg / temp) - 1)

    def set_planckrad(self, temp: float):
        """
        Get the Planck radiance
        @param temp  The temperature
        """
        self.temp = temp
        self.btemp = self.top / (np.exp(self.exp_arg / temp) - 1)

        # Units:
        # [rad] = cm*s-1*J*s*cm3*s-3*cm-3*m-2*s2 = W m-2 cm

    def get_planckrad(self, temp: float) -> np.ndarray:
        """
        Return the Planck radiance
        @param temp  The temperature
        @return  The Planck function radiance (W m-2 cm)
        """
        self.set_planckrad(temp)
        return self.btemp


def plancknu(nu_icm, temp: float):
    """
    Get the Planck function of wavenumber for a given temperature
    @param nu_icm  Wavenumber
    @param temp  Temperature
    @return  The Planck function temperature as a function of wavenumber
    """

    # function f = plancknu(nu_icm,T);
    #
    # spectral Planck function as function of wavenumbers (cm-1)
    #
    # [h]    = J*s
    # [c]    = m/s
    # [cbar] = cm/s
    # [k]    = J*K-1
    #
    #
    #    Note: LBLRTM uses ...
    #    Constants from NIST 01/11/2002
    #
    #      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
    #     *     CLIGHT / 2.99792458E+10 /,
    #
    # h    = 6.62606896e-34				# J s;  CODATA 2006
    # c    = 2.99792458e8				# m/s;  NIST
    # k    = 1.3806504e-23				# J K-1; CODATA 2006
    # cbar = 2.99792458e10			    # cm/s
    #
    # h_2_cbar4_over_c2 = 1.1910427584934559e-08
    # h_cbar_over_k = 1.4387751601679204
    #
    # Units:
    # [B] = cm*s-1*J*s*cm3*s-3*cm-3*m-2*s2
    # [B]= W m-2 cm

    # QC: nu_icm must be float or numpy array
    if not isinstance(nu_icm, (float, np.float32, np.float64, np.ndarray)):
        raise TypeError("plancknu: nu_icm must be np array or float")

    # QC: temperature must be float
    if not isinstance(temp, (float, np.float32)):
        raise TypeError("plancknu: temp must be float")

    h_2_cbar4_over_c2 = 1.1910428681415877e-08
    h_cbar_over_k = 1.4387769599838156
    bottom = np.exp(h_cbar_over_k * nu_icm / temp) - 1

    # Catch where wavenumber is 0 and set to small number
    if np.any(nu_icm == 0):
        ized = np.where(nu_icm == 0)[0]
        bottom[ized] = 1  # Does not matter b/c numerator is 0

    return h_2_cbar4_over_c2 * nu_icm**3 / bottom


def plancknu_temps(nu_icm: np.array, temp: np.array):
    """
    Get the Planck function of wavenubmer and temperature
    @param nu_icm  Wavenumbers
    @param temp  Temperatures numpy array or list only
    @return  The Planck function as a function of wavenumber or temp
    """

    # function f = plancknu(nu_icm,temp);
    #
    # spectral Planck function as function of wavenumbers (cm-1)
    #
    # [h]    = J*s
    # [c]    = m/s
    # [cbar] = cm/s
    # [k]    = J*K-1
    #
    #
    #    Note: LBLRTM uses ...
    # c    Constants from NIST 01/11/2002
    #
    #      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
    #     *     CLIGHT / 2.99792458E+10 /,
    #
    # Old contstants
    # h    = 6.62606896e-34				# J s;  CODATA 2006
    # c    = 2.99792458e8				# m/s;  NIST
    # k    = 1.3806504e-23				# J K-1; CODATA 2006
    # cbar = 2.99792458e10			    # cm/s

    # QC: nu_icm must by float or numpy array
    if not isinstance(nu_icm, np.ndarray):
        raise TypeError("plancknu_temps: nu_icm must be np.ndarray")
    # QC: temp must by numpy array or list
    if not isinstance(temp, np.ndarray) and not isinstance(temp, list):
        raise TypeError("plancknu_temps: temp must be np.ndarray or list")

    # Catch where wavenumber is 0 and set to small number
    fix_wnum_of_0 = False
    if np.any(nu_icm == 0):
        fix_wnum_of_0 = True
        ized = np.where(nu_icm == 0)[0]

    # Constants
    h_2_cbar4_over_c2 = 1.1910428681415877e-08
    h_cbar_over_k = 1.4387769599838156

    top = h_2_cbar4_over_c2 * nu_icm**3
    exp_arg = h_cbar_over_k * nu_icm

    ntemps = len(temp)
    btemp = np.zeros((len(nu_icm), ntemps))
    for i in range(ntemps):
        bottom = np.exp(exp_arg / temp[i]) - 1
        if fix_wnum_of_0:
            bottom[ized] = 1e-4
        btemp[:, i] = top / bottom

    # [B] = cm*s-1*J*s*cm3*s-3*cm-3*m-2*s2
    # [B] = W m-2 cm

    return btemp


def plancknu_slow(nu_icm, temp: float):
    """
    Get the Planck function of wavenumber *or* temperature, starting from
    fundamental constants (slower method)
    @param nu_icm  Wavenumber numpy array or float
    @param temp  Temperature  numpy array or float
    @return  The Planck function as a function of wavenumber or temp
    Notes:
       1) Only nu_icm or temp can be a vector/np array, not both!
          The other must be a float
       2) Constants updated from
          h = 6.62606896e-34				# J s;  CODATA 2006
          k = 1.3806504e-23				# J K-1; CODATA 2006
    """

    # QC: nu_icm cannot be a list
    if not isinstance(nu_icm, np.ndarray) and not isinstance(nu_icm, float):
        raise TypeError("plancknu_slow: nu_icm must be np array or float")

    # QC: temperature must be float
    if not isinstance(temp, float):
        raise TypeError("plancknu_slow: temp must be float")

    # Constants
    h = 6.62606957e-34  # J s;  CODATA 2010
    k = 1.3806488e-23  # J K-1; CODATA 2010
    c = 2.99792458e8  # m/s;  NIST
    cbar = 100 * c  # cm/s

    top = h * cbar**3 * 2 * nu_icm**3
    bottom = c**2 * (np.exp(h * cbar * nu_icm / (k * temp)) - 1)

    # Catch where wavenumber is 0 and set to small number
    if np.any(nu_icm == 0):
        ized = np.where(nu_icm == 0)[0]
        bottom[ized] = 1e-4

    btemp = cbar * top / bottom

    return btemp
