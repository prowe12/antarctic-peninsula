#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Built-in modules
import numpy as np


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
    albedo_data = np.loadtxt(surface_emissivity_file) #, comments="%")
    surf_albedo_data = np.array([albedo_data[:, 1], 1 - albedo_data[:, 2]]).T
    return surf_albedo_data
