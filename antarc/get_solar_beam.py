#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Built-in modules
import numpy as np


def get_solar_beam(wnum, solarbeam_file):
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
