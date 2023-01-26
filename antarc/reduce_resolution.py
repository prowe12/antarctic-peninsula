# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 19:55:13 2015

@author: prowe
"""

# Dependencies
import numpy as np

# My modules
from antarc.ft_stuff import ft, ift


def reduce_resolution(wnum_mono, calc_mono, dnu, N, n):
    """
    def reduce_resolution(wnum_mono,calc_mono,dnu,N,n):
    #
    #
    # function [wnc,calc] = reduceRes3(wnum_mono,calc_mono,dnu,N,n)
    #
    # Reduce resolution according to dnu, assuming N pts in interfogram
    #   and using n points to zero-pad the spectrum
    #
    #
    # by Penny M. Rowe, based on aeri.m by Von P. Walden and others
    # 2013/06/06
    #
    # inputs:
    #   wnum_mono
    #   calc_mono
    #   dnu: desired resolution
    #   N: number of points in interferogram
    #   n: large number of points, for interpolation
    #
    #
    # updates:
    #   Modified to take as input dnu, and N instead of laserWn, and n
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    """

    # so the final wavenumber vector needs to be:
    # N=16384;                # number of points in a single-sided interferogram
    # xmax = N/laserWn ;
    # dnu = 1/2/xmax;
    # numax = dnu*(N-1);
    # wnf=(0:dnu:numax);

    # we want to put the calculated interferogram on the same spacing
    # as the measured interferogram so we can chop it at the same
    # place

    # The calculated spectrum needs to begin at zero wavenumbers
    #   (it typically has a spacing of 0.0003 cm-1)
    #   and it needs to end at numax to give dx properly
    bwn = wnum_mono[0]
    ewn = wnum_mono[-1]

    if n == 0 or np.isnan(n):
        n = 2**22  # large number of points

    numax = dnu * N  # dnu = lnu/N/2
    # so lnu = dnu*N*2
    # and numax=lnu/2 or dnu*N

    wnbr = np.linspace(0, numax, n + 1)  # e.g. 0.030s
    wnbr = wnbr.T

    # .. Create new radiance vector with tapered ends.
    #    Create np.zeros array (0.019 s)
    calc_tmp = np.zeros(n + 1)
    ind = np.where(np.logical_and(wnbr >= bwn, wnbr <= ewn))[0]
    # i1 = int(ind[1])

    # Interpolate onto the new grid (0.014 s)
    # this should probably use the cubic interpolation, but there
    # is not a good one in Python and it would take even longer.
    # What are the resulting errors?
    # linear interpolation: 0.013 s for first layer
    # cubic spline: 0.07 to 0.12 s for first - 37th layer
    # Cubic spline adds about 5s
    # However, the improvement is not noticeable.
    # start = time.time()
    calc_tmp[ind] = np.interp(wnbr[ind], wnum_mono, calc_mono)
    # spl = interpolate.splrep(wnum_mono, calc_mono)
    # calc_tmp[ind] = interpolate.splev(wnbr[ind], spl)
    # print('time to interpolate onto new grid: ' + str(time.time() - start))

    # start = time.time()
    # Low-wavenumber roll-off
    # set up index vectors (0.012 s)
    ind_lwn = np.arange(ind[1])
    # print(ind_lwn.shape)
    # ind_lwn = np.arange(1e5); print(ind[1])
    # ind_hwn = np.arange(ind[-1], n+1); print(ind[-1])
    ind_hwn = np.arange(ind[-1], ind[-1] + 3000)
    # print(ind_hwn.shape) #n+1); print(ind[-1])

    # Compute and paste on the low-wavenumber roll-offs (0.13 s)
    hwn = ind_hwn - ind[-1] + 1

    # ind_lwn
    # i1-100000:i1
    calc_tmp[ind_lwn] = (
        ((np.cos(ind_lwn * np.pi / ind[0]))[::-1] + 1) / 2 * calc_tmp[ind[0]]
    )
    calc_tmp[ind_hwn] = (
        (np.cos(hwn * np.pi / (len(hwn))) + 1) / 2 * calc_tmp[ind_hwn[0]]
    )
    # print('time for low-wavenumber roll-off: ' + str(time.time() - start))

    # IFT: slow step (e.g. 1.1 s)
    # start = time.time()
    x, ifg = ift(wnbr, calc_tmp, 1, 1)
    # print('time for ift: ' + str(time.time() - start))

    # FT chop down to N+1 (0.002 s)
    wnc, calc = ft(x[: N + 1], ifg[: N + 1], 1, 1)

    #
    # .. Extract the portion of the spectrum between the original
    #    wavenumber limits. (0.0002 s)
    ind = np.where(np.logical_and(wnc >= bwn, wnc <= ewn))[0]
    wnc = wnc[ind]
    calc = calc[ind]

    return (wnc, calc)
