#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 13:44:01 2023

@author: prowe
"""

"""
Created on 2019/01/10

@author: prowe

#
#  Copyright 2017 by Penny M. Rowe and NorthWest Research Associates.
#  All rights reserved.

"""


import numpy as np

from antarc.humidRS80 import humidRS80


def hypsometric(alt, temp, rhw, press0):
    """
    Hypsometric equation
    @param alt altitudes
    @param temp temperatures
    @param rhw relative humidities wrt water
    @param press0 surface pressure

    # code written by Steve H for hypsometric determination
    # of altitude (in meters).
    #
    # Modified by PMR for Python, etc
    #
    # Warning: this is not very accurate, just approximate.

    # The hydrostatic equation:
    #
    # P1/P2 = exp[ g/RT * (z2-z1) ]
    # P2    = P1 * exp[ g/RT * (z2-z1) ]
    """

    # constants
    gas_const = 287.058
    g_const = 9.80665
    eps = 0.622

    # input is alt, output needs to be press
    press = np.zeros(len(alt)) + press0

    # Get the partial pressure of water vapor (e)
    press_h2o = rhw / 100 * esw(temp)  # I assume wrt water, but should be?

    for i in range(len(temp) - 1):
        wrat = eps * (press_h2o / (press - press_h2o))
        temp_v = temp * ((1 + (wrat / eps)) / (1 + wrat))
        temp_vbar = (temp_v[i] + temp_v[i + 1]) / 2
        dalt = alt[i + 1] - alt[i]
        press[i + 1 :] = press[i] * np.exp(
            -dalt / (gas_const * temp_vbar / g_const)
        )  # disp(press(i+1));

        # dP = 1e6
        # while np.abs(dP) > .1:
        #     Pprev = press[i]
        #     w = eps * (press_h2o/(press-press_h2o))
        #     temp_v = temp * ((1+(w/eps))/(1+w))
        #     temp_vbar = (temp_v[i] + temp_v[i+1]) / 2
        #     dz = alt[i+1] - alt[i]
        #     press[i+1:] = press[i] * np.exp(-dz / (gas_const*temp_vbar/g_const))
        #     dP = press[i+1] - Pprev

    return press


def hypsometric_for_z(press, temp, rhw, alt0):
    """
    Use the hypsometric equation to compute altitude
    @param press: Pressure in mb
    @param temp: temperature in Kelvin
    @param rhw: relative humidity (%)
    @param alt0: starting height in m

    Note: the first element of press, temp, rhw must correspond to alt0

    The hydrostatic equation:
    P1/P2 = exp[ g/RT * (z2-z1) ]
    P2    = P1 * exp[ g/RT * (z2-z1) ]
    """

    # constants
    gas_const = 287.058
    g_const = 9.80665
    eps = 0.622

    # input is press, output needs to be alt
    alt = np.zeros(len(press))
    alt[0] = alt0

    # Get the partial pressure of water vapor (e or press_h2o)
    press_h2o = humidRS80(rhw, temp, "rhw", "Pa", press) / 100

    for i in range(len(temp) - 1):
        temp_v = temp / (1 - press_h2o / press * (1 - eps))
        temp_vbar = (temp_v[i] + temp_v[i + 1]) / 2

        # dz = alt[i+1] - alt[i]
        # press[i+1:] = press[i] * np.exp(-dz / (gas_const*temp_vbar/g_const))  #disp(press(i+1));
        # dP = np.exp(-dz / (gas_const*temp_vbar/g_const))
        # ln(dP) = -dz / (gas_const*temp_vbar/g_const)
        # dz = -ln(dP)*gas_const*temp_vbar/g_const
        dalt = (
            -gas_const * temp_vbar / g_const * np.log(press[i + 1] / press[i])
        )
        alt[i + 1] = alt[i] + dalt

    return alt


def esw(temp):
    """
    input temperature in kelvin.  output of esw in millibars
    """

    # using intergrated clausius-clapyeron equation for h2o mixing ratio
    cona = 17.27
    conb = 35.86
    return (3.8 / 0.62197) * np.exp(cona * ((temp - 273) / (temp - conb)))
