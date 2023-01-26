"""
Created on 2019/01/10

@author: prowe

#
#  Copyright 2017 by Penny M. Rowe and NorthWest Research Associates.
#  All rights reserved.

"""


import numpy as np

def hypsometric(Z, T, RH, Po):

    # Inputs
    # Z: height in m
    # T: temperature in Kelvin
    # RH: relative humidity (%)
    # Po: surface pressure in mb
    # inisZ_P (character string, either 'z' or 'P'), to tell whether
    #   in is height or pressure
    #
    #
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
    
    
    # constants
    Rd = 287.058
    g = 9.80665
    eps = 0.622
        
    # input is z, output needs to be P
    P = np.zeros(len(Z)) + Po 
    
    e = RH/100 * esw(T)          # I assume wrt water, but should be?
    
    for i in range(len(T)-1):
        w = eps * (e/(P-e))
        Tv = T * ((1+(w/eps))/(1+w))
        TvBar = (Tv[i] + Tv[i+1]) / 2
        dz = Z[i+1] - Z[i]
        P[i+1:] = P[i] * np.exp(-dz / (Rd*TvBar/g))  #disp(P(i+1));
        
        '''
        dP = 1e6
        while np.abs(dP) > .1:
            Pprev = P[i]
            w = eps * (e/(P-e))
            Tv = T * ((1+(w/eps))/(1+w))
            TvBar = (Tv[i] + Tv[i+1]) / 2
            dz = Z[i+1] - Z[i]
            P[i+1:] = P[i] * np.exp(-dz / (Rd*TvBar/g))  #disp(P(i+1));
            dP = P[i+1] - Pprev
        '''
            
            
    return P
    

def hypsometric_for_z(P, T, RH, zo):

    # Inputs
    # Z: height in m
    # T: temperature in Kelvin
    # RH: relative humidity (#)
    # Po: surface pressure in mb
    # inisZ_P (character string, either 'z' or 'P'), to tell whether
    #   in is height or pressure
    #
    #
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
    
    
    # constants
    Rd = 287.058
    g = 9.80665
    eps = 0.622
        
    # input is P, output needs to be z
    z = np.zeros(len(P)) + zo
    
    e = RH/100 * esw(T)          # I assume wrt water, but should be?
    
    for i in range(len(T)-1):
        w = eps * (e/(P-e))
        Tv = T * ((1+(w/eps))/(1+w))
        TvBar = (Tv[i] + Tv[i+1]) / 2
        
        # dz = Z[i+1] - Z[i]
        # P[i+1:] = P[i] * np.exp(-dz / (Rd*TvBar/g))  #disp(P(i+1));
        # dP = np.exp(-dz / (Rd*TvBar/g))
        # ln(dP) = -dz / (Rd*TvBar/g)
        # dz = -ln(dP)*Rd*TvBar/g
        dz = -np.log(P[i+1] / P[i]) * Rd * TvBar / g
        z[i+1] = z[i] + dz
        
        '''
        dP = 1e6
        while np.abs(dP) > .1:
            Pprev = P[i]
            w = eps * (e/(P-e))
            Tv = T * ((1+(w/eps))/(1+w))
            TvBar = (Tv[i] + Tv[i+1]) / 2
            dz = Z[i+1] - Z[i]
            P[i+1:] = P[i] * np.exp(-dz / (Rd*TvBar/g))  #disp(P(i+1));
            dP = P[i+1] - Pprev
        '''
            
            
    return z
    
    
    
def esw(T):
    '''
    input temperature in kelvin.  output of esw in millibars
    '''
    
    #using intergrated clausius-clapyeron equation for h2o mixing ratio
    a = 17.27
    b = 35.86
    es = (3.8 / 0.62197) * np.exp(a * ((T-273)/(T-b)))
    
    return es

