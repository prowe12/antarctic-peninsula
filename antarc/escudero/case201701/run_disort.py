"""
run_disort.py

Copyright 2018-2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

Purpose:
  A code for calling disort_driver_py for running DISORT.
"""

# Built in modules
import numpy as np
from scipy.interpolate import interp1d

# My modules
from antarc.disort_driver_py import disort_driver



def run_disort(nu, dtau_gas, cldlyr,
    reff_liq, total_od_vis_liq,
    reff_ice, total_od_vis_ice,
    itemp_liq, od_vis_liq, itemp_ice, od_vis_ice,
    iliq_layer, iice_layer,
    ssp_liq, ssp_ice, temper, nstr_in, umu0, umu,
    albedo_nu, fbeam_nu, delv_nu, iobs):
    """
    Purpose: Run DISORT for multi-layer liquid and ice clouds.

    Inputs:
     nu:               wavenumber vector
     dtau_gas:         layer optical depths of gas
     total_od_vis_liq: Total optical depth in visible for liquid
     total_od_vis_ice: Total optical depth in visible for liquid
     itemp_liq:        Inds to ssp files that are used, for liquid
     itemp_ice:        Inds to ssp files that are used, for liquid
     od_vis_liq:       od in visible for liquid, per layer, per ssp
     od_vis_ice:       od in visible for ice, per layer, per ssp set
     reff_liq:         effective radius of liquid (scalar)
     reff_ice:         effective radius of ice (scalar)
     iliq_layer:
     iice_layer:
     cldlyr_liq:       vector for number of layers
     cldlyr_ice:       vector, for number of layers
     iobs:             observer height index

    Input variable types and sizes:
     nu:               (vector)
     dtau_gas:         (nus x layers)
     total_od_vis_liq: (scalar)
     total_od_vis_ice: (scalar)
     itemp_liq:        (x layers)
     itemp_ice:        (x layers)
     od_vis_liq:       (layers x ssps)
     od_vis_ice:       (layers x ssps)

    Disort parameters (asterisks indicate variables that change with nu):
    From DISORT.txt:
      ANGLE CONVENTION:
        Polar (zenith) angles are measured from the upward direction:
        straight up is zero degrees and straight down is 180 degrees.
        There is a small inconsistency in that, for historical reasons,
        the cosine of the incident beam angle (umu0) is positive,
        whereas according to this convention it should be negative.
      nlyr                     Number layers, INGEGER*
      dtauc                    Layer optical depths, REAL [1:MAXCLY]*
      ssalb                    Single scattering albedo, REAL [1:MAXCLY]*
      nmom                     Number of moments, INTEGER*
      ncldlyr                  Number of cloud layers, INTEGER
      cldlyr+1                 cloud layers TOA to surf, INTEGER [1:MAXCLY]
      Pmom_cld                 Legendre Poly, REAL [0:maxmom,1:MAXCLY]*
      temper = (input)         bndry Ts TOA-to-surface REAL [0:MAXCLY]
      WVNUMLO                  REAL*
      WVNUMHI,                 REAL*
      usrtau = 1               1=>specify taus, LOGICAL
      ntau   = 1               Number output levels?  INTEGER
      UTAU   =                 Output levels, as sum(dtauc), REAL [MAXULV]*
      nstr   = (input)         Number of streams*
      USRANG = 1               True => specify observation zenith angle
      NUMU   = 1               Number of observation zenith angles
      umu    = (input)         Cosine of output angles in increasing order:
                               neg (downward) to pos (upward) values, [1:MAXUMU]
      NPHI   = 1               Number of observation azimuthal angles
      PHI    = [0.]            Observation azimuthal angle
      IBCND  = 0               0 => General boundary condition flag
      FBEAM  =                 solar, depends on fbeam_nu*
      umu0   = (input)         Cosine of incident beam angle, REAL
      PHI0   = 0.              Azimuth angle of incident beam
      FISOT  = 0.              Isotropic incident flux
      LAMBER = 1               True => The surface is Lambertian
      ALBEDO =                 ALBEDO_long[nu], from input, REAL*
      BTEMP                    Bottom T (e.g. surface T), REAL
      TTEMP =                  temper[i1], Top T (NA if TEMIS = 0), REAL*
      TEMIS  = 0.              Top emissivity, REAL
      plank  = 1               True => thermal emission included, INTEGER
      ONLYFL = 0               False => Return azimuthally averaged
                               intensities (not just fluxes), INTEGER
      HEADER = '...'           unused string, CHAR STRING
      MAXCLY =                 len(dtauc), INTEGER
      MAXULV =                 len(UTAU), number output levels, INTEGER
      MAXUMU =                 len(umu), number umu?, INTEGER
      MAXPHI =                 len(PHI), number PHI, INTEGER

      So we should have:
      dtauc_nu
      ssalb_nu
      Pmom_cld_nu
      utau_nu
      albedo_nu (input)
      fbeam_nu (input)
    """

    # # # # # # # # #     Debugging     # # # # # # # # # # #
    # If the following flag is set to True, the disort inputs
    #    will be printed to a file for each wavenumber bin
    #    to allow debugging.
    #
    #    Note: do not leave this set to True if not debugging, as it
    #    will slow the code down!
    debug_flag = False          # If True, will print DISORT inputs to log
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



    # # # # # # # # #     Set parameters     # # # # # # # # # # #

    # For each nu and Temp, there can be a different number of Pmoms
    #    So we need to set the maximum allowed and also set the maximum
    #    ever used (starting with 0, here)
    disort_warning = False      # Will become true if DISORT returns a message
    maxpmom = 500               # Maximum number of Pmoms allowed

    # Number of layers and frequencies
    nlyr_long, nnus = dtau_gas.shape

    # DISORT inputs that are hard-wired
    plank  = 1              # True => thermal emission included
    usrtau = 1              # True => we want to specify taus
    ntau   = 1              # There will be one tau
    USRANG = 1              # True => specify observation zenith angle
    NUMU   = 1              # Number of observation zenith angles
    NPHI   = 1              # Number of observation azimuthal angles
    PHI    = [0.]           # Observation azimuthal angle
    IBCND  = 0              # 0 => General boundary condition flag
    PHI0   = 0.             # Azimuth angle of incident beam
    LAMBER = 1              # True => The surface is Lambertian
    ONLYFL = 0              # False => Return azimuthally averaged intensities
                            # (not just fluxes)
    FISOT  = 0.             # Isotropic incident flux
    TEMIS  = 0.             # Top emissivity
    HEADER = 'Run DISORT from Python'  # char string, not used by disort.disort
    MAXPHI = 1
    MAXULV = 1
    MAXUMU = len(umu)
    BTEMP  = temper[nlyr_long]   # Bottom temperature


    # Number of temperatures to interpolate SSPs over
    iliq = np.unique(itemp_liq)
    iice = np.unique(itemp_ice)
    nssp_sets_liq = len(iliq)
    nssp_sets_ice = len(iice)

    # Remove zero-d out layers
    if total_od_vis_liq <= 0:
        if total_od_vis_ice <= 0:
            nssp_sets_liq = 0
            nssp_sets_ice = 0
            maxpmom = 5
            nstr_in = 4
            # We have to have at least one cloud layer, put it at the surface
            cldlyr = np.array([nlyr_long-1])
            ncldlyr = 1
        else:
            nssp_sets_liq = 0
            cldlyr = np.array(cldlyr)[iice_layer]
            od_vis_ice = od_vis_ice[iice_layer,:]
    elif total_od_vis_ice <= 0:
        nssp_sets_ice = 0
        cldlyr = np.array(cldlyr)[iliq_layer]
        od_vis_liq = od_vis_liq[iliq_layer,:]


    # Preallocate
    #    The array cldyTau_* holds the cloud od (in visible region) for
    #    liquid and ice for each possible ssp index and each layer
    #    so it is large and mostly full of zeros
    radiance = np.zeros(nnus)
    flux_down = np.zeros(nnus)

    # Preallocate vectors and matrices that depend on number cloud layers
    ncldlyr = len(cldlyr)
    ssalb = np.zeros(nlyr_long)
    w0dtau = np.zeros((ncldlyr, nnus))
    pmom_w0dtau = np.zeros((ncldlyr, nnus, maxpmom))
    pmom = np.zeros((maxpmom, nnus, ncldlyr))
    dtau = 1 * dtau_gas

    # Require at least nstr+1 pmoms    
    # npmom_nu = 5 * np.ones(nnus).astype('int32')
    npmom_nu = (nstr_in + 3) * np.ones(nnus).astype('int32')

    # Loop over the sspfiles we are using (depends on temps and layers)
    # Get liquid cloud component based on weights to ssp files
    # There may be multiple liquid cloud layers, and for each
    # there can be one or two ssp files; so interpolate
    # Interpolate liquid optical properties for frequencies & reff
    # NPmom_max = nstr  # Number of Pmoms must be at least nstr
    for i in range(nssp_sets_liq):
        j = iliq[i] # index to ssp set
        current_ssp = ssp_liq[j]
        qext = current_ssp.fun_qext(reff_liq)
        w0_liq = current_ssp.fun_w0(reff_liq)
        npmom_nui = (np.floor(current_ssp.fun_npmom(reff_liq))).astype('int32')
        npmom_nu = (np.max([npmom_nu, npmom_nui], axis=0)).astype('int32')
        # NPmom_max_i = np.max(npmom_nui)                # scalar
        # NPmom_max = np.max([NPmom_max, NPmom_max_i])    # scalar

        # Loop over wavenumbers
        for inu in range(nnus):
            interp_fun = interp1d(current_ssp.reff, current_ssp.pmom[:,inu,:npmom_nui[inu]].T)
            pmom0 = interp_fun(reff_liq)
            # pmom0 = ssp_liq[j].fPmom[inu](reff_liq)[:npmom_nui[inu]]
            dtau_liq_i_inu = qext[inu] * od_vis_liq[:,j] / 2
            w0dtau0_liq_i_inu = w0_liq[inu] * dtau_liq_i_inu
            dtau[cldlyr,inu] += dtau_liq_i_inu
            w0dtau[:,inu] += w0dtau0_liq_i_inu

            for icld in range(len(cldlyr)):
                pmom_w0dtau[icld, inu, :npmom_nui[inu]] += \
                  w0dtau0_liq_i_inu[icld] * pmom0
                #if np.any(np.isnan(pmom_w0dtau)):
                #    raise ValueError('Nan in pmom')

                # If no ice, then
                #    on last ssp, to get pmom, divide by (dtau*ssalb)
                if (nssp_sets_ice==0) and (i == nssp_sets_liq-1):
                    pmom[:,inu,icld] = pmom_w0dtau[icld,inu,:]/w0dtau[icld,inu]
                    #if np.any(np.isnan(pmom)):
                    #    raise ValueError('Nan in pmom')

    for i in range(nssp_sets_ice):
        j = iice[i]                                     # index to ssp set
        current_ssp = ssp_ice[j]
        qext = current_ssp.fun_qext(reff_ice)
        w0_ice = current_ssp.fun_w0(reff_ice)
        npmom_nui = (np.floor(current_ssp.fun_npmom(reff_ice))).astype('int32')
        npmom_nu = (np.max([npmom_nu, npmom_nui], axis=0)).astype('int32')
        # NPmom_max_i = np.max(npmom_nui)
        # NPmom_max = np.max([NPmom_max, NPmom_max_i])  # scalar

        # Loop over wavenumbers
        for inu in range(nnus):
            interp_fun = interp1d(current_ssp.reff, current_ssp.pmom[:,inu,:npmom_nui[inu]].T)
            pmom0 = interp_fun(reff_ice)
            # pmom0 = ssp_ice[j].fPmom[inu](reff_ice)[:npmom_nui[inu]]
            dtau_ice_i_inu = qext[inu] * od_vis_ice[:,j]/2
            w0dtau0_ice_i_inu = w0_ice[inu] * dtau_ice_i_inu
            dtau[cldlyr,inu] += dtau_ice_i_inu
            w0dtau[:,inu] += w0dtau0_ice_i_inu


            for icld in range(len(cldlyr)):
                pmom_w0dtau[icld, inu, :npmom_nui[inu]] += \
                  w0dtau0_ice_i_inu[icld] * pmom0 #[:npmom_nui[inu]]
                #if np.any(np.isnan(pmom_w0dtau)):
                #    raise ValueError('Nan in pmom')


                # On last ssp, to get pmom, divide by (dtau*ssalb)
                if i == nssp_sets_ice-1:
                    pmom[:,inu,icld] = pmom_w0dtau[icld,inu,:]/w0dtau[icld,inu]
                    #if np.any(np.isnan(pmom)):
                    #    raise ValueError('Nan in pmom')


    # To get ssalb, divide by dtau
    w0 = w0dtau / dtau[cldlyr,:]

    # Solar contribution only if umu0 > 0, otherwise
    #    sun is below the horizon, so turn it off
    if umu0 <= 0:
        FBEAM_nu = 0*nu
        umu0 = 0
    elif umu0 > 0:
        FBEAM_nu = fbeam_nu
    else:
        raise NameError('Bad value for umu0')


    # Wavenumbers
    WVNUMLOs = nu - delv_nu/2
    WVNUMHIs = nu + delv_nu/2


    # Make sure there are at least as many Pmoms as strings
    nstr = 1 * nstr_in
    #if np.any(npmom_nu < nstr+1):
    #    print('got one')
    npmom_nu[npmom_nu<nstr+1] = nstr+1

    #  Loop over the frequencies
    izm = np.zeros((nnus,1))
    for inu in range(nnus):
        uu = np.nan # This logic because disort doesn't always converge

        # Values for this wavenumber
        WVNUMLO = WVNUMLOs[inu]
        WVNUMHI = WVNUMHIs[inu]
        FBEAM = FBEAM_nu[inu]
        dtauc = dtau[:,inu]                             # x nlyr
        ssalb[cldlyr] = w0[:,inu]                       # x ncldlyr
        PMOM_CLD = pmom[:npmom_nu[inu],inu,:]       # nmom x cldlyr

        #if np.any(np.isnan(PMOM_CLD)):
        #    raise ValueError('Nan in pmom')

        i2 = len(dtauc)

        # For an observer layer very close to the surface (within a few m),
        #    The optical depth in-between may be < 1e-5 because the layer
        #    is so thin. Furthermore, the layer between the surface and the
        #    observer plays a minior role for IR downwelling radiance,
        #    as it is only impacts radiation scattered off the surface and
        #    then scattered back down again to the observer. So, if the
        #    optical depth of this layer is too small for single precision
        #    accuracy, it's best to just remove it.
        #
        #    Here we check if:
        #    1) the observer layer is specified: iobs > 0,
        #    2) the observer layer is above the: iobs <= len(dtauc),
        #    3) The optical depths for layers below iobs to the surface layer
        #     are of the order of the precision: (np.sum(dtauc[iobs:i2])<1e-5)
        #
        #    If 1-3 all true, remove the surface layers by setting i2 = iobs
        # if (iobs > 0) and (iobs < len(dtauc)) \
        #   and (np.sum(dtauc[iobs:i2])<1e-5):
        if (0 < iobs < len(dtauc)) and (np.sum(dtauc[iobs:i2])<1e-5):
            raise ValueError('Increase optical thickness of surface layer')

        # Some of the upper layers may also have optical depths that
        # are on the order of the precision (od<5e-7 or 1e-5). Remove those
        # by increasing i1 from 0.
        izeds = np.where(dtauc[:i2]<1e-5)[0]
        if len(izeds)==0:
            i1 = 0
        elif (izeds[0] != 0) or (np.any(np.diff(izeds)!=1)):
          # print('Warning: Redesign your layers so od gets smaller going up.');
            i1 = np.max(izeds)+1
        else:
            i1 = np.max(izeds)+1
            
        if i1 == i2:
            msg = f"i1=i2={i1}, so no atmosphere, likely because too many " \
                  "optical depths are <=1e-6. \n Near-surface optical " \
                  f"depths are {dtauc[-1]} {dtauc[-2]}  {dtauc[-3]} " \
                  f"{dtauc[-4]}.\n nu = {nu[inu]}"
            raise NameError(msg)
        izm[inu] = i1
        nlyr = i2-i1
        try:
            CLDLYR_inu = np.array(cldlyr) - i1 + 1
        except:
            raise ValueError('This should not happen')
            # CLDLYR_inu = np.array(cldlyr-i1) + 1

        # Height of observer set with index iobs (index to layer),
        # iobs = len(dtauc) => surface, cumulative od of atmosphere
        # iobs = 0 => TOA, UTAU = 0
        if (iobs > 0) and (iobs <= len(dtauc)):
              UTAU = [np.sum(dtauc[i1:iobs])]
        elif iobs == 0.:
            UTAU =[0.]
        else:
            raise NameError('Option iobs = ' + str(iobs) + ' not allowed.')

        # Old maxes not used, b/c sizes must be exact
        # e.g. MAXCLY = 120; MAXPHI = 1; MAXULV = 2; MAXUMU = 10
        MAXCLY = nlyr
        maxmom = PMOM_CLD.shape[0]-1
        nmom = maxmom
        #nstr = 1 * nstr_in
        #if nstr > nmom:
        #    if nmom % 2 == 0:
        #        nstr = nmom
        #    else:
        #        nstr = nmom - 1
        #    if nstr<=3:
        #        raise NameError('nstr  is <=3')
        # We want at least 16 streams
        if nstr < 16:
            msg = 'Wavenumber: ' + str(WVNUMLO) + ', Re, liq: ' \
                   + str(reff_liq)  + ', Re, ice: ' +str(reff_ice) \
                   + ', streams: ' + str(nstr)
            raise ValueError(msg)

        if debug_flag: # and inu==0:
            dlog = open("disort_output.txt", "w")
            print("WVNUMLO, WVNUMHI", WVNUMLO,WVNUMHI, file = dlog)
            print("umu0 = ", umu0, file = dlog)
            print("MAXULV, MAXUMU, MAXPHI = ",MAXULV,MAXUMU,MAXPHI,file=dlog)
            print("PHIP, FISOT, LAMBER = ", PHI0, FISOT, LAMBER, file=dlog)
            print("TEMIS, plank, ONLYFL = ", TEMIS,plank,ONLYFL, file = dlog)
            print("HEADER = ", HEADER, file = dlog)
            print("ALBEDO = ", albedo_nu[inu], file = dlog)
            print("FBEAM = ", FBEAM, file = dlog)

            print("nlyr = ", nlyr, file = dlog)
            print("MAXCLY = ", MAXCLY, file = dlog)
            print("temper.shape = ", temper[i1:i2+1].shape, file = dlog)
            print("temper = ", temper[i1:i2+1], file = dlog)
            print("TTEMP = ", temper[i1], file = dlog)
            print("BTEMP =", BTEMP, file = dlog)

            print("usrtau,ntau,UTAU = ", usrtau,ntau,UTAU, file = dlog)
            print("USRANG,NUMU,umu = ", USRANG,NUMU,umu, file = dlog)
            print("NPHI,PHI,IBCND = ", NPHI,PHI,IBCND, file = dlog)
            print("DTAU = ", dtauc[i1:i2], file = dlog)
            print("DTAU.shape = ", dtauc[i1:i2].shape, file = dlog)
            print("ssalb = ", ssalb[i1:i2], file = dlog)
            print("ssalb.shape = ", ssalb[i1:i2].shape, file = dlog)
            print("ncldlyr = ", ncldlyr, file = dlog)
            print("pmom.shape[1] = ", PMOM_CLD.shape[1], file = dlog)
            print("cldlyr.shape = ", CLDLYR_inu.shape, file = dlog)
            print("cldlyr = ", CLDLYR_inu, file = dlog)

            print("nstr =", nstr, file = dlog)
            print("nmom = ", nmom, file = dlog)
            print("maxmom = ", maxmom, file = dlog)
            print("pmom.shape[0] = ", PMOM_CLD.shape[0], file = dlog)
            print("pmom = ", PMOM_CLD, file = dlog)
            dlog.close()

            # Quality control
            if np.any(PMOM_CLD==np.nan) or np.any(ssalb[i1:i2]==np.nan) \
              or np.any(dtauc[i1:i2]==np.nan) or np.any(dtauc[i1:i2]<=0):
                  raise ValueError('Bad inputs to DISORT!')


        while np.isnan(uu):
            # Call disort
            # We must add 1 to cldlyr below because python indexes from
            # zero but in the fortran code it is indexed from one
            # PMR, 2016/08/28
            _, rfldn, flup,  _, _, uu,  _, _, errmsg, errflag = disort_driver(
                nlyr, dtauc[i1:i2], ssalb[i1:i2], nmom, ncldlyr, CLDLYR_inu, \
                PMOM_CLD, temper[i1:i2+1],  WVNUMLO, WVNUMHI, usrtau, ntau, \
                UTAU, nstr, USRANG, NUMU, umu,  NPHI, PHI, IBCND, FBEAM,
                umu0, PHI0, FISOT, LAMBER,  albedo_nu[inu], BTEMP, \
                temper[i1], TEMIS, plank,  ONLYFL, HEADER, MAXCLY, \
                MAXULV, MAXUMU, MAXPHI, maxmom)

            # If DISORT returns a fatal error (indicated by errflag), quit
            if errflag:
                raise ValueError(errmsg)


        #    If DISORT returns a warning, print the wavenumber and warning
        #    to the log file and set output variable disort_warning to true
        errmsg = (errmsg.decode("utf-8")).rstrip()
        if errmsg != '':
            disort_warning = True
            print(str(round(WVNUMLO,2)) + ' cm-1: ' + errmsg)


        radiance[inu] = uu
        flux_down[inu] = rfldn

    radiance = radiance*1e3 # radiance units (milliwatts)


    return radiance, izm, flux_down, flup, disort_warning
