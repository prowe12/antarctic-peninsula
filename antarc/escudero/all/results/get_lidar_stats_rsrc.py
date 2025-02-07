#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 08:04:34 2024

@author: prowe
"""

import numpy as np
import matplotlib.colors as mcolors
from netCDF4 import Dataset
import os
import datetime as dt
import itertools


def get_cloud_phase_stats(phs, ncloud, label: str = "", verbose=False):
    """
    Get stats for cloud height and altitude observations
    @param phs  cloud base temperatures, for each month and year
    @param ncloud  Numpy array of year x month with number of clouds
    """
    totphs = 0
    totcld = 0
    for imo in range(np.shape(ncloud)[1]):
        sumy = 0
        ncld = 0
        for iyr in range(np.shape(ncloud)[0]):
            sumy += len(phs[imo][iyr])
            ncld += ncloud[iyr, imo]
        totcld += ncld
        totphs += sumy
        if verbose:
            print(label)
            print(sumy)
            print(ncld)
            print(sumy / ncld * 100)
            print()

    return totphs, totcld


def organize_for_hist(myvar):
    data = [[] for i in range(4)]
    # DJF
    data[0] = (
        list(itertools.chain.from_iterable(myvar[11]))
        + list(itertools.chain.from_iterable(myvar[0]))
        + list(itertools.chain.from_iterable(myvar[1]))
    )
    # MAM
    data[1] = (
        list(itertools.chain.from_iterable(myvar[2]))
        + list(itertools.chain.from_iterable(myvar[3]))
        + list(itertools.chain.from_iterable(myvar[4]))
    )
    # JJA
    data[2] = (
        list(itertools.chain.from_iterable(myvar[5]))
        + list(itertools.chain.from_iterable(myvar[6]))
        + list(itertools.chain.from_iterable(myvar[7]))
    )
    # SON
    data[3] = (
        list(itertools.chain.from_iterable(myvar[8]))
        + list(itertools.chain.from_iterable(myvar[9]))
        + list(itertools.chain.from_iterable(myvar[10]))
    )

    return data


def get_files(direc, fmt, date1, date2):
    samplefile = date1.strftime(fmt)
    allfiles = os.listdir(direc)
    allfiles.sort()
    fnames = [
        x
        for x in allfiles
        if len(x) == len(samplefile)
        and (
            dt.datetime.strptime(x, fmt) >= date1
            and dt.datetime.strptime(x, fmt) < date2
        )
    ]
    return fnames


def plotsinglepanel(
    fig,
    ax,
    panelnames,
    time,
    alt,
    val,
    ymax,
    ylabel,
    cbar,
    cmaps,
    ticks,
    ticklabels,
    cbarfrac,
    titlestr,
    labelcolors="gray",
):
    """
    @param direc  Directory to save figure to
    altitudes
    calibration
    counts
    datestring
    options
    pulse_info
    summary
    """

    ymin = 0
    fontsize = 12

    # Discrete pixels
    i = 0
    bounds = cmaps[i][0]
    mymap = get_cmap(cmaps[i][1])
    norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=mymap.N)
    pcm = ax.pcolor(time, alt, val[i].T, norm=norm, cmap=mymap)
    #                    vmin=cbar[i,0], vmax=cbar[i,1],
    cbar = fig.colorbar(pcm, ax=ax, fraction=cbarfrac, location="right")
    cbar.set_ticks(ticks[i])
    cbar.set_ticklabels(ticklabels[i], fontsize=fontsize)

    # Other formatting
    ax.set_xlim([0, 24])
    ax.set_xticks([0, 3, 6, 9, 12, 15, 18, 21, 24])
    ax.set_position((0.11, 0.15, 0.75, 0.75))
    ylim = ax.set_ylim([ymin, ymax[i]])
    ax.set_ylabel(ylabel[i], fontsize=fontsize)
    ax.text(
        0.5,
        (ylim[1] - ylim[0]) * 10 / 12 + ylim[0],
        panelnames[i],
        fontsize=fontsize,
        fontweight="bold",
        color=labelcolors,
    )
    # ax.tick_params(axis="x", rotation=rot)


def get_cmap(rgb_list, float_list=None):
    """creates and returns a color map that can be used in heat map figures.
    If float_list is not provided, colour map graduates linearly between each color in hex_list.
    If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

    Parameters
    ----------
    rgb_list: list of colors
    float_list: list of floats between 0 and 1, same length as rgb_list. Must start with 0 and end with 1.

    Returns
    ----------
    colour map"""
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(["red", "green", "blue"]):
        col_list = [
            [float_list[i], rgb_list[i][num], rgb_list[i][num]]
            for i in range(len(float_list))
        ]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap("my_cmp", segmentdata=cdict, N=256)

    return cmp


def get_column_type(fname):

    with Dataset(fname) as ncid:
        # First get the combined mask
        # 0.6 => clear (0.5-1.5)
        # 1.6 => unknown cloud (1.5-2.5)
        # 3.0 => liq (2.5-3.5)
        # 4.0 => ice (3.5-4.5)
        # combo = get_combined_mask(ncid)

        # ntime = ncid.dimensions["Time"].size

        # Stillwell match - First use Stillwell's column types
        column_type_s = ncid["ColumnType"][:].data  # stillwell's
        iliq = np.where(column_type_s == 3)[0]

        for i in iliq:
            # Is there anything classified liquid that is neither nan nor 1?
            dff = [x != 1 and not np.isnan(x) for x in ncid["CloudMask"][i, :]]
            if np.any(dff):
                idiff = np.where(dff)[0]
                print(ncid["CloudMask"][i, idiff])
                # If there is clear sky, fine
                if np.any(ncid["CloudMask"][i, idiff]) != 0:
                    print(ncid["CloudMask"][i, idiff])

        # for itime in range(ntime):
        #     if np.all(np.isnan(combo[itime, :])):
        #         column_type_s[itime] = 0
        #     elif np.any(combo[itime, :] == 3):
        #         column_type_s[itime] = 3  # liquid
        #     elif np.any(combo[itime, :] == 4):
        #         column_type_s[itime] = 4  # ice
        #     elif np.any(np.isclose(combo[itime, :], 1.6)):
        #         column_type_s[itime] = 2  # unk cloud
        #     elif np.any(np.isclose(combo[itime, :], 0.6)):
        #         column_type_s[itime] = 1  # clear
        #     else:
        #         raise ValueError("Bad mask value")

        # Now loop over time and get the column
        # for Stillwell's classification:
        #   Long_Name: Column Data Mask (
        #     0 = obscured column,
        #     1 = clear column,
        #     2 = unidentified cloud column,
        #     3 = liquid column,
        #     4 = ice column)
        # for itime in range(ntime):
        #     if np.all(np.isnan(combo[itime, :])):
        #         column_type_s[itime] = 0
        #     elif np.any(combo[itime, :] == 3):
        #         column_type_s[itime] = 3  # liquid
        #     elif np.any(combo[itime, :] == 4):
        #         column_type_s[itime] = 4  # ice
        #     elif np.any(np.isclose(combo[itime, :], 1.6)):
        #         column_type_s[itime] = 2  # unk cloud
        #     elif np.any(np.isclose(combo[itime, :], 0.6)):
        #         column_type_s[itime] = 1  # clear
        #     else:
        #         raise ValueError("Bad mask value")

        return column_type_s

        # Start from the bottom and get column type.
        # Liquid only => liquid
        # Liquid/Ice/nan => liquid
        # Liquid/Ice/clear => liquid


def get_combined_mask(ncid):
    """
    Combine the cloud mask (cloud/no cloud/ nan)
    with the phase mask (liq/ice/unknown)
    """
    # Create a combined mask:
    #  PhaseMask: nan=>CloudMask, 1(liq)=>3, 2(ice)=>4
    #  CloudMask: nan=>?, 0(no cloud()=>0.6, 1(cloud)=>1.6
    # Thus in the combined mask:
    # nan => no data
    # 0.6 => clear (0.5-1.5)
    # 1.6 => unknown cloud (1.5-2.5)
    # 3.0 => liq (2.5-3.5)
    # 4.0 => ice (3.5-4.5)
    combo = -999 * np.ones(ncid["PhaseMask"].shape)
    for i in range(ncid["PhaseMask"].shape[1]):
        combined = ncid["PhaseMask"][:, i] + 2
        # TODO: what is the following for?
        # if not np.all(ncid["CloudMask"][np.isnan(combined), i]):
        #     raise ValueError("Found cloud where NaN expected")
        inan = np.where(np.isnan(combined))[0]
        combined[inan] = ncid["CloudMask"][inan, i] + 0.6
        combo[:, i] = combined
    return combo


def get_cloud_rs(height, mask, temperature):
    """Require five pixels in a row to call it the cloud base"""

    # Get indices to bases
    ibases = np.where(mask == 1)[0]

    # Ignore the first 4 levels, which are unreliable
    ibases = ibases[ibases > 3]

    # Find the base height for the first 5 consecutive pixes

    # If only one level exists, use it
    # otherwise, pick the first of any 2 consecutive
    # otherwise if none consecutive, use the first > 5
    # or use 5 if len = 3
    if len(ibases) < 1:
        raise ValueError("Clouds indicated but not found")
    elif len(ibases) == 1:
        ibase = ibases[0]
    else:
        ibase_consec = ibases[:-1][np.diff(ibases) == 1]
        # Find first case of two pixels in a row
        if len(ibase_consec) > 0:
            ibase = ibase_consec[0]
        elif ibases[0] == 5 and len(ibases) == 2:
            ibase = ibases[1]
        else:
            ibase = ibases[0]
    base = height[ibase]
    # if np.abs(alt_km[i] - base) > 5000:
    #     print("pause")
    # basediff.append(alt_km[i] - base)

    # basetemp = ncid["Temperature"][i, ncid["Range"][:].data == base].data[0]
    basetemp = temperature[height == base]

    # cloud top
    itop = np.where(mask == 1)[0][-1]
    topalt = height[itop]
    toptemp = temperature[itop]

    # Find first clear-sky, first nan, and first cloud after cloud base
    inan = np.where(np.isnan(mask[ibase:]))[0]
    ized = np.where(mask[ibase:] == 0)[0]
    iuppercld = np.where(mask[ibase:] == 1)[0]

    upperlayer = 0
    seethru = 0

    # If we get another 1 after a zero or nan, we have an upper layer
    # and we must be able to see through the lower layer
    if (
        len(iuppercld) > 0
        and (len(inan) > 0 and np.any(iuppercld > inan[0]))
        or (len(ized) > 0 and np.any(iuppercld > ized[0]))
    ):
        upperlayer = 1
        seethru = 1

    # If we get a zero before a nan, we can see through
    if len(ized) > 0 and len(inan) > 0 and ized[0] < inan[0]:
        seethru = 1

    return base, basetemp, topalt, toptemp, seethru, upperlayer, ibase


def get_cloud_complex(height, mask, temperature):
    ibases = np.where(mask == 1)[0]

    # Ignore the first levels, which are unreliable
    ibases = ibases[ibases > 3]

    # If only one level exists, use it
    # otherwise, pick the first of any 2 consecutive
    # otherwise if none consecutive, use the first > 5
    # or use 5 if len = 3
    if len(ibases) < 1:
        raise ValueError("Clouds indicated but not found")
    elif len(ibases) == 1:
        ibase = ibases[0]
    else:
        ibase_consec = ibases[:-1][np.diff(ibases) == 1]
        # Find first case of two pixels in a row
        if len(ibase_consec) > 0:
            ibase = ibase_consec[0]
        elif ibases[0] == 5 and len(ibases) == 2:
            ibase = ibases[1]
        else:
            ibase = ibases[0]
    base = height[ibase]
    # if np.abs(alt_km[i] - base) > 5000:
    #     print("pause")
    # basediff.append(alt_km[i] - base)

    # basetemp = ncid["Temperature"][i, ncid["Range"][:].data == base].data[0]
    basetemp = temperature[height == base]

    # cloud top
    itop = np.where(mask == 1)[0][-1]
    topalt = height[itop]
    toptemp = temperature[itop]

    # Find first clear-sky, first nan, and first cloud after cloud base
    inan = np.where(np.isnan(mask[ibase:]))[0]
    ized = np.where(mask[ibase:] == 0)[0]
    iuppercld = np.where(mask[ibase:] == 1)[0]

    upperlayer = 0
    seethru = 0

    # If we get another 1 after a zero or nan, we have an upper layer
    # and we must be able to see through the lower layer
    if (
        len(iuppercld) > 0
        and (len(inan) > 0 and np.any(iuppercld > inan[0]))
        or (len(ized) > 0 and np.any(iuppercld > ized[0]))
    ):
        upperlayer = 1
        seethru = 1

    # If we get a zero before a nan, we can see through
    if len(ized) > 0 and len(inan) > 0 and ized[0] < inan[0]:
        seethru = 1

    if ibase <= 3:
        print(ibase)
    return base, basetemp, topalt, toptemp, seethru, upperlayer, ibase


def get_cloud_lowest_after_4th_range(height, mask, temperature):
    """Get the cloud as the lowest base after the first 4 ranges"""

    ibases = np.where(mask == 1)[0]

    # Ignore the first levels, which are unreliable
    ibases = ibases[ibases > 3]
    ibase = ibases[0]
    base = height[ibase]
    basetemp = temperature[height == base]

    # cloud top
    itop = np.where(mask == 1)[0][-1]
    topalt = height[itop]
    toptemp = temperature[itop]

    # Find first clear-sky, first nan, and first cloud after cloud base
    inan = np.where(np.isnan(mask[ibase:]))[0]
    ized = np.where(mask[ibase:] == 0)[0]
    iuppercld = np.where(mask[ibase:] == 1)[0]

    upperlayer = 0
    seethru = 0

    # If we get another 1 after a zero or nan, we have an upper layer
    # and we must be able to see through the lower layer
    if (
        len(iuppercld) > 0
        and (len(inan) > 0 and np.any(iuppercld > inan[0]))
        or (len(ized) > 0 and np.any(iuppercld > ized[0]))
    ):
        upperlayer = 1
        seethru = 1

    # If we get a zero before a nan, we can see through
    if len(ized) > 0 and len(inan) > 0 and ized[0] < inan[0]:
        seethru = 1

    if ibase <= 3:
        raise ValueError("We should not get a base index <= 3")
    return base, basetemp, topalt, toptemp, seethru, upperlayer, ibase


def get_cloud_lowest_after_2nd_range(height, mask, temperature):
    """Get the cloud as the lowest base after the first 4 ranges"""

    ibases = np.where(mask == 1)[0]

    # Ignore the first levels, which are unreliable
    ibases = ibases[ibases > 1]
    ibase = ibases[0]
    base = height[ibase]
    basetemp = temperature[height == base]

    # cloud top
    itop = np.where(mask == 1)[0][-1]
    topalt = height[itop]
    toptemp = temperature[itop]

    # Find first clear-sky, first nan, and first cloud after cloud base
    inan = np.where(np.isnan(mask[ibase:]))[0]
    ized = np.where(mask[ibase:] == 0)[0]
    iuppercld = np.where(mask[ibase:] == 1)[0]

    upperlayer = 0
    seethru = 0

    # If we get another 1 after a zero or nan, we have an upper layer
    # and we must be able to see through the lower layer
    if (
        len(iuppercld) > 0
        and (len(inan) > 0 and np.any(iuppercld > inan[0]))
        or (len(ized) > 0 and np.any(iuppercld > ized[0]))
    ):
        upperlayer = 1
        seethru = 1

    # If we get a zero before a nan, we can see through
    if len(ized) > 0 and len(inan) > 0 and ized[0] < inan[0]:
        seethru = 1

    if ibase <= 1:
        raise ValueError("We should not get a base index <= 3")
    return base, basetemp, topalt, toptemp, seethru, upperlayer, ibase


# def get_layers(height, mask, temperature):
#     """Get the various cloud layers"""

#     foga = 0 # Fog in first two layers
#     fogb = 0 # Fog in third layer
#     fogc = 0 # Fog in first 5 layers

#     return fog


# def get_mo_data(inp):
#     # Loop over years and get months
#     data = [[] for i in range(len(months))]
#     for iyr in range(len(years)):
#         for imo in range(len(months)):
#             data[imo] += inp[iyr][imo]
#     return data


# def plot_day(fname):
#     fig, axno = plt.subplots(1, 1, constrained_layout=True)

#     with Dataset(fname) as ncid:
#         # dict_keys(['DataTime', 'Range', 'CloudMask', 'PhaseMask',
#         #   'CloudBaseAltitude', 'CloudBaseTemperature', 'CloudTopAltitude',
#         #   'CloudTopTemperature', 'ColumnType', 'Temperature', 'CountRate',
#         #   'Depolarization', 'BackscatterRatio'])
#         panelnames = ["a)"]  # ["Data Mask "]

#         #% This colortable makes the following values in the data mask a single color
#         #% Take the log10 of the mask, you will get data that is colored as:
#         #% Good data = green, clipped data = red, missing = blue, SNR filtered = white
#         #% and density filtered = yellow. If multiple filters are tripped, the
#         #% highest number is plotted
#         colortable = [
#             [0, 1, 1],  # Clear
#             [0.75, 0.75, 0.75],  # Cloud
#             [0, 0, 1],  # Liquid
#             [1, 0, 0],  # Ice
#         ]

#         alt = ncid["Range"][:] / 1000
#         timestamps = ncid["DataTime"][:]

#         combo = get_combined_mask(ncid)

#         val = [combo]

#     ymax = [2.5]
#     ylabel = ["Altitude (km)"]

#     cbarbnds = np.array([[0, 10]])
#     cmaps = [[np.arange(0.5, 5, 1), colortable]]
#     cticks = [cmaps[0][0][:-1] + 0.5]
#     cticklabels = [["Clear", "Cloud", "Liquid", "Ice"]]

#     cbarfrac = 0.03

#     # Create the figure
#     plotsinglepanel(
#         fig,
#         axno,
#         panelnames,
#         timestamps,
#         alt,
#         val,
#         ymax,
#         ylabel,
#         cbarbnds,
#         cmaps,
#         cticks,
#         cticklabels,
#         cbarfrac,
#         titlestr,
#     )

#     plt.text(3, 4, "a)")

#     # Save the figure
#     # outfig = f"KGI_{datestring}_PolarizationPrelim.png"
#     # fig.savefig(lidar_dir + outfig)
#     # plt.close(fig=fig)


# def get_mean_and_std(inp):
#     data = np.nan + np.zeros([12, 2])
#     box = get_mo_data(inp)
#     for imo in range(len(months)):
#         data[imo, 0] = np.nanmean(np.array(box[imo]))
#         data[imo, 1] = np.nanstd(np.array(box[imo]))
#     return data


# def sample_plot():  # fig, axno, lidar_dir):

#     meas_dir = "/Users/prowe/Sync/measurements/Escudero/mpl/"
#     lidar_dir = meas_dir + "reprocessing_stillwell/NetCDF/"

#     datestring = "20220516"
#     datestring = "20211225"

#     # titlestr = "King George Island MPL Data Processed " + datestring

#     fname = f"KGI_MPLData_{datestring}.cdf"
#     plot_day(lidar_dir + fname)
