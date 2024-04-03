#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 11:32:46 2023

@author: prowe
"""

# Built-in modules
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import matplotlib.colors as mcolors
from netCDF4 import Dataset


def plot_measured_and_clear(
    swd_date, swd, swd_clr, lwd_date, lwd, lwd_clr_date, lwd_clr, savef, fname
):
    rot = 0
    datefmt = "%d"
    plt.figure(num=1)
    ax1 = plt.subplot(211)
    ax3 = plt.subplot(212, sharex=ax1)
    wid = 0.85
    hit = 0.4
    ax1.set_position([0.12, 0.57, wid, hit])
    ax3.set_position([0.12, 0.10, wid, hit])

    ax1.plot(swd_date, swd, "r", label="Measured")
    ax1.plot(swd_date, swd_clr, "k--", label="Clear")
    ax1.legend()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
    ax1.tick_params(axis="x", rotation=rot)
    ax1.set_ylabel("SWD (W/m$^{2}$)")

    ax3.plot(lwd_date, lwd, color="blue", label="Measured")
    ax3.plot(lwd_clr_date, lwd_clr, "o:", color="blue", label="Clear")
    ax3.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_xlabel("Day of 2022/05")
    ax3.set_ylabel("LWD (W/m$^{2}$)")
    ax3.legend(loc="upper left")

    if savef:
        plt.savefig(fname)


def plot_meas_clear_forcings(
    swd_date,
    swd,
    swd_clr,
    swd_force,
    lwd_date,
    lwd,
    lwd_clr_date,
    lwd_clr,
    lwd_force,
    savef,
    fname,
    start_date,
    final_date,
):

    rot = 40
    plt.figure(num=2, figsize=[6.5, 5.6])
    plt.clf()
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312, sharex=ax1)
    ax3 = plt.subplot(313, sharex=ax1)
    wid = 0.84
    ax1.set_position([0.13, 0.67, wid, 0.30])
    ax2.set_position([0.13, 0.45, wid, 0.18])
    ax3.set_position([0.13, 0.09, wid, 0.33])

    ax1.plot(swd_date, swd, color="red", label="Meas")
    ax1.plot(swd_date, swd_clr, "k--", label="Clear")
    ax1.set_ylim([-2, 120])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax1.tick_params(axis="x", rotation=rot, labelcolor="white")
    ax1.set_xlim(start_date, final_date)
    ax1.set_ylabel("SWD (W/m$^{2}$)")
    ax1.legend()

    ax2.plot(lwd_date, lwd, color="blue", label="Measured")
    ax2.plot(lwd_clr_date, lwd_clr, "o:", color="blue", label="Clear")
    ax2.set_ylim([150, 350])
    ax2.legend()
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax2.tick_params(axis="x", rotation=rot, labelcolor="white")
    ax2.set_xlabel("Day and Hour in 2022")
    ax2.set_ylabel("LWD (W/m$^{2}$)")

    tot_force = swd_force + lwd_force
    ax3.plot(lwd_date, lwd_force, color="blue", linewidth=2, label="LWD")
    ax3.plot(swd_date, swd_force, color="red", label="SWD")
    ax3.plot(swd_date, tot_force, color="orange", label="Total")
    ax3.plot(lwd_date, np.zeros(len(lwd_date)), "k:")
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
    ax3.set_ylim([-100, 150])
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_ylabel("Forcing (W/m$^{2}$)")
    ax3.legend(loc=[0.01, 0.62])

    if savef:
        plt.savefig(fname)


def plot_forcings(
    swd_date,
    swd,
    swd_clr,
    swd_force,
    lwd_date,
    lwd,
    lwd_clr_date,
    lwd_clr,
    lwd_force,
    start_date,
    final_date,
    datefmt: str = "%d",
    xlab: str = "Day of 2022/05",
    savef: bool = False,
    fname: str = "",
):

    rot = 0
    plt.figure(num=3)  # , figsize=[6.5, 5.6])
    plt.clf()
    ax3 = plt.subplot(111)
    ax3.set_position([0.13, 0.11, 0.84, 0.84])

    tot_force = swd_force + lwd_force
    ax3.plot(lwd_date, lwd_force, color="orange", linewidth=3, label="LWD")
    ax3.plot(swd_date, swd_force, color="red", label="SWD")
    ax3.plot(swd_date, tot_force, color="blue", label="Total")
    ax3.plot(lwd_date, np.zeros(len(lwd_date)), "k:")
    ax3.legend()
    ax3.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_xlabel(xlab)

    ax3.set_xlim(start_date, final_date)
    ax3.set_ylim([-75, 150])
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_ylabel("Forcing (W/m$^{2}$)")

    if savef:
        plt.savefig(fname)


def plot_forcings_and_lidar(
    swd_date,
    swd,
    swd_clr,
    swd_force,
    lwd_date,
    lwd,
    lwd_clr_date,
    lwd_clr,
    lwd_force,
    start_date,
    final_date,
    lidar_dir: str,
    datefmt: str = "%H",
    xlab: str = "Hor on 2022/05/16",
    savef: bool = False,
    fname: str = "",
):
    """
    lidar_dir = f"{params.MEAS_DIR}Escudero/mpl/reprocessing/ResultNetCDF/"
    """

    fig, axs = plt.subplots(2, 1, constrained_layout=True, num=4)
    # plt.figure(num=3)  # , figsize=[6.5, 5.6])
    # plt.clf()
    # figsize=(plotwidth, plotheight)
    rot = 0

    # Plot the miniMPL data above
    run_plotmultipanel(fig, axs[0], lidar_dir)

    ax3 = axs[1]

    tot_force = swd_force + lwd_force
    ax3.plot(lwd_date, lwd_force, color="orange", linewidth=3, label="LWD")
    ax3.plot(swd_date, swd_force, color="red", label="SWD")
    ax3.plot(swd_date, tot_force, color="blue", label="Total")
    ax3.plot(lwd_date, np.zeros(len(lwd_date)), "k:")
    ax3.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_xlabel(xlab)

    ax3.set_xlim(start_date, final_date)
    ax3.set_ylim([-75, 150])
    ax3.tick_params(axis="x", rotation=rot)
    ax3.set_ylabel("Forcing (W/m$^{2}$)")
    # ax3.legend(bbox_to_anchor=(1., 1), loc="upper left", borderaxespad=0.0)
    ax3.legend(loc="upper center", bbox_to_anchor=(1.093, 0.8))
    ax3.set_position((0.11, 0.1, 0.75, 0.43))
    # ax3.legend(loc=[1.3, 0.62])

    if savef:
        plt.savefig(fname)


def plotmultipanel(
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
    ax.set_position((0.11, 0.6, 0.75, 0.35))
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


def run_plotmultipanel(fig, axno, lidar_dir):

    datestring = "20220516"
    fname = f"KGI_MPLDataV2{datestring}.nc"
    titlestr = "King George Island MPL Data Processed " + datestring

    with Dataset(lidar_dir + fname) as ncid:
        # dict_keys(['DataTime', 'Range', 'CloudMask', 'PhaseMask',
        #   'CloudBaseAltitude', 'CloudBaseTemperature', 'CloudTopAltitude',
        #   'CloudTopTemperature', 'ColumnType', 'Temperature', 'CountRate',
        #   'Depolarization', 'BackscatterRatio'])
        panelnames = ["a)"]  # ["Data Mask "]

        #% This colortable makes the following values in the data mask a single color
        #% Take the log10 of the mask, you will get data that is colored as:
        #% Good data = green, clipped data = red, missing = blue, SNR filtered = white
        #% and density filtered = yellow. If multiple filters are tripped, the
        #% highest number is plotted
        colortable = [
            [0, 1, 1],  # Clear
            [0.75, 0.75, 0.75],  # Cloud
            [0, 0, 1],  # Liquid
            [1, 0, 0],  # Ice
        ]

        alt = ncid["Range"][:] / 1000
        timestamps = ncid["DataTime"][:]

        # Create a combined mask
        combo = -999 * np.ones(ncid["PhaseMask"].shape)
        for i in range(ncid["PhaseMask"].shape[1]):
            combined = ncid["PhaseMask"][:, i] + 2
            # TODO: what is the following for?
            # if not np.all(ncid["CloudMask"][np.isnan(combined), i]):
            #     raise ValueError("Found cloud where NaN expected")
            inan = np.where(np.isnan(combined))[0]
            combined[inan] = ncid["CloudMask"][inan, i] + 0.6
            combo[:, i] = combined

        val = [combo]

    ymax = [2.5]
    ylabel = ["Altitude (km)"]

    cbarbnds = np.array([[0, 10]])
    cmaps = [[np.arange(0.5, 5, 1), colortable]]
    cticks = [cmaps[0][0][:-1] + 0.5]
    cticklabels = [["Clear", "Cloud", "Liquid", "Ice"]]

    cbarfrac = 0.03

    # Create the figure
    plotmultipanel(
        fig,
        axno,
        panelnames,
        timestamps,
        alt,
        val,
        ymax,
        ylabel,
        cbarbnds,
        cmaps,
        cticks,
        cticklabels,
        cbarfrac,
        titlestr,
    )

    plt.text(3, 4, "a)")

    # Save the figure
    # outfig = f"KGI_{datestring}_PolarizationPrelim.png"
    # fig.savefig(lidar_dir + outfig)
    # plt.close(fig=fig)


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
