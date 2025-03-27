#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 14:19:28 2023

@author: prowe

Plot the cloud mask from the mpl
"""

# Built-in modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates
from netCDF4 import Dataset
import datetime as dt
import matplotlib.colors as mcolors

# Import params containing personal root directory
# (this file must be created by each user)
from antarc import params


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


def plotmultipanel(
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
    labelcolors=["gray", "gray", "gray", "gray", "gray"],
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

    # Set the output figure resolution:
    # plt.rcParams["figure.dpi"] = 300
    # plt.rcParams["savefig.dpi"] = 300

    ymin = 0
    npanels = len(val)
    timespan = 0.01 * (time[-1] - time[0]).total_seconds()
    panelx = time[0] + dt.timedelta(seconds=timespan)

    ### Setting parameters for figure and panel size
    plotheight = 12.0  # Total plot height
    plotwidth = 8.0  # Total plot width
    leftmargin = 0.08  # Margin on the left side of all subplots
    botmarg = 0.06  # Margin on the bottom side of all subplots
    fontsize = 12
    panelwidth = 0.86
    vertpad = 0.05
    toppad = 0.03
    panelheight = (1 - (botmarg + (npanels - 1) * vertpad + toppad)) / npanels

    panelbottom = []
    for i in range(npanels):
        fac = npanels - 1 - i
        panelbottom.append(botmarg + fac * vertpad + fac * panelheight)

    fig, axs = plt.subplots(
        npanels,
        1,
        num=1,
        clear=True,
        figsize=(plotwidth, plotheight),
        constrained_layout=True,
    )

    for i in range(npanels):
        ax = axs[i]
        ax.set_position([leftmargin, panelbottom[i], panelwidth, panelheight])

        # Create the pcolor plot
        if not isinstance(val[i], str) and isinstance(cmaps[i], str):
            # Contour plot
            pcm = ax.pcolor(
                time,
                alt,
                val[i].T,
                shading="auto",
                vmin=cbar[i, 0],
                vmax=cbar[i, 1],
                cmap=plt.get_cmap(cmaps[i]),
            )
            fig.colorbar(pcm, ax=ax, fraction=cbarfrac, location="right")
            ax.xaxis.set_minor_locator(mdates.HourLocator(interval=6))

        elif not isinstance(val[i], str):
            if np.shape(val[i])[1] == 1:
                # No altitude data; time series only - assume cloud column mask
                col = val[i][:, 0]
                clrs = cmaps[i]
                # ax.plot(time[col == -1], col[col == -1], ".", color=clrs[0])
                ax.plot(time[col == 0], col[col == 0], ".", color=clrs[1])
                ax.plot(time[col == 1], col[col == 1], ".", color=clrs[2])
                ax.plot(time[col == 2], col[col == 2], ".", color=clrs[3])
                ax.plot(time[col == 3], col[col == 3], ".", color=clrs[4])
                ax.plot(time[col == 4], col[col == 4], ".", color=clrs[5])
                ax.set_yticks(ticks[i])
                ax.set_yticklabels(ticklabels[i])
                ax.set_position(
                    [
                        leftmargin,
                        panelbottom[i],
                        panelwidth * 0.9,
                        panelheight,
                    ]
                )
                ymin = -0.5
            else:
                # Discrete pixels
                bounds = cmaps[i][0]
                mymap = get_cmap(cmaps[i][1])
                norm = colors.BoundaryNorm(boundaries=bounds, ncolors=mymap.N)
                pcm = ax.pcolor(time, alt, val[i].T, norm=norm, cmap=mymap)
                #                    vmin=cbar[i,0], vmax=cbar[i,1],
                cbar = fig.colorbar(
                    pcm, ax=ax, fraction=cbarfrac, location="right"
                )
                cbar.set_ticks(ticks[i])
                cbar.set_ticklabels(ticklabels[i], fontsize=fontsize)

        # Major ticks every half year, minor ticks every month,
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))

        # Title
        if i == 0:
            ax.set_title(titlestr, fontweight="bold", fontsize=fontsize)

        # Other formatting
        ax.set_xlim([time[0] - dt.timedelta(minutes=20), time[-1]])
        ylim = ax.set_ylim([ymin, ymax[i]])
        ax.set_ylabel(ylabel[i], fontsize=fontsize)
        ax.text(
            panelx,
            (ylim[1] - ylim[0]) * 10.8 / 12 + ylim[0],
            panelnames[i],
            fontsize=fontsize,
            fontweight="bold",
            color=labelcolors[i],
        )

        # Xlabel
        if i == npanels - 1:
            ax.set_xlabel("Month/Day in 2022 [UTC] ")

    fmt = {
        "plotheight": plotheight,
        "plotwidth": plotwidth,
        "leftmargin": leftmargin,
        "bottommargin": botmarg,
        "fontsize": fontsize,
        "panelwidth": panelwidth,
        "verticalpad": vertpad,
        "toppad": toppad,
        "panelheight": panelheight,
        "panelbottom": panelbottom,
    }

    return fig, axs, fmt


direc = f"{params.MEAS_DIR}Escudero/mpl/reprocessing_stillwell/v3/"
SAVE_FIGURE = False


datestring = "20220516"
datestrings = ["20221228", "20221229", "20221230", "20221231"]

panelnames = [
    # "Relative Backscatter [arb]  ",
    "a) Depolarization  ",
    "b) Backscatter Ratio  ",
    "c) Cloud Mask ",
    "d) Column Cloud Mask ",
]

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
ymax = [10.0, 10.0, 10.0, 5.6]
ylabel = [
    "Altitude (km)",
    "Altitude (km)",
    "Altitude (km)",
    "",
]

# cbarbnds = np.array([[-2, 6], [0, 0.4], [0, 10], [0.5, 4.5]])
cbarbnds = np.array([[0, 0.4], [0, 10], [0.5, 4.5]])
cmaps = [
    "viridis",
    "viridis",
    [np.arange(0.5, 5, 1), colortable],
    ["brown", "black", "cyan", "gray", "blue", "red"],
]

cticks = [
    "none",
    "none",
    cmaps[-2][0][:-1] + 0.5,
    [0, 1, 2, 3, 4],
]

cticklabels = [
    "none",
    "none",
    ["Clear", "Cloud", "Liquid", "Ice"],
    ["Obsc.", "Clear", "Cloud", "Liquid", "Ice"],
]

cbarfrac = 0.05

timestamps = []

firsttime = True

for datestring in datestrings:
    # fname = f"KGI_MPLDataV2{datestring}.nc"
    fname = f"KGI_MPLData_{datestring}.cdf"
    dtime = dt.datetime.strptime(datestring, "%Y%m%d")

    with Dataset(direc + fname) as ncid:
        # dict_keys(['DataTime', 'Range', 'CloudMask', 'PhaseMask',
        #   'CloudBaseAltitude', 'CloudBaseTemperature', 'CloudTopAltitude',
        #   'CloudTopTemperature', 'ColumnType', 'Temperature', 'CountRate',
        #   'Depolarization', 'BackscatterRatio'])

        # def plot_mask(
        #     direc,
        #     altitudes,
        #     counts,
        #     datestring,
        #     mask,
        #     options,
        #     polarization,
        #     pulse_info,
        #     summary,
        # ):
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

        # first = np.where(~np.isnan(altitudes.Integrated[:, 0]))[0][0]
        # alt = altitudes.Integrated[first, :] / 1e3
        alt = ncid["Range"][:] / 1000
        # timestamps = pulse_info.TimeStamps[:, 0]
        seconds = ncid["DataTime"][:].data * 3600
        timestamps += [dtime + dt.timedelta(seconds=float(x)) for x in seconds]

        # Create a combined mask
        combo = -999 * np.ones(ncid["PhaseMask"].shape)
        for i in range(ncid["PhaseMask"].shape[1]):
            combined = ncid["PhaseMask"][:, i] + 2
            # if not np.all(ncid["CloudMask"][np.isnan(combined), i]):
            #     raise ValueError("Found cloud where NaN expected")
            inan = np.where(np.isnan(combined))[0]
            combined[inan] = ncid["CloudMask"][inan, i] + 0.6
            combo[:, i] = combined

        # Pad mask.CombinedMask with nan so dim 1 is same as alt
        # cols = np.shape(mask.CombinedMask)[1]
        # rowdiff = len(alt) - np.shape(mask.CombinedMask)[0]
        # combined_mask = np.vstack(
        #     [mask.CombinedMask, np.nan * np.ones((rowdiff, cols))]
        # )
        # val = [
        #     np.real(np.log10(counts.SpeckleFiltered2[-1, 0])).T,
        #     polarization.Depol[0, 0].T,
        #     polarization.SigmaR[0, 0],
        #     combined_mask,
        #     mask.ColumnType,
        # ]

        if firsttime:
            val = [
                ncid["Depolarization"][:],
                ncid["BackscatterRatio"][:],
                combo,
                np.reshape(ncid["ColumnType"][:], (-1, 1)),
            ]
            firsttime = False
        else:
            val[0] = np.vstack([val[0], ncid["Depolarization"][:]])
            val[1] = np.vstack([val[1], ncid["BackscatterRatio"][:]])
            val[2] = np.vstack([val[2], combo])
            val[3] = np.vstack(
                [val[3], np.reshape(ncid["ColumnType"][:], (-1, 1))]
            )


# Create the figure
fig, axs, _ = plotmultipanel(
    panelnames,
    np.array(timestamps),
    alt,
    val,
    ymax,
    ylabel,
    cbarbnds,
    cmaps,
    cticks,
    cticklabels,
    cbarfrac,
    "",
    ["black", "black", "black", "black", "black"],
)

# Save the figure
if SAVE_FIGURE:
    outfig = f"KGI_lidar_{datestring}.png"
    fig.savefig(direc + "figures/" + outfig, format="png", dpi=600)
