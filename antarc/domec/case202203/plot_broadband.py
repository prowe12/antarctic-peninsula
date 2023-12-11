#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 09:53:04 2022

@author: prowe

Purpose: Plot surface energy balance components at DomeC from
         observations, and optionally add PWRF, ERA5, and SNOWPACK
"""

import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

# from read_raw import RadRaw
from broadband_reader import Rad
from broadband_plotter import plotnet  # , plotnet2
from seb_reader import Pwrf, Era5, Snowpack

# from era5_reader import Era5
# from pwrf_reader import Pwrf


def dup(arr: np.ndarray):
    """
    Duplicate x by vertically stacking x on top of x
    @param x  A numpy array of values
    """
    return np.reshape((np.vstack([arr, arr])).T, -1)


def dupday(arr: np.ndarray):
    """
    Duplicate x by vertically stacking x on top of x,
    then removing the first value, and then appending
    the incremented last value
    @param x  A numpy array of values
    """
    arr = dup(arr)[1:]
    return np.hstack([arr, arr[-1] + 1])


def plotnet_obs():
    """Plot the total net, with daily aves"""

    _, axs_net = plt.subplots(2, 1, figsize=(6, 6))
    ax = axs_net[0]
    ax.plot(obs.dtime, tot, label="Observed")
    ax.plot(obs.dtime, np.zeros(len(obs.dtime)), "k:")

    # dayave = dupday(obs.dayofmonth)
    ax = axs_net[1]
    ax.plot(dayave, dup(totave), label="Observed", linewidth=3)

    ax = axs_net[0]
    ax.set_ylim([-80, 120])
    ax.legend()
    xlim = [dt.datetime(2022, 3, 12), dt.datetime(2022, 3, 26)]
    # xlim = [dt.datetime(2022, 3, 16), dt.datetime(2022, 3, 20, 12, 0, 0)]
    ax.set_xlim(xlim)
    # Xtick = day of month
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    ax.set_ylabel("Net radiation (W/m$^{2}$)")
    axs_net[0].set_position([0.12, 0.56, 0.83, 0.42])

    # dayave = dupday(obs.dayofmonth)
    ax = axs_net[1]
    ax.plot(dayave, np.zeros(len(dayave)), "k:")
    ax.legend()
    ax.set_xlim([12, 26])
    ax.set_ylim([-80, 120])
    ax.set_xlabel("Day on 2022/03")
    ax.set_ylabel("Daily ave net radiation (W/m$^{2}$)")
    ax.set_position([0.12, 0.082, 0.83, 0.42])
    ax.set_xticks(list(range(12, 27)))

    return axs_net


def plotnet_era5(axs):
    """
    Add ERA5 data to the existing plotof total net from observations
    @param axs  The axes to add to, consisting of two subplots
    """
    ax = axs[0]
    ax.plot(era5.dtime, etot, "c", label="ERA5")
    ax.legend()

    ax = axs[1]
    ax.plot(dupday(eday), dup(etotave), "c", label="ERA5")
    ax.legend()


def plotnet_pwrf(axs):
    """
    Add Polar WRF data to the existing plotof total net from observations
    @param axs  The axes to add to, consisting of two subplots
    """
    ax = axs[0]
    ax.plot(pwrf.dtime, ptot, color="orange", label="PWRF")
    ax.legend()

    ax = axs[1]
    ax.plot(dupday(pday), dup(ptotave), color="orange", label="PWRF")
    ax.legend()


def plotnet_snowpack(axs):
    """
    Add SNOWPACK data to the existing plotof total net from observations
    @param axs  The axes to add to, consisting of two subplots
    """
    ax = axs[0]
    ax.plot(spack.dtime, spacktot, "m", label="SNOWPACK")
    ax.text(spack.dtime[273], 97, "a)", fontsize=16)
    ax.legend()

    ax = axs[1]
    ax.plot(dupday(spackday), dup(spacktotave), "m", label="SNOWPACK")
    ax.text(12.4, 97, "b)", fontsize=16)
    ax.legend()


def plot_each():
    """Plot each of the surface energy balance components from obs"""
    _, axs = plt.subplots(3, 2, figsize=(10, 10))
    axs = [axs[0][0], axs[0][1], axs[1][0], axs[1][1], axs[2][0], axs[2][1]]
    cols = ["SWD", "LWD", "SWU", "LWU", "SWnet", "LWnet"]

    for icol, ax in enumerate(axs):
        locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)

        seb_component = obs.data[cols[icol]][: obs.ndata]
        ax.plot(obs.dtime, seb_component, color=crs[icol], label="Observed")
        if cols[icol] == "SWD":
            label = "SWD (W/m$^2$)"
        elif cols[icol] == "SWU":
            label = "SWU (W/m$^2$)"
        elif cols[icol] == "LWD":
            label = "LWD (W/m$^2$)"
        elif cols[icol] == "LWU":
            label = "LWU (W/m$^2$)"
        elif cols[icol] == "SWnet":
            label = "SW Net (W/m$^2$)"
        elif cols[icol] == "LWnet":
            label = "LW Net (W/m$^2$)"
        ax.set_ylabel(label)

    # Restrict the x-axis
    xlim = [dt.datetime(2022, 3, 12), dt.datetime(2022, 3, 26)]
    # xlim = [dt.datetime(2022, 3, 16), dt.datetime(2022, 3, 20, 12, 0, 0)]
    axs[0].set_xlim(xlim)
    axs[0].set_ylim([-20, 500])
    axs[1].set_xlim(xlim)
    # axs[1].set_ylim([120, 310])
    axs[1].set_ylim([-20, 500])
    axs[2].set_xlim(xlim)
    axs[2].set_ylim([-20, 500])
    axs[3].set_xlim(xlim)
    # axs[3].set_ylim([120, 310])
    axs[3].set_ylim([-20, 500])
    axs[4].set_xlim(xlim)
    axs[4].set_ylim([-75, 125])
    axs[5].set_xlim(xlim)
    axs[5].set_ylim([-75, 125])

    # Set the xtick
    # axs[0].xaxis.set_major_locator(mdates.DayLocator())
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    # axs[2].xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    # axs[3].xaxis.set_major_formatter(mdates.DateFormatter("%d"))

    # Set the positions
    lft = 0.07
    wid = 0.41
    bot = 0.06
    mid = 0.38
    top = 0.7
    hit = 0.29
    axs[0].set_position([lft, top, wid, hit])
    axs[1].set_position([0.56, top, wid, hit])
    axs[2].set_position([lft, mid, wid, hit])
    axs[3].set_position([0.56, mid, wid, hit])
    axs[4].set_position([lft, bot, wid, hit])
    axs[5].set_position([0.56, bot, wid, hit])

    return axs


def add_pwrf(axs):
    """
    Add the PWRF components
    @param axs  The axes to add to
    """
    axs[0].plot(pwrf.dtime, pwrf.swd, ".", color="orange", label="PWRF")
    axs[1].plot(pwrf.dtime, pwrf.lwd, ".", color="orange", label="PWRF LWD")
    axs[2].plot(pwrf.dtime, pwrf.swu, ".", color="orange", label="PWRF LWU")
    axs[3].plot(pwrf.dtime, pwrf.lwu, ".", color="orange", label="PWRF SWU")
    axs[4].plot(pwrf.dtime, pwrf.swnet, ".", color="orange", label="PWRF")
    axs[5].plot(pwrf.dtime, pwrf.lwnet, ".", color="orange", label="PWRF")

    axs[0].legend()
    axs[1].legend()


def add_era5(axs):
    """
    Add the ERA5 components
    @param axs  The axes to add to
    """
    axs[0].plot(era5.dtime, era5.swd, "c", label="ERA5")
    axs[1].plot(era5.dtime, era5.lwd, "c", label="ERA5")
    axs[2].plot(era5.dtime, era5.swu, "c", label="ERA5")
    axs[3].plot(era5.dtime, era5.lwu, "c", label="ERA5")
    axs[4].plot(era5.dtime, era5.swnet, "c", label="ERA5")
    axs[5].plot(era5.dtime, era5.lwnet, "c", label="ERA5")

    zeds = np.zeros(len(era5.dtime))
    for i in range(6):
        axs[i].plot(era5.dtime, zeds, "k:")

    axs[0].legend()
    axs[1].legend()


def add_snowpack(axs):
    """
    Add the snowpack components
    @param axs  The axes to add to
    """
    axs[0].plot(spack.dtime, spack.swd, "m--", label="SNOWPACK")
    axs[1].plot(spack.dtime, spack.lwd, "m--", label="SNOWPACK")
    axs[2].plot(spack.dtime, spack.swu, "m--", label="SNOWPACK")
    axs[3].plot(spack.dtime, spack.lwu, "m--", label="SNOWPACK")
    axs[4].plot(spack.dtime, spack.swnet, "m--", label="SNOWPACK")
    axs[5].plot(spack.dtime, spack.lwnet, "m--", label="SNOWPACK")

    # Add labels
    ind = 275
    fsz = 16
    axs[0].text(spack.dtime[ind], 440, "a)", fontsize=fsz)
    axs[1].text(spack.dtime[ind], 440, "b)", fontsize=fsz)
    axs[2].text(spack.dtime[ind], 440, "c)", fontsize=fsz)
    axs[3].text(spack.dtime[ind], 440, "d)", fontsize=fsz)
    axs[4].text(spack.dtime[ind], -50, "e)", fontsize=fsz)
    axs[5].text(spack.dtime[ind], 98, "f)", fontsize=fsz)

    axs[0].legend()
    axs[1].legend()


# Inputs
main_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/EastAntarc_202203/"
file = main_dir + "observations/DomeC_broadband/202203.domeC_mod.bsrn"
era5_file = main_dir + "ERA5/ERA5_SEB_Concordia_Mar_12_25_hourly.nc"
pwrf_file = main_dir + "PWRF/PWRF_D02_Concordia_20220316_20_LWSW_hrly.nc"
spack_file = main_dir + "snowpack/DOMEC_DFLT_202203.txt"

# Which figures to plot
plotqc = False
ploteach = False
plotallzoom = False
ploteach_with_era5_pwrf_snowpack = True  # Fig. S5
plotnet_with_era5_pwrf_snowpack = False  # Fig. S6
plot_misc = False
# plotall = True

# Get the observations, including net radiations, and daily average
obs = Rad(file)
obs.clean()
obs.getnet()
tot = obs.data["SWnet"][: obs.ndata] + obs.data["LWnet"][: obs.ndata]
totave = obs.swdave + obs.lwdave - obs.swuave - obs.lwuave
dayave = dupday(obs.dayofmonth)


pwrf = Pwrf(pwrf_file)
ptot = pwrf.getnet()
pday, ptotave = pwrf.getnet_dailyave()

era5 = Era5(era5_file)
etot = era5.getnet()
eday, etotave = era5.getnet_dailyave()

spack = Snowpack(spack_file)
spacktot = spack.getnet()
spackday, spacktotave = spack.getnet_dailyave()


if plotqc:
    obs.plotqc()


if ploteach:
    crs = ["orange", "blue", "black", "cyan", "blue", "orange"]
    plot_each()

    fname = "fig_s5_broadband_with_era5_pwrf_snowpack.pdf"
    plt.savefig(fname, format="pdf")


if ploteach_with_era5_pwrf_snowpack:
    # Plot the Dome C components
    crs = ["blue", "blue", "blue", "blue", "blue", "blue"]
    axs_each = plot_each()

    # Add ERA5, PWRF, and SNOWPACK:
    add_era5(axs_each)
    add_pwrf(axs_each)
    add_snowpack(axs_each)


if plotnet_with_era5_pwrf_snowpack:
    axes_net = plotnet_obs()
    plotnet_era5(axes_net)
    plotnet_pwrf(axes_net)
    plotnet_snowpack(axes_net)
    # plotnet_finish(axs_net)

    fname = "fig_s6_net_with_era5_pwrf_snowpack.pdf"
    # plt.savefig(fname, format="pdf")


# if plotall:
#     davetime = np.array(
#         [
#             dt.datetime(obs.dtime[0].year, obs.dtime[0].month, d)
#             for d in obs.dayofmonth
#         ]
#     )

#     # axs1 = plotnet(4, np.arange(0, len(dtime)), 16)
#     axs2 = plotnet2(
#         6, np.arange(0, len(obs.dtime)), np.arange(len(davetime)), 16
#     )

#     # plotaves(axs1, np.arange(len(davetime)))
#     # plotaves2(6, np.arange(len(davetime)), 16)


if plotallzoom:
    loday = 12
    hiday = 24

    i1 = np.where(np.array([x.day for x in obs.dtime]) >= loday)[0][0]
    i2 = np.where(np.array([x.day for x in obs.dtime]) <= hiday)[0][-1]
    axs1 = plotnet(obs.dtime, obs.data, 5, np.arange(i1, i2 + 1), 5)

    inds = [
        i
        for i, m in enumerate(obs.dayofmonth)
        if loday <= m <= hiday + 1  # m >= loday and m <= hiday + 1
    ]
    davetime = np.array(
        [
            dt.datetime(obs.dtime[0].year, obs.dtime[0].month, d)
            for d in obs.dayofmonth
        ]
    )
    # plotaves(axs2, np.array(inds))

    # # axs1 = plotnet(5, np.arange(0, len(dtime)), 16)
    # axs2 = plotnet2(
    #     obs.dtime,
    #     obs.swdave,
    #     obs.swuave,
    #     obs.lwdave,
    #     obs.lwuave,
    #     obs.data,
    #     davetime,
    #     7,
    #     np.arange(i1, i2 + 1),
    #     np.array(inds),
    #     12,
    # )

if plot_misc:
    plt.figure()
    plt.subplot(211)
    plt.plot(spack.dtime, spack.lwd, label="lwd, snowpack")
    plt.plot(spack.dtime, spack.lwu, label="lwu, snowpack")
    plt.plot(obs.dtime, obs.data["LWD"], "--", color="blue", label="lwd, obs")
    plt.plot(
        obs.dtime, obs.data["LWU"], "--", color="orange", label="lwu, obs"
    )
    plt.legend()
    plt.ylabel("Longwave radiation (W/m2)")

    plt.figure()
    plt.plot(
        obs.dtime,
        obs.data["LWD"] - obs.data["LWU"],
        "--",
        color="cyan",
        label="lwnet, obs",
    )
    plt.plot(spack.dtime, spack.lwnet, label="lwnet snowpack")
    plt.plot(obs.dtime, np.zeros(np.shape(obs.dtime)), "k:")
    plt.legend()
    plt.ylabel("Longwave radiation (W/m2)")

    plt.figure()
    plt.plot(spack.dtime, spack.swd, label="swd, snowpack")
    plt.plot(spack.dtime, spack.swu, label="swu, snowpack")
    plt.plot(obs.dtime, obs.data["SWD"], "--", color="blue", label="swd, obs")
    plt.plot(
        obs.dtime, obs.data["SWU"], "--", color="orange", label="swu, obs"
    )
    plt.legend()
    plt.ylabel("Shortwave radiation (W/m2)")

    plt.figure()
    plt.plot(spack.dtime, spack.lwnet, label="lwnet snowpack")
    plt.plot(spack.dtime, spack.swnet, label="swnet, snowpack")
    plt.plot(
        obs.dtime,
        obs.data["LWD"] - obs.data["LWU"],
        "--",
        color="cyan",
        label="lwnet, obs",
    )
    plt.plot(
        obs.dtime,
        obs.data["SWD"] - obs.data["SWU"],
        "--",
        color="black",
        label="swnet, obs",
    )
    plt.legend()
    plt.ylabel("Shortwave radiation (W/m2)")
