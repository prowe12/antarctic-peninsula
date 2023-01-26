#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 09:53:04 2022

@author: prowe

"""

import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

# import copy


def dup(val, off=0):
    return [val + off, val + off]


def plotaves(axs, inds):
    """
    For the average over each day, plot:
    2) The net broadband: SWD - SWU, LWD - LWU
    3) The net broadband: SWD - SWU + LWD - LWU
    @param figno  Desired figure number
    @param inds  Indices to selected times
    @param minticks  Minimum number of ticks (days)
    """

    swnet = swdave - swuave
    lwnet = lwdave - lwuave
    net = swdave - swuave + lwdave - lwuave

    # dtime = np.ravel(np.vstack([dtime, dtime]), order='F')
    # swnet = np.ravel(np.vstack([swnet, swnet]), order='F')
    # lwnet = np.ravel(np.vstack([lwnet, lwnet]), order='F')
    # net = np.ravel(np.vstack([net, net]), order='F')

    plt.figure()
    for i in inds:
        day = [davetime[i], davetime[i] + dt.timedelta(days=1)]
        plt.plot(day, dup(swnet[i]), color="orange")
        plt.plot(day, dup(lwnet[i]), "b")

        plt.plot(day, dup(net[i]), "g")


def plotaves2(figno, inds, minticks):
    """
    Plot:
    1) The broadband components: SWD, SWU, LWD, LWU
    2) The net broadband: SWD - SWU, LWD - LWU
    3) The net broadband: SWD - SWU + LWD - LWU
    @param figno  Desired figure number
    @param inds  Indices to selected times
    @param minticks  Minimum number of ticks (days)
    @return The axes
    """

    swnet = swdave - swuave
    lwnet = lwdave - lwuave
    net = swdave - swuave + lwdave - lwuave

    fmt = mdates.DateFormatter("%d")

    fig, axs = plt.subplots(3, 1, figsize=(10, 15), num=figno)
    ht = 0.29
    wid = 0.9
    lft = 0.08

    ax = axs[0]
    ax.set_position([lft, 0.69, wid, ht])
    i = inds[0]
    day = [davetime[i], davetime[i] + dt.timedelta(days=1)]
    ax.plot(day, dup(swdave[i]), color="orange", label="SWD, daily ave")
    ax.plot(day, dup(swuave[i]), "k--", label="SWU, daily ave")
    ax.plot(day, dup(lwdave[i]), "b", label="LWD, daily ave")
    ax.plot(day, dup(lwuave[i]), "c--", label="LWU, daily ave")
    for i in inds:
        day = [davetime[i], davetime[i] + dt.timedelta(days=1)]
        ax.plot(day, dup(swdave[i]), color="orange")
        ax.plot(day, dup(swuave[i]), "k--")
        ax.plot(day, dup(lwdave[i]), "b")
        ax.plot(day, dup(lwuave[i]), "c--")
    ax.set_ylabel("Broadband Radiation (W/m$^2$)")
    ax.margins(x=0)
    ax.set_ylim([0, 500])
    ax.legend()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=20)
    # formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    # ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_formatter(fmt)
    ax.grid()

    ax = axs[1]
    ax.set_position([lft, 0.37, wid, ht])
    i = inds[0]
    day = [davetime[i], davetime[i] + dt.timedelta(days=1)]
    ax.plot(day, dup(swnet[i]), color="orange", label="SWD - SWU, daily ave")
    ax.plot(day, dup(lwnet[i]), "b", label="LWD - LWU, daily ave")
    for i in inds:
        day = [davetime[i], davetime[i] + dt.timedelta(days=1)]
        ax.plot(day, dup(swnet[i]), color="orange")
        ax.plot(day, dup(lwnet[i]), "b")
    ax.set_ylabel("Net Broadband Radiation (W/m$^2$)")
    ax.margins(x=0)
    ax.set_ylim([0, 200])
    ax.grid()
    ax.legend()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=20)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(fmt)
    # ax.set_xticklabels('')

    ax = axs[2]
    ax.set_position([lft, 0.05, 0.9, ht])
    i = inds[0]
    day = [davetime[i], davetime[i] + dt.timedelta(days=1)]
    ax.plot(day, dup(net[i]), "g", label="SWD - SWU + LWD - LWU, daily ave")
    for i in inds:
        day = [davetime[i], davetime[i] + dt.timedelta(days=1)]
        ax.plot(day, dup(net[i]), "g")
    ax.set_ylabel("Net Broadband Radiation (W/m$^2$)")
    ax.margins(x=0)
    ax.set_ylim([0, 200])
    ax.grid()
    ax.legend()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=20)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(fmt)

    return axs


def plotnet(dtime, data, figno, inds, minticks):
    """
    Plot:
    1) The broadband components: SWD, SWU, LWD, LWU
    2) The net broadband: SWD - SWU, LWD - LWU
    3) The net broadband: SWD - SWU + LWD - LWU
    @param figno  Desired figure number
    @param inds  Indices to selected times
    @param minticks  Minimum number of ticks (days)
    @return The axes
    """

    fmt = mdates.DateFormatter("%d")

    # Date range
    dtm = [dtime[i] for i in inds]

    fig, axs = plt.subplots(3, 1, figsize=(10, 15), num=figno)
    ht = 0.29
    wid = 0.9
    lft = 0.08

    ax = axs[0]
    ax.set_position([lft, 0.69, wid, ht])
    ax.plot(dtm, data["SWD"][inds], color="orange", label="SWD")
    ax.plot(dtm, data["SWU"][inds], "k--", label="SWU")
    ax.plot(dtm, data["LWD"][inds], "b", label="LWD")
    ax.plot(dtm, data["LWU"][inds], "c--", label="LWU")
    ax.set_ylabel("Broadband Radiation (W/m$^2$)")
    ax.margins(x=0)
    ax.set_ylim([-10, 500])
    ax.legend()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=20)
    # formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    # ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_formatter(fmt)
    ax.grid()

    ax = axs[1]
    ax.set_position([lft, 0.37, wid, ht])
    ax.plot(
        dtm,
        data["SWD"][inds] - data["SWU"][inds],
        color="orange",
        label="SWD - SWU",
    )
    ax.plot(dtm, data["LWD"][inds] - data["LWU"][inds], "b", label="LWD - LWU")
    ax.plot(dtm, data["SWD"][inds] - data["SWU"][inds], color="orange")
    # ax.plot(dtm, np.zeros(len(dtm)), 'k:')
    ax.set_ylabel("Net Broadband Radiation (W/m$^2$)")
    ax.margins(x=0)
    ax.set_ylim([-100, 150])
    ax.grid()
    ax.legend()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=20)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(fmt)
    # ax.set_xticklabels('')

    ax = axs[2]
    ax.set_position([lft, 0.05, 0.9, ht])
    ax.plot(
        dtm,
        data["SWD"][inds]
        - data["SWU"][inds]
        + data["LWD"][inds]
        - data["LWU"][inds],
        "g",
        label="SWD + LWD - SWU - LWU",
    )
    # ax.plot(dtm, np.zeros(len(dtm)), 'k:')
    ax.set_ylabel("Net Broadband Radiation (W/m$^2$)")
    ax.margins(x=0)
    ax.set_ylim([-100, 150])
    ax.grid()
    ax.legend()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=20)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(fmt)

    return axs


def plotnet2(
    dtime,
    swdave,
    swuave,
    lwdave,
    lwuave,
    data,
    davetime,
    figno,
    inds,
    inds2,
    minticks,
):
    """
    Plot:
    1) The broadband components: SWD, SWU, LWD, LWU
    2) The net broadband: SWD - SWU, LWD - LWU
    3) The net broadband: SWD - SWU + LWD - LWU
    @param figno  Desired figure number
    @param inds  Indices to selected times
    @param minticks  Minimum number of ticks (days)
    @return The axes
    """

    fmt = mdates.DateFormatter("%d")
    swnet = swdave - swuave
    lwnet = lwdave - lwuave
    net = swdave - swuave + lwdave - lwuave

    # Date range
    dtm = [dtime[i] for i in inds]

    fig, axs = plt.subplots(3, 1, figsize=(10, 15), num=figno)
    ht = 0.29
    wid = 0.9
    lft = 0.08
    alph = 0.1

    ax = axs[0]
    ax.set_position([lft, 0.69, wid, ht])
    ax.plot(dtm, data["SWD"][inds], color="orange", alpha=alph)
    ax.plot(dtm, data["SWU"][inds], "k--", alpha=alph)
    ax.plot(dtm, data["LWD"][inds], "b", alpha=alph)
    ax.plot(dtm, data["LWU"][inds], "c", alpha=alph)
    # Plot the averages
    i = inds2[0]
    day = [davetime[i], davetime[i] + dt.timedelta(days=0.999)]
    ax.plot(day, dup(swdave[i]), color="orange", label="SWD, daily ave")
    ax.plot(day, dup(swuave[i]), "k", label="SWU, daily ave")
    ax.plot(day, dup(lwdave[i]), "b", label="LWD, daily ave")
    ax.plot(day, dup(lwuave[i]), "c", label="LWU, daily ave")
    for i in inds2:
        day = [davetime[i], davetime[i] + dt.timedelta(days=0.999)]
        ax.plot(day, dup(swdave[i]), color="orange")
        ax.plot(day, dup(swuave[i]), "k")
        ax.plot(day, dup(lwdave[i]), "b")
        ax.plot(day, dup(lwuave[i]), "c")
    ax.plot(dtm, np.zeros(len(dtm)), "k:")
    ax.set_ylabel("Broadband Radiation (W/m$^2$)")
    ax.margins(x=0)
    ax.set_ylim([-10, 500])
    ax.legend()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=20)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(fmt)

    ax = axs[1]
    ax.set_position([lft, 0.37, wid, ht])
    ax.plot(
        dtm, data["SWD"][inds] - data["SWU"][inds], color="orange", alpha=alph
    )
    ax.plot(dtm, data["LWD"][inds] - data["LWU"][inds], "b", alpha=alph)
    ax.plot(
        dtm, data["SWD"][inds] - data["SWU"][inds], color="orange", alpha=alph
    )
    # Plot the averages
    i = inds2[0]
    day = [davetime[i], davetime[i] + dt.timedelta(days=0.999)]
    ax.plot(day, dup(swnet[i]), color="orange", label="SWD - SWU, daily ave")
    ax.plot(day, dup(lwnet[i]), "b", label="LWD - LWU, daily ave")
    for i in inds2:
        day = [davetime[i], davetime[i] + dt.timedelta(days=0.999)]
        ax.plot(day, dup(swnet[i]), color="orange")
        ax.plot(day, dup(lwnet[i]), "b")
    ax.plot(dtm, np.zeros(len(dtm)), "k:")

    ax.set_ylabel("Net Broadband Radiation (W/m$^2$)")
    ax.margins(x=0)
    ax.set_ylim([-100, 150])
    ax.legend()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=20)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(fmt)
    # ax.set_xticklabels('')

    ax = axs[2]
    ax.set_position([lft, 0.05, 0.9, ht])
    ax.plot(
        dtm,
        data["SWD"][inds]
        - data["SWU"][inds]
        + data["LWD"][inds]
        - data["LWU"][inds],
        "g",
        alpha=alph,
    )
    # Plot the averages
    i = inds2[0]
    day = [davetime[i], davetime[i] + dt.timedelta(days=0.999)]
    ax.plot(day, dup(net[i]), "g", label="SWD - SWU + LWD - LWU, daily ave")
    for i in inds2:
        day = [davetime[i], davetime[i] + dt.timedelta(days=0.999)]
        ax.plot(day, dup(net[i]), "g")
    ax.plot(dtm, np.zeros(len(dtm)), "k:")

    ax.set_ylabel("Net Broadband Radiation (W/m$^2$)")
    ax.margins(x=0)
    ax.set_ylim([-100, 150])
    ax.legend()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=20)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(fmt)
    ax.set_xlabel("Day on March 2022")

    axs[0].text(dtm[400], 450, "a)", fontsize=18)
    axs[1].text(dtm[400], 120, "b)", fontsize=18)
    axs[2].text(dtm[400], 120, "c)", fontsize=18)

    return axs
