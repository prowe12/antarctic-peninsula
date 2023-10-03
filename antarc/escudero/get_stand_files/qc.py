#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 14:30:22 2023

@author: prowe

Vaisala installed in 2020
2020: A lot of gaps
2021: No personnel in winter - need to retrieve memory card
2022: have data
2023: Transferring data to lidar computer - csv file every minute
Edgardo code reads those files and creates json files
temperature and humidity sensor were reading incorrect values
Fernanda will determine dates 2020/2021 old sensor; 2022 new sensor
humidity was really bad, always supersaturated
temperature may have been incorrect too

Army station next to Escudero has weather station
Before launch were comparing our temperature to Army station values

For gap from 2023/04 - 2023/08 use 15 minute data
"""


import os
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import bisect


from antarc.hypsometric import hypsometric_for_z, hypsometric


def get_files(direc: str, sampfname: str, ipfx=3) -> list[str]:
    """
    Get files from directory with names similar to sampfname
    @param direc  The directory
    @param sampfname  The sample file names
    @param ipfx  Number of starting characters to match
    @returns  List of file names
    """
    allfiles = os.listdir(direc)
    fnames = []
    for fname in allfiles:
        if len(fname) == len(sampfname) and fname[:5] == sampfname[:5]:
            fnames.append(fname)
    fnames.sort()
    return fnames


def get_sonde_data():
    """
    Get radiosonde data
    """
    # Radiosonde Data
    sonde["files"] = get_files(
        sonde["direc"], sonde["sampfname"], sonde["pfx"]
    )
    dtime = []
    alt0 = []
    alt1 = []
    alt2 = []
    alt3 = []
    temp0 = []
    temp1 = []
    temp2 = []
    temp3 = []
    press0 = []
    press1 = []
    press2 = []
    press3 = []
    wdir0 = []
    wdir1 = []
    wspd0 = []
    wspd1 = []
    for fname in sonde["files"]:
        df = pd.read_csv(sonde["direc"] + fname, na_values="nan")
        df.columns = df.columns.str.replace(" ", "")
        for col in df.columns:
            if col != "Date":
                df[col] = df[col].astype(float)
        temp0.append(df["Temp(C)"][0])
        temp1.append(df["Temp(C)"][1])
        temp2.append(df["Temp(C)"][2])
        temp3.append(df["Temp(C)"][3])
        press0.append(df["P(hPa)"][0])
        press1.append(df["P(hPa)"][1])
        press2.append(df["P(hPa)"][2])
        press3.append(df["P(hPa)"][3])
        alt0.append(df["Alt(m)"][0])
        alt1.append(df["Alt(m)"][1])
        alt2.append(df["Alt(m)"][2])
        alt3.append(df["Alt(m)"][3])
        wdir0.append(df["wdir(deg)"][0])
        wdir1.append(df["wdir(deg)"][1])
        wspd0.append(df["wspd(kn)"][0])
        wspd1.append(df["wspd(kn)"][1])
        dtime0 = dt.datetime.strptime(df["Date"][0], "%Y-%m-%d_%H:%M:%S")
        dtime.append(dtime0)
    snd_dtime = np.array(dtime)
    snd_alt = np.vstack(
        [np.array(alt0), np.array(alt1), np.array(alt2), np.array(alt3)]
    )
    snd_temp = np.vstack(
        [np.array(temp0), np.array(temp1), np.array(temp2), np.array(temp3)]
    )
    snd_press = np.vstack(
        [
            np.array(press0),
            np.array(press1),
            np.array(press2),
            np.array(press3),
        ]
    )
    snd_wdir = np.vstack([np.array(wdir0), np.array(wdir1)])
    snd_wspd = np.vstack([np.array(wspd0), np.array(wspd1)])

    return snd_dtime, snd_alt, snd_temp, snd_press, snd_wdir, snd_wspd


def get_frei_201701():
    """
    Get Frei data from 2017/01
    Notes:
    Frei met data: 2017/01, every hour
    QFE is Station pressure
    QFF is Mean Sea level pressure
    """
    met = frei_201701
    met["files"] = get_files(met["direc"], met["sampfname"], met["pfx"])
    dtime = []
    temp = []
    press = []
    slp = []
    wdir = []
    wspd = []
    rhw = []
    # TODO: for now, just get the 201701 file. Eventually we will use the
    # merged files
    for fname in met["files"][:1]:
        dfm = pd.read_csv(met["direc"] + fname)
        temp += list(dfm["temperature (C)"])
        press += list(dfm["QFE (hPa)"])
        slp += list(dfm["QFF (hPa)"])
        wdir += list(dfm["Wdir (deg)"])
        wspd += list(dfm["Wspd (knot)"])
        rhw += list(dfm["RH (%)"])
        date0 = list(dfm["Date"])
        dtime0 = [dt.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in date0]
        dtime += dtime0
    met_dtime = np.array(dtime)
    met_temp = np.array(temp)
    met_press = np.array(press)
    met_slp = np.array(slp)
    met_wdir = np.array(wdir)
    met_wspd = np.array(wspd)
    met_rh = np.array(rhw)

    return met_dtime, met_temp, met_press, met_slp, met_wdir, met_wspd, met_rh


def get_frei_15():
    """
    Frei met data: every 15 minutes
    """
    frei["files"] = get_files(frei["direc"], frei["sampfname"], frei["pfx"])
    dtime = []
    temp = []
    press = []
    slp = []
    wdir = []
    wspd = []
    rhw = []
    for fname in frei["files"]:
        dfm = pd.read_csv(frei["direc"] + fname)
        temp += list(dfm["temperature (C)"])
        press += list(dfm["Pressure (hPa)"])
        slp += list(dfm["Sea level pressure (hPa)"])
        wdir += list(dfm["Wind direction (degree)"])
        wspd += list(dfm["Wind speed (knot)"])
        rhw += list(dfm["RH (%)"])
        date0 = list(dfm["Date"])
        dtime0 = [dt.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in date0]
        dtime += dtime0

    return (
        np.array(dtime),
        np.array(temp),
        np.array(press),
        np.array(slp),
        np.array(wdir),
        np.array(wspd),
        np.array(rhw),
    )


def get_tarp_met():
    """
    TARP AWS data
    """
    aws["files"] = get_files(aws["direc"], aws["sampfname"], aws["pfx"])
    dtime = []
    temp = []
    press = []
    wdir = []
    wspd = []
    rhw = []

    # frei["files"] = get_files(frei["direc"], frei["sampfname"], frei["pfx"])
    for fname in aws["files"]:
        dfm = pd.read_csv(aws["direc"] + fname)
        temp += list(dfm["Temperature(C)"])
        press += list(dfm["Pressure(hPa)"])
        wdir += list(dfm["Wdir(deg)"])
        wspd += list(dfm["Wspd"])
        rhw += list(dfm["RH(%)"])
        date0 = list(dfm["Date"])
        dtime0 = [dt.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in date0]
        dtime += dtime0

    return (
        np.array(dtime),
        np.array(temp),
        np.array(press),
        np.array(wdir),
        np.array(wspd),
        np.array(rhw),
    )


def get_frei_missing():
    """Determine what Frei met data is missing by interpolating
    it to sounding times"""

    missing_dtime = []
    imet = []
    for dtime in snd_dtime:
        ind = bisect.bisect_right(met_dtime, dtime)
        if ind == 0:
            missing_dtime.append(dtime)
            imet.append(0)
            continue
        elif ind >= len(met_dtime) - 1:
            missing_dtime.append(dtime)
            raise ValueError("Not expecting this")

        # Time difference for nearest point
        dtbef = (dtime - met_dtime[ind - 1]).total_seconds() / 60
        dtaft = (met_dtime[ind] - dtime).total_seconds() / 60
        if dtbef < 0 or dtaft < 0:
            raise ValueError("Not expecting this")
        # Check if we have a sounding within 35 minutes
        if dtbef > 35 and dtaft > 35:
            missing_dtime.append(dtime)

        if dtbef < dtaft:
            imet.append(ind - 1)
        else:
            imet.append(ind)
    imet = np.array(imet)

    # Plot time series of all temperature and pressure data
    plt.figure()
    plt.subplot(211)
    plt.plot(met_dtime, met_temp, color="gray", label="Frei")
    plt.plot(met_dtime[imet], met_temp[imet], "s", label="Frei")
    plt.plot(aws_dtime, aws_temp, label="TARP AWS")
    plt.plot(snd_dtime, snd_temp[0, :], "r.", label="sonde")
    plt.legend()
    plt.ylabel("Temperature (C)")
    # plt.xlim([snd_dtime[0] - dday, snd_dtime[-1] + dday])

    plt.subplot(212)
    plt.plot(met_dtime, met_slp, "k", label="Frei SLP")
    plt.plot(met_dtime, met_press, color="gray", label="Frei Station")
    plt.plot(met_dtime[imet], met_press[imet], "s", label="Frei Station")
    plt.plot(aws_dtime, aws_press, label="TARP AWS")
    plt.plot(snd_dtime, snd_press[0, :], "r.", label="sonde1")
    plt.plot(snd_dtime, snd_press[1, :], ".", color="pink", label="sonde2")
    plt.plot(snd_dtime, snd_press[2, :], ".", color="pink", label="sonde2")
    plt.plot(snd_dtime, snd_press[3, :], "c.", label="sonde2")
    plt.ylabel("Pressure (hPa)")

    ddtime = met_dtime[imet] - snd_dtime
    dmin = np.array([x.total_seconds() / 60 for x in ddtime])

    plt.figure()
    plt.semilogy(snd_dtime, abs(dmin), ".")

    plt.figure()
    plt.plot(met_temp[imet], snd_temp[0, :], ".")
    plt.plot([-10, 7], [-10, 7])


# Root measurement directory
meas_dir = "/Users/prowe/Sync/measurements/"
# TARP AWS data
aws = {
    "direc": meas_dir + "Escudero/aws/std/",
    "sampfname": "tarp02_met_202304.txt",
    "pfx": 10,
}
# Frei met data: 2017/01
frei_201701 = {
    "direc": meas_dir + "Escudero/frei_met/std/",
    "sampfname": "esc_met_201701.txt",
    "pfx": 7,
}
# Frei met data: post 2017/01
frei = {
    "direc": meas_dir + "Escudero/frei_met/by_month_every_15min/",
    "sampfname": "esc_met_201803.txt",
    "pfx": 7,
}
# Radiosondes
sonde = {
    "direc": meas_dir + "Escudero/radiosondes/graw/",
    "sampfname": "esc_sonde_v0_2023081018.txt",
    "pfx": 12,
}


snd_dtime, snd_alt, snd_temp, snd_press, snd_wdir, snd_wspd = get_sonde_data()

# Get and combine Frei met data in different formats from 2017/01 and after
dtime1, temp1, press1, slp1, wdir1, wspd1, rhw1 = get_frei_201701()
dtime, temp, press, slp, wdir, wspd, rhw = get_frei_15()

met_dtime = np.hstack([dtime1, dtime])
met_temp = np.hstack([temp1, temp])
met_press = np.hstack([press1, press])
met_slp = np.hstack([slp1, slp])
met_wdir = np.hstack([wdir1, wdir])
met_wspd = np.hstack([wspd1, wspd])
met_rh = np.hstack([rhw1, rhw])

# Get TARP weather station data
aws_dtime, aws_temp, aws_press, aws_wdir, aws_wspd, aws_rh = get_tarp_met()


# Date ranges for field campaigns
dt_ranges = [
    dt.datetime(2017, 1, 1),
    dt.datetime(2017, 7, 1),
    dt.datetime(2018, 7, 1),
    dt.datetime(2019, 7, 1),
    dt.datetime(2021, 7, 1),
    dt.datetime(2023, 1, 1),
    dt.datetime(2023, 9, 1),
]

# Create figures
dday = dt.timedelta(5)
ifig = 0
for i, dt_lo in enumerate(dt_ranges[:-1]):
    ifig += 1
    # Temperature
    plt.figure(ifig)
    plt.clf()
    # met data
    inds = np.where(met_dtime >= dt_lo)[0]
    inds = np.intersect1d(inds, np.where(met_dtime < dt_ranges[i + 1])[0])
    plt.subplot(211)
    plt.plot(met_dtime[inds], met_temp[inds], label="Frei")
    plt.subplot(212)
    plt.plot(met_dtime[inds], met_press[inds], label="Frei")
    # TARP aws data
    inds = np.where(aws_dtime >= dt_lo)[0]
    inds = np.intersect1d(inds, np.where(aws_dtime < dt_ranges[i + 1])[0])
    if len(inds) > 0:
        plt.subplot(211)
        plt.plot(aws_dtime[inds], aws_temp[inds], label="TARP AWS")
        plt.subplot(212)
        plt.plot(aws_dtime[inds], aws_press[inds] - 3, label="TARP AWS:P-3")
    # sounding data
    inds = np.where(snd_dtime >= dt_lo)[0]
    inds = np.intersect1d(inds, np.where(snd_dtime < dt_ranges[i + 1])[0])
    plt.subplot(211)
    plt.plot(snd_dtime[inds], snd_temp[0, inds], "r.", label="sonde")
    # plt.plot(snd_dtime[inds], snd_temp[3, inds], "k.", label="sonde 4")
    plt.legend()
    plt.ylabel("Temperature (C)")
    plt.xlim([snd_dtime[inds[0]] - dday, snd_dtime[inds[-1]] + dday])
    plt.subplot(212)
    plt.plot(snd_dtime[inds], snd_press[0, inds], "r.", label="sonde")
    # plt.plot(snd_dtime[inds], snd_press[3, inds], "k.", label="sonde 4")
    plt.ylabel("Pressure (hPa)")
    plt.xlim([snd_dtime[inds[0]] - dday, snd_dtime[inds[-1]] + dday])
    plt.legend()


# Notes regarding TARP02 Vaisala weather station
# There is no data for either json_daily_backup or 15 min backup between:
# 2023-01-18 20:32:07,986.5,3.8,89.0,295.0,17.0, and
# 2023-01-21 14:34:07,994.7,3.9,92.0,321.0,15.0,
#
# 2023-03-26 07:47:07,971.3,1.2,95.0,322.0,10.0,  and
# 2023-04-08 14:56:07,982.6,1.0,84.0,284.0,13.0,

# There is no data for json_daily_backup *but there is* for the 15 min backup
# between the following (and thus 15 min backup is pasted into "std" files)
# 2023-05-25 12:26:07,983.3,1.3,89.0,60.0,17.0, and
# 2023-05-31 23:59:07,981.9,0.8,101.0,294.0,16.0,
#
# 2023-06-19 12:15:07,973.6,-0.8,92.0,nan,nan, and
# 2023-06-21 11:55:07,997.6,-5.1,93.0,120.0,12.0,
#
# 2023-08-01 00:00:09,978.9,-1.0,95.0,321.0,10.0, and
# 2023-08-18 03:00:07,1013.3,-1.0,96.0,nan,nan,
#
# 2023-09-12 23:56:07,1007.2,-6.5,72.0,62.0,3.0, and
# 2023-09-15 22:44:07,990.0,-10.7,85.0,119.0,27.0,

# Bad soundings determined visually:
# 2017/01/19 11: T

# 2017/12/08 18: T

# 2018/12: ?

# 2019/09/13 14: T
# 2019/09/15 14: T
# 2019/09/19 14: T
# 2019/09/23 14: T
# 2019/09/24 15: T
# 2019/09/27 03: T

# 2022/05/07 12: T
# 2022/05/10 12: T

# Determine what frei data is missing
get_frei_missing()


# Find TARP AWS data near in time to soundings
# For this date range, compare TARP, Frei, and soundings
# Determine the effective Frei station pressure and alt relative to soundings
# and interpolate met data (temp, rh, winds) to height of Frei.
# Compare these values

# Find TARP AWS data near in time to soundings
imet = []
itarp = []
isnd = []
for i, dtime in enumerate(snd_dtime):
    # Skip any sounding times before or after TARP AWS was in operation
    ind = bisect.bisect_right(aws_dtime, dtime)
    if ind == 0:
        continue
    elif ind >= len(aws_dtime) - 1:
        continue

    # Time difference for nearest point
    dtbef = (dtime - aws_dtime[ind - 1]).total_seconds() / 60
    dtaft = (aws_dtime[ind] - dtime).total_seconds() / 60
    if dtbef < 0 or dtaft < 0:
        raise ValueError("Not expecting this")
    # Check if we have a sounding within 35 minutes
    if dtbef > 35 and dtaft > 35:
        print(
            f"Using dummy value for {dtime}: time diff from TARP met is "
            + f"{round(min(dtaft,dtbef))} minutes"
        )
        isnd.append(i)
        itarp.append(0)
    else:
        if dtbef < dtaft:
            itarp.append(ind - 1)
        else:
            itarp.append(ind)
        isnd.append(i)

    # As above but for Frei met data
    ind = bisect.bisect_right(met_dtime, dtime)
    if ind == 0 or ind >= len(aws_dtime) - 1:
        raise ValueError("not expecting this")
    dtbef = (dtime - met_dtime[ind - 1]).total_seconds() / 60
    dtaft = (met_dtime[ind] - dtime).total_seconds() / 60
    if dtbef < 0 or dtaft < 0:
        raise ValueError("Not expecting this")
    # Check if we have a sounding within 35 minutes
    if dtbef > 35 and dtaft > 35:
        print(
            f"Using dummy value for {dtime}: time diff from Frei met is "
            + f"{round(min(dtaft,dtbef))} min"
        )
        imet.append(0)
    elif dtbef < dtaft:
        imet.append(ind - 1)
    else:
        imet.append(ind)

isnd = np.array(isnd)
itarp = np.array(itarp)
imet = np.array(imet)

# Remove dummy values
met_temp[0] = np.nan
met_press[0] = np.nan
aws_temp[0] = np.nan
aws_press[0] = np.nan


# Scatter plots comparing temepratures
# excluding points where TARP AWS or Frei are missing
plt.figure()
plt.subplot(221)
plt.plot(met_temp[imet], snd_temp[0, isnd], ".")
plt.plot([-10, 7], [-10, 7])
plt.xlabel("T at Frei (C)")
plt.ylabel("T from sounding at surface (C)")
plt.subplot(222)
plt.plot(aws_temp[itarp], snd_temp[0, isnd], ".")
plt.plot([-10, 7], [-10, 7])
plt.xlabel("T from Tarp weather station (C)")
plt.ylabel("T from sounding at surface (C)")
plt.subplot(223)
plt.plot(met_temp[imet], aws_temp[itarp], ".")
plt.plot([-10, 7], [-10, 7])
plt.xlabel("T at Frei (C)")
plt.ylabel("T from Tarp weather station (C)")


# Time series of pressures in 2022
# plt.figure()
# plt.plot(snd_dtime[isnd], aws_press[itarp], "s", label="TARP")
# plt.plot(snd_dtime[isnd], met_press[imet], "o", label="Frei")
# plt.plot(snd_dtime[isnd], snd_press[0, isnd], ".", label="Sonde")
# plt.legend()

# Make sure all TARP AWS points are visible
aws_temp[0] = 20.0
aws_press[0] = 1020.0

# Time series of pressure differences
plt.figure()
plt.plot(
    snd_dtime[isnd], aws_press[itarp] - met_press[imet], "o", label="TARP-Frei"
)
plt.plot(
    snd_dtime[isnd],
    aws_press[itarp] - snd_press[0, isnd],
    ".",
    label="TARP-sonde",
)
plt.plot([snd_dtime[isnd[0]], snd_dtime[isnd[-1]]], [0, 0], "k:")
plt.legend()
plt.ylabel("Pressure difference (hPa)")

# Time series of temperature differences
plt.figure()
plt.plot(
    snd_dtime[isnd], aws_temp[itarp] - met_temp[imet], "o", label="TARP-Frei"
)
plt.plot(
    snd_dtime[isnd],
    aws_temp[itarp] - snd_temp[0, isnd],
    ".",
    label="TARP-sonde",
)
plt.plot([snd_dtime[isnd[0]], snd_dtime[isnd[-1]]], [0, 0], "k:")
plt.legend()
plt.ylabel("Temperature difference (C)")

# As shown from the following:
# 1) The pressures at TARP and the Frei met station are not very well
#    related by the hypsometric equation, which is only within 10 m.
# 2) Temperatures above TARP are higher at the ~altitude of the Frei station
#
# Caveats:
# 1) The Frei pressure is probably for the airfield not the station
# 2) We don't know the Frei station or airfield altitudes
# 3) The accuracy of the launch altitude of the sounding is not known
#
# Based on this analysis, it seems reasonable to keep all soundings
# with values that are within some threshold for each variable,
# and redo those outside the threshold.
# We can first make sure the surface values are as expected.

# Use the hypsometric equation to determine the relative altitudes
htdiff = np.zeros(len(imet))
snd_temp_met_alt = np.zeros(len(imet))
snd_temp_met_press = np.zeros(len(imet))
qff = np.zeros(len(imet))
for i in range(len(imet)):
    temp = np.array([aws_temp[itarp[i]], met_temp[imet[i]]]) + 273.15
    press = np.array([aws_press[itarp[i]], met_press[imet[i]]])
    rhw = np.array([aws_rh[itarp[i]], met_rh[imet[i]]])

    alts = hypsometric_for_z(press, temp, rhw, 30.0)  # m
    htdiff[i] = np.diff(alts)[0]  # m

    # Get the temperature at the frei height and pressure
    ht_frei = 44.0  # snd_alt[0, isnd[i]] + htdiff[i]
    snd_temp_met_alt[i] = np.interp(
        ht_frei, snd_alt[:, isnd[i]], snd_temp[:, isnd[i]]
    )
    snd_temp_met_press[i] = np.interp(
        -met_press[imet[i]], -snd_press[:, isnd[i]], snd_temp[:, isnd[i]]
    )

    # Calculate QFF
    temp = np.array([met_temp[imet[i]], met_temp[imet[i]]]) + 273.15
    rhw = np.array([met_rh[imet[i]], met_rh[imet[i]]])
    alts = np.array([44.0, 0.0])
    qff_vec = hypsometric(alts, temp, rhw, met_press[imet[i]])
    qff[i] = qff_vec[1]

# Plot qff calculated here vs by Frei met staff
plt.figure()
plt.plot(met_dtime[imet], met_slp[imet] - qff, ".")
plt.ylabel("My QFF - Frei QFF")
# plt.plot(met_dtime[imet], met_slp[imet]-met_press[imet],'.')
# plt.plot([960, 1020],[960,1020])

# Plot the differences in height calculated from the pressure differences
plt.figure()
plt.plot(snd_dtime[isnd], 20.0 + htdiff, ".")
plt.ylabel("Frei Height determined from pressure (m)")


plt.figure()
plt.plot(snd_dtime[isnd], snd_temp[:, isnd].T, ".", color="gray")
plt.plot(snd_dtime[isnd], snd_temp_met_alt, "o", label="interp for z")
plt.plot(snd_dtime[isnd], snd_temp_met_press, "s", label="interp for P")
plt.plot(snd_dtime[isnd], met_temp[imet], "+", label="Frei")
plt.legend()

plt.figure()
plt.plot(snd_dtime[isnd], snd_press[:, isnd].T, ".", color="gray")
plt.plot(snd_dtime[isnd], aws_press[itarp], "s", label="TARP")
plt.plot(snd_dtime[isnd], met_press[imet], "+", label="Frei")

# e = RH / 100 * esw(temp)  # I assume wrt water, but should be?
# w = eps * (e / (P - e))
# tvirt = T * ((1 + (w / eps)) / (1 + w))
# TvBar = (tvirt[i] + tvirt[i + 1]) / 2
# dz = -np.log(met_press[imet] / aws_press[itarp]) * gas_const * TvBar / grav
# hypsometric_for_z(P, T, RH, zstart)
