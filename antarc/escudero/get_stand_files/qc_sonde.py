#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 09:36:52 2023

@author: prowe
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


def plot_sounding():
    # Plot it
    plt.figure()
    plt.subplot(221)
    plt.plot(dfm["P(hPa)"], dfm["Alt(m)"])
    plt.plot(dfm["P(hPa)"][:iburst], dfm["Alt(m)"][:iburst], "s")
    plt.plot(dfm["P(hPa)"][:isnd], dfm["Alt(m)"][:isnd], ".")
    plt.ylabel("Pressure")
    plt.subplot(222)
    plt.plot(dfm["Temp(C)"][:iburst], dfm["Alt(m)"][:iburst])
    plt.ylabel("Temperature")
    plt.subplot(223)
    plt.plot(dfm["RH(%)"][:iburst], dfm["Alt(m)"][:iburst])
    plt.ylabel("RH(%)")
    plt.subplot(224)
    plt.plot(dfm["wspd(kn)"][:iburst], dfm["Alt(m)"][:iburst])
    plt.ylabel("wspd(kn)")
    return


# Root measurement directory
meas_dir = "/Users/prowe/Sync/measurements/"


# Radiosondes
sonde = {
    "direc": meas_dir + "Escudero/radiosondes/graw/",
    "sampfname": "esc_sonde_v0_2023081018.txt",
    "pfx": 12,
}
sonde["files"] = get_files(sonde["direc"], sonde["sampfname"], sonde["pfx"])

# Load in the radiosondes
for fname in sonde["files"]:
    # fname = sonde["files"][0]
    dfm = pd.read_csv(sonde["direc"] + fname, na_values="nan")

    dfm.columns = dfm.columns.str.replace(" ", "")
    for col in dfm.columns:
        if col != "Date":
            dfm[col] = dfm[col].astype(float)

    alt = np.array(dfm["Alt(m)"][:])
    alt[alt > 1000000] = -999

    # # burst height
    # iburst = np.argmax(alt)
    # alt[alt<0] == np.nan

    # # 50 m below burst height
    # isnd = np.where(alt >= alt[iburst]-50)[0][0]

    # plot_sounding()

    # Quality control
    # 1) Check for surface wind speed of zero
    if dfm["wspd(kn)"][0] <= 0.1:
        print(f"{fname}: Wind speed is zero")
    # 2) Check if near-surface pressures increse
    if np.any(np.diff(dfm["P(hPa)"][:10]) >= 0):
        print(f"{fname}: Pressure change >= 0 within first 10 points")
    # 3) Check if aloft pressures decrease
    # elif np.any(np.diff(dfm["P(hPa)"][:100]) >= 0):
    #     print(f"{fname}: Pressure change >=0 win 100 pts {dfm['P(hPa)'][100]}")
    # 5) And higher
    # elif np.any(np.diff(dfm["P(hPa)"][:2000]) >= 0):
    #     print(f"{fname}: Pressure change >=0 win 1000 pts {dfm['P(hPa)'][1000]}")
