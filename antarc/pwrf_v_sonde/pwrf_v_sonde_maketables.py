#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:25 2022

@author: prowe

Compare Polar WRF and radiosonde temperature, RH, and wind speed
at stations in the Antarctic Peninsula: Rothera and Escudero
During Dec. 8 2018 and Feb. 6 2022 foehn events
"""

import datetime as dt
from matplotlib import pyplot as plt

from load_pwrf import LoadPwrf, Pwrf
from load_sonde import IgraSonde, DenialSonde, DenialSondePR


# Base directories
meas_dir = "/Users/prowe/Sync/measurements/"
main_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/"

plt.close("all")


snd_files = [
    "esc_sonde_dd_2018120523.txt",
    "esc_sonde_dd_2018120612.txt",
    "esc_sonde_dd_2018120614.txt",
    "esc_sonde_dd_2018120617.txt",
    "esc_sonde_dd_2018120623.txt",
    "esc_sonde_dd_2018120723.txt",
    "esc_sonde_dd_2018120823.txt",
]
date_strs = [
    "2018/12/05 23 UT",
    "2018/12/06 12 UT",
    "2018/12/06 14 UT",
    "2018/12/06 17 UT",
    "2018/12/06 23 UT",
    "2018/12/07 23 UT",
    "2018/12/08 23 UT",
]


# Upwind: Either Escudero or Rotheraf
cases = []
for snd_file, date_str in zip(snd_files, date_strs):
    cases.append(
        {
            "loc": "Escudero",
            "pwrf_dir": "Dec_2018/pwrf_stations/",
            "pwrf_date": "20181205_08",
            "pwrf_date_fmt": "%Y%m%d_08",
            "snd_dir": meas_dir + "Escudero/radiosondes/datadenial/",
            "snd_file": snd_file,
            "sndfile_fmt": "esc_sonde_dd_%Y%m%d%H.txt",
            "date_str": date_str,
            "snd_class": DenialSondePR,
        },
    )

cases += [
    {
        "loc": "Rothera",
        "pwrf_dir": "Dec_2018/pwrf_stations/",
        "pwrf_date": "20181205_08",
        "pwrf_date_fmt": "%Y%m%d_08",
        "snd_dir": meas_dir + "Rothera/radiosonde/",
        "snd_file": "Rothera_2018120511.dat",
        "sndfile_fmt": "Rothera_%Y%m%d%H.dat",
        "date_str": "2018/12/05 12 UT",
        "snd_class": DenialSonde,
    },
    {
        "loc": "Rothera",
        "pwrf_dir": "Dec_2018/pwrf_stations/",
        "pwrf_date": "20181205_08",
        "pwrf_date_fmt": "%Y%m%d_08",
        "snd_dir": meas_dir + "Rothera/radiosonde/",
        "snd_file": "Rothera_2018120611.dat",
        "sndfile_fmt": "Rothera_%Y%m%d%H.dat",
        "date_str": "2018/12/06 12 UT",
        "snd_class": DenialSonde,
    },
    {
        "loc": "Rothera",
        "pwrf_dir": "Dec_2018/pwrf_stations/",
        "pwrf_date": "20181205_08",
        "pwrf_date_fmt": "%Y%m%d_08",
        "snd_dir": meas_dir + "Rothera/igra/derived/",
        "snd_file": "Rothera_igra_20181207_12.txt",
        "sndfile_fmt": "Rothera_igra_%Y%m%d_%H.txt",
        "date_str": "2018/12/07 12 UT",
        "snd_class": IgraSonde,
    },
    {
        "loc": "Marambio",
        "pwrf_dir": "Dec_2018/pwrf_stations/",
        "pwrf_date": "20181205_08",
        "pwrf_date_fmt": "%Y%m%d_08",
        "snd_dir": meas_dir + "Marambio/radiosonde/",
        "snd_file": "Marambio_2018120511.dat",
        "date_str": "2018/12/05 12 UT",
        "sndfile_fmt": "Marambio_%Y%m%d%H.dat",
        "snd_class": DenialSonde,
    },
]


dhour = 1.0

for idate in range(len(cases)):
    # idate = 0

    case = cases[idate]
    pwrfdir = main_dir + case["pwrf_dir"]  # pwrfdirs[idate]
    pwrf_date_fmt = case["pwrf_date_fmt"]  # pwrf_date_fmts[idate]
    pwrf_date = case["pwrf_date"]  # pwrf_dates[idate]
    loc = case["loc"]

    # Files and dir
    snd_dir = case["snd_dir"]
    snd_file = case["snd_file"]
    sndfile_fmt = case["sndfile_fmt"]
    date_str = case["date_str"]
    upwndSndClass = case["snd_class"]

    # Create the PWRF paths for the case location
    starttime = dt.datetime.strptime(case["pwrf_date"], case["pwrf_date_fmt"])
    tc_file = pwrfdir + "PWRF_D03_" + loc + "_tc_" + pwrf_date + "_hly.nc"
    rh_file = pwrfdir + "PWRF_D03_" + loc + "_rh_" + pwrf_date + "_hly.nc"
    wnd_file = pwrfdir + "PWRF_D03_" + loc + "_uvmet_" + pwrf_date + "_hly.nc"

    # The radiosonde data
    sonde = upwndSndClass(snd_dir, snd_file, sndfile_fmt)

    # The PWRF data for Rothera, including the index to the PWRF time of interest
    allPwrf = LoadPwrf(tc_file, rh_file, wnd_file, starttime, dhour)
    pwrf = Pwrf(allPwrf, sonde.dtime)
    pwrf.trim(sonde.lev)

    # Get difference statistics: max, min, mean, rms difference
    dtemp, drh, dwind = sonde.get_diffs(pwrf)

    # Uncomment one of the following blocks to print the table

    # RH
    # maxy = max(abs(drh["max"]), abs(drh["min"]))
    # print(f"\n{case['loc'][:7]}  {date_str} ", end="")
    # print(f"{round(maxy,1)} ", end="")
    # print(f"{round(drh['mean'],1)} ", end="")
    # print(f"{round(drh['rms'],1)} ", end="")

    # Temperature and wind speed
    maxt = max(abs(dtemp["max"]), abs(dtemp["min"]))
    maxws = max(abs(dwind["max"]), abs(dwind["min"]))
    print(f"\n{case['loc'][:7]}  {date_str} ", end="")
    print(f"{round(maxt,1)} ", end="")
    print(f"{round(dtemp['mean'],1)} ", end="")
    print(f"{round(dtemp['rms'],1)}     ", end="")
    print(f"{round(maxws,1)} ", end="")
    print(f"{round(dwind['mean'],1)} ", end="")
    print(f"{round(dwind['rms'],1)} ", end="")
