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
import numpy as np

# Parameter modules
from antarc import params
from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import atm_prof_params as esc_atm_prof_params
from antarc.rothera.parameters import roth_params, roth_atm_prof_params
from antarc.marambio.parameters import mbio_atm_prof_params, mbio_params

from antarc.pwrf_v_sonde.load_pwrf import LoadPwrf, Pwrf
from antarc.pwrf_v_sonde.load_sonde import IgraSonde, DenialSonde
from antarc.pwrf_v_sonde.load_sonde import DenialSondePR

from antarc.pwrf_v_sonde.get_era_profs import get_era


# def sigf(val, digits):
#     if digits != 2:
#         raise ValueError("Digits must be 2!")
#     # Round original value to 2 sig figs
#     if val >= 10:
#         return f"{round(val, digits-2):0.0f}"
#     if val >= 1:
#         return f"{round(val, digits-1):0.1f}"
#     return f"{round(val, digits):.2f}"


# Select which table(s) to print out
# options: pwrf or era5
# rh or temp
table_model = "era5"
table_type = "rh"

print(table_model)
print(table_type)

# Base directories
meas_dir = params.MEAS_DIR
main_dir = params.PROJECT_DIR + "case_studies/"


snd_files = [
    "esc_sonde_dd_v0_2018120523.txt",
    "esc_sonde_dd_v0_2018120612.txt",
    "esc_sonde_dd_v0_2018120614.txt",
    "esc_sonde_dd_v0_2018120617.txt",
    "esc_sonde_dd_v0_2018120623.txt",
    "esc_sonde_dd_v0_2018120723.txt",
    "esc_sonde_dd_v0_2018120823.txt",
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


era_pfx = {
    "Escudero": "era5_esc_",
    "Rothera": "era5_roth_",
    "Marambio": "era5_mbio_",
}
# Location
station_params = {
    "Escudero": [esc_params, esc_atm_prof_params],
    "Rothera": [roth_params, roth_atm_prof_params],
    "Marambio": [mbio_params, mbio_atm_prof_params],
}
erafmt = "%Y/%m/%d %H UT"

# Upwind: Either Escudero or Rothera
cases = []
for snd_file, date_str in zip(snd_files, date_strs):
    cases.append(
        {
            "loc": "Escudero",
            "pwrf_dir": "Dec_2018/pwrf_stations/",
            "pwrf_date": "20181205_08",
            "pwrf_date_fmt": "%Y%m%d_08",
            "era_dir": meas_dir + "Escudero/era5/",
            "snd_dir": meas_dir + "Escudero/radiosondes/datadenial/v0/",
            "snd_file": snd_file,
            "sndfile_fmt": "esc_sonde_dd_v0_%Y%m%d%H.txt",
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
        "era_dir": meas_dir + "Rothera/era5/",
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
        "era_dir": meas_dir + "Rothera/era5/",
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
        "era_dir": meas_dir + "Rothera/era5/",
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
        "era_dir": meas_dir + "Rothera/era5/",
        "snd_dir": meas_dir + "Marambio/radiosonde/",
        "snd_file": "Marambio_2018120511.dat",
        "date_str": "2018/12/05 12 UT",
        "sndfile_fmt": "Marambio_%Y%m%d%H.dat",
        "snd_class": DenialSonde,
    },
]


dhour = 1.0


sum_max_dt = 0
sum_mean_dt = 0
sum_rms_dt = 0
sum_max_dws = 0
sum_mean_dws = 0
sum_rms_dws = 0
sum_max = 0
sum_mean = 0
sum_rms = 0
npts = 0

for idate in range(len(cases)):
    # idate = 0

    case = cases[idate]

    # Sonde files and dir
    snd_dir = case["snd_dir"]
    snd_file = case["snd_file"]
    sndfile_fmt = case["sndfile_fmt"]
    date_str = case["date_str"]
    upwndSndClass = case["snd_class"]

    # The radiosonde data
    sonde = upwndSndClass(snd_dir, snd_file, sndfile_fmt)

    if table_model == "pwrf":
        pwrfdir = main_dir + case["pwrf_dir"]
        pwrf_date_fmt = case["pwrf_date_fmt"]
        pwrf_date = case["pwrf_date"]
        loc = case["loc"]

        # Create the PWRF paths for the case location
        starttime = dt.datetime.strptime(pwrf_date, pwrf_date_fmt)
        tc_file = pwrfdir + "PWRF_D03_" + loc + "_tc_" + pwrf_date + "_hly.nc"
        rh_file = pwrfdir + "PWRF_D03_" + loc + "_rh_" + pwrf_date + "_hly.nc"
        wnd_file = (
            pwrfdir + "PWRF_D03_" + loc + "_uvmet_" + pwrf_date + "_hly.nc"
        )

        # PWRF data, including index to the PWRF time of interest
        allPwrf = LoadPwrf(tc_file, rh_file, wnd_file, starttime, dhour)
        pwrf = Pwrf(allPwrf, sonde.dtime)
        # pwrf.trim(sonde.lev)

        # Get difference statistics: max, min, mean, rms difference
        dtemp, drh, dwind = sonde.get_diffs(
            pwrf.level, pwrf.temp, pwrf.rh, pwrf.wind
        )

    elif table_model == "era5":

        direc = main_dir + case["era_dir"]
        # date_fmt = case["pwrf_date_fmt"]
        loc = case["loc"]
        pfx = era_pfx[loc]

        sta_params, prof_params = station_params[loc]

        era = get_era(pfx, sta_params, prof_params, case["date_str"], erafmt)

        # Get difference statistics: max, min, mean, rms difference
        wspd = np.sqrt(era.uwind**2 + era.vwind**2)
        temp = era.T - 273.15
        dtemp, drh, dwind = sonde.get_diffs(era.P, temp, era.rh, wspd)

    if table_type == "rh":
        maxy = max(abs(drh["max"]), abs(drh["min"]))
        print(f"\n{case['loc']}, {date_str[:10]},{date_str[10:13]},  ", end="")
        print(f"{round(maxy,1)}, ", end="")
        print(f"{round(drh['mean'],1)}, ", end="")
        print(f"{round(drh['rms'],1)}, ", end="")

        sum_max += maxy
        sum_mean += drh["mean"]
        sum_rms += drh["rms"]
        npts += 1

    elif table_type == "temp":
        sif = 1
        maxt = max(abs(dtemp["max"]), abs(dtemp["min"]))
        maxws = max(abs(dwind["max"]), abs(dwind["min"]))
        print(f"\n{case['loc']}, {date_str[:10]},{date_str[10:13]}, ", end="")
        print(f"{round(maxt,sif)}, ", end="")
        print(f"{round(dtemp['mean'],sif)}, ", end="")
        print(f"{round(dtemp['rms'],sif)}, ", end="")
        print(f"{round(maxws,sif)}, ", end="")
        print(f"{round(dwind['mean'],sif)}, ", end="")
        print(f"{round(dwind['rms'],sif)} ", end="")

        sum_max_dt += maxt
        sum_mean_dt += dtemp["mean"]
        sum_rms_dt += dtemp["rms"]
        sum_max_dws += maxws
        sum_mean_dws += dwind["mean"]
        sum_rms_dws += dwind["rms"]
        npts += 1

if table_type == "rh":
    print("\nAverage, , , ", end="")
    print(f"{round(sum_max/npts,sif)}, ", end="")
    print(f"{round(sum_mean/npts,sif)}, ", end="")
    print(f"{round(sum_rms/npts,sif)}, ", end="")
elif table_type == "temp":
    print("\nAverage, , , ", end="")
    print(f"{round(sum_max_dt/npts,sif)}, ", end="")
    print(f"{round(sum_mean_dt/npts,sif)}, ", end="")
    print(f"{round(sum_rms_dt/npts,sif)}, ", end="")
    print(f"{round(sum_max_dws/npts,sif)}, ", end="")
    print(f"{round(sum_mean_dws/npts,sif)}, ", end="")
    print(f"{round(sum_rms_dws/npts,sif)} ", end="")
