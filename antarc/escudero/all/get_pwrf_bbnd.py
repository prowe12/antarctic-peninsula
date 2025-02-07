#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:03:39 2024

@author: prowe

Make sure my method for getting the Polar WRF nearest gridpoint works

"""

# Dependencies
from netCDF4 import Dataset
import datetime as dt
import numpy as np


# My modules
from antarc import params
from antarc.pwrf_rsrc import get_inds_to_latlon
from antarc.escudero.parameters import esc_params


def get_pwrf_key(fname: str, lat: float, lon: float, key: str):
    with Dataset(fname) as nci:
        i, j = get_inds_to_latlon(
            nci["XLAT"][:].data, nci["XLONG"][:].data, lat, lon
        )
        res = nci[key][:, i, j]
    return res


def get_pwrf_bbnd(direc: str, dt0, fmt, lat, lon):
    fname_lwd = f"{direc}wrfout_LWD_d03_{dt0.strftime(fmt)}.nc"
    fname_swd = f"{direc}wrfout_SWD_d03_{dt0.strftime(fmt)}.nc"

    dtday = [dt.datetime(dt0.year, dt0.month, dt0.day, x) for x in range(24)]
    lwd = get_pwrf_key(fname_lwd, lat, lon, "GLW")
    swd = get_pwrf_key(fname_swd, lat, lon, "SWDOWN")
    return dtday, lwd, swd


def get_pwrf_bbnd_daterange(start_date, end_date) -> dict[list]:
    # Directories
    MEAS_DIR = params.MEAS_DIR
    dir_pwrf = MEAS_DIR + "Escudero/pwrf/broadband/"

    lat = esc_params.LATITUDE
    lon = esc_params.LONGITUDE
    fmt = "%Y%m%d"
    days = np.round((end_date - start_date).total_seconds() / 60 / 60 / 24)
    days = int(days)

    dtime = []
    lwd = []
    swd = []
    for day in range(days):
        dtday = start_date + dt.timedelta(days=day)

        dtime0, lwd0, swd0 = get_pwrf_bbnd(dir_pwrf, dtday, fmt, lat, lon)

        dtime += dtime0
        lwd += list(lwd0)
        swd += list(swd0)

    pwrf = {"date": dtime, "lwd": lwd, "swd": swd}
    pwrf["tot"] = [x + y for x, y in zip(swd, lwd)]
    return pwrf


# dir_lwdclear = MEAS_DIR + "Escudero/lwd_clear/"
# out_dir = params.PROJECT_DIR + "case_studies/May_2022/figures/"
# era_broadband_dir = f"{params.MEAS_DIR}Escudero/era5/broadband/"


# hour = 23
# filename = f"wrfout_d03_2022-05-14_{hour}_00_00.nc"
# ncall = Dataset(dir_pwrf + filename)
# print(ncall.variables.keys())
# jall, iall = ll_to_xy(ncall, lat, lon, meta=False)

# plt.figure(num=2, clear=True)
# plt.subplot(211)
# plt.plot([x.hour for x in dtime], swd)
# plt.plot(hour, ncall["SWDOWN"][0, iall, jall], "s")
# plt.plot(hour, ncall["GSW"][0, iall, jall], "*")
# plt.subplot(212)
# plt.plot([x.hour for x in dtime], lwd)
# plt.plot(hour, ncall["GLW"][0, iall, jall], "*")

# ncall.close()
