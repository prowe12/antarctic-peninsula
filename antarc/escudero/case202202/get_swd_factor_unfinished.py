#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 16:12:26 2019

@author: prowe
"""

# Dependencies
import numpy as np
from datetime import datetime
import pysolar
import matplotlib.pyplot as plt
import pytz
import pandas as pd
import datetime as dt

# My modules
from antarc.getfilenames import get_filenames_and_dates

# Parameter modules
from antarc.escudero.parameters import esc_params, swd_params  # , esc202202


# Location and time-dependent parameters
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
dir_out = swd_params.dir_out
pyr_dir = swd_params.STAND_DIR
pyr_fileformat = swd_params.PYR_FILEFORMAT
outfile_sw_factor = swd_params.outfile_sw_factor
# samplefilename = swd_params.samplefilename

# Solar constant
Q = 1361  # Solar constant, W/m2

# Inputs
start_year = 2022
pyr_dir += str(start_year) + "/"


# Files
outfile_sw_factor = dir_out + outfile_sw_factor
outfile_npz = "pyranometer_" + str(start_year) + ".npz"
outfile_txt = "pyranometer_" + str(start_year) + ".txt"


# Set up for datetimes referenced to UTC
utc = pytz.utc

# Get filenames
files, dates = get_filenames_and_dates(pyr_dir, pyr_fileformat)

# Load the pyranometer data
pir = pd.read_csv(pyr_dir + files[0])

# Get solar zenith angles for all pyranometer data (slow)
# dtimes = pir[["Date", "Time"]].apply(pd.to_datetime)
# dtimes = dtimes["Time"].values
ntimes = len(pir)
fmt = "%Y-%m-%d%H:%M:%S"
dtimes = [
    dt.datetime.strptime(pir["Date"][i] + pir["Time"][i], fmt).replace(
        tzinfo=utc
    )
    for i in range(ntimes)
]

pir_sza = np.zeros(len(dtimes))
for i, date in enumerate(dtimes):
    pir_sza[i] = pysolar.solar.get_altitude(lat, lon, date)


# Get direct/diffuse radiation from pysolar
# i = np.arange(0, len(pir["date"]), 20)
# rad_dates = pir["date"][i]
# N_dates = len(rad_dates)

rad = np.zeros(ntimes)
direct = 0 * np.zeros(ntimes)
diffuse = 0 * np.zeros(ntimes)
sza = np.zeros(ntimes)
# sglobal=0* np.zeros(N_dates)
for i, date in enumerate(dtimes):
    # date = rad_dates[i]
    # -datetime.timedelta(3, hours) #pir['date'][i] #, tzinfo=datetime.timezone.utc)
    # date = datetime.datetime(2017, 11, 8, 17, 13, 1, tzinfo=datetime.timezone.utc)

    alt = pysolar.solar.get_altitude(lat, lon, date)

    # If you get the radiation when the sun is below the horizon, you get nonsense
    if alt > 0:
        rad[i] = pysolar.solar.radiation.get_radiation_direct(date, alt)
        diffuse[i] = pysolar.util.diffuse_underclear(lat, lon, date)
        direct[i] = pysolar.util.direct_underclear(lat, lon, date)
        # sglobal[i] = pysolar.util.global_irradiance_clear(latitude_deg, longitude_deg, thisdate)
    sza[i] = alt

plt.figure()
plt.plot(dtimes, sza)
plt.plot(dtimes, rad)
plt.plot(dtimes, diffuse)
plt.plot(dtimes, 1e-18 * direct)

# Sort and parametrize
isort = np.argsort(sza)
sza_sort = sza[isort]
diffuse_sort = diffuse[isort]
Qcos_theta = (Q * np.cos(np.deg2rad(90 - sza)))[isort]

ipos = np.where(sza_sort >= 0)[0]
p = np.polyfit(sza_sort[ipos], diffuse_sort[ipos], 1)
p_Q = np.polyfit(sza_sort[ipos], Qcos_theta[ipos], 1)
x = np.linspace(0, max(sza), 100)

plt.figure(0)
plt.cla()
plt.plot(sza_sort, diffuse_sort, label="diffuse")
# plt.plot(sza_sort, rad[isort], label = 'rad')
# plt.plot(sza_sort, direct[isort], label = 'direct')
plt.plot(x, np.polyval(p, x), label="fit, pysolar diffuse")
plt.plot(x, np.polyval(p_Q, x), label="fit, Qcos(theta)")
plt.plot(
    sza_sort[ipos],
    np.polyval(p, sza_sort[ipos]) - diffuse_sort[ipos],
    label="polyfit to diffuse",
)
plt.legend()

pir_model = np.polyval(p, pir_sza)
pir_model_Q = np.polyval(p_Q, pir_sza)

plt.figure(0)
plt.cla()
plt.plot(pir_sza, pir["Radiation"] / pir_model, ".", markersize=0.5)
plt.plot(pir_sza, pir["Radiation"] / pir_model_Q, ".", markersize=0.5)
plt.axis([0, 60, 0, 1])

# Plot as function of sza
plt.figure(1)
plt.cla()
plt.plot(sza, diffuse, ".", label="model, diffuse")
# plt.plot(sza, rad, '+', label = 'model, rad')
plt.plot(sza, Q * np.cos(np.deg2rad(90 - sza)), ".", label="Qcos(theta)")
plt.plot(pir_sza, pir["sw_radiation"], ".", markersize=0.3, label="meas")
plt.legend()


# .. Get pyranometer szas for times manually id'd as nearly clear
d1 = utc.localize(datetime(2017, 12, 2, 22, 53))
d2 = utc.localize(datetime(2017, 12, 2, 23, 8))
inds = np.where(np.logical_and(pir["date"] > d1, pir["date"] < d2))
iclear = inds
d1 = utc.localize(datetime(2017, 12, 2, 23, 22))
d2 = utc.localize(datetime(2017, 12, 2, 23, 24))
inds = np.where(np.logical_and(pir["date"] > d1, pir["date"] < d2))
iclear = np.union1d(iclear, inds)
d1 = utc.localize(datetime(2017, 12, 2, 23, 37))
d2 = utc.localize(datetime(2017, 12, 2, 23, 39))
inds = np.where(np.logical_and(pir["date"] > d1, pir["date"] < d2))
iclear = np.union1d(iclear, inds)
d1 = utc.localize(datetime(2017, 12, 7, 7, 46))
d2 = utc.localize(datetime(2017, 12, 7, 7, 48))
inds = np.where(np.logical_and(pir["date"] > d1, pir["date"] < d2))
iclear = np.union1d(iclear, inds)
d1 = utc.localize(datetime(2017, 12, 7, 8, 47))
d2 = utc.localize(datetime(2017, 12, 7, 9, 2))
inds = np.where(np.logical_and(pir["date"] > d1, pir["date"] < d2))
iclear = np.union1d(iclear, inds)
d1 = utc.localize(datetime(2017, 12, 7, 10, 48))
d2 = utc.localize(datetime(2017, 12, 7, 11, 18))
inds = np.where(np.logical_and(pir["date"] > d1, pir["date"] < d2))
iclear = np.union1d(iclear, inds)
d1 = utc.localize(datetime(2017, 12, 12, 20, 30))
d2 = utc.localize(datetime(2017, 12, 12, 20, 45))
inds = np.where(np.logical_and(pir["date"] > d1, pir["date"] < d2))
iclear = np.union1d(iclear, inds)

# Zeros and below
# d1 = utc.localize(datetime(2017,12, 3, 1, 54))
# d2 = utc.localize(datetime(2017,12, 3, 5, 25))
# inds = np.where(np.logical_and(pir['date'] > d1, pir['date'] < d2))
# iclear = np.union1d(iclear,inds)


pir_sza_clear = np.zeros(len(pir["date"][iclear]))
for i, thisdate in enumerate(pir["date"][iclear]):
    pir_sza_clear[i] = pysolar.solar.get_altitude(
        latitude_deg, longitude_deg, thisdate
    )

# Plot as function of sza
plt.figure(3)
plt.cla()
plt.plot(sza, diffuse, ".", label="model, diffuse")
plt.plot(sza_sort, Qcos_theta, ".", label="Qcos(theta)")
plt.plot(pir_sza, pir["sw_radiation"], ".", markersize=0.3, label="meas")
plt.plot(
    pir_sza[iclear],
    pir["sw_radiation"][iclear],
    "o",
    fillstyle="none",
    label="~clear",
)
# plt.plot(sza, rad, '+', label = 'model, rad')
plt.legend()


pir_model_clear = np.polyval(p, pir_sza[iclear])
pir_model_clear_Q = np.polyval(p_Q, pir_sza[iclear])

# Fit data by eye
x = [2, 5.5, 9.0, 14.0, 21.0, 30.0, 40.0, 51.0]
y = [0.56, 0.66, 0.75, 0.83, 0.89, 0.92, 0.91, 0.88]
y_Q = [0.46, 0.56, 0.65, 0.73, 0.79, 0.82, 0.81, 0.78]

# Fit to quartic
pfac = np.polyfit(x, y, 4)
pfac_Q = np.polyfit(x, y, 4)

plt.figure(4)
plt.cla()
plt.plot(
    pir_sza,
    pir["sw_radiation"] / pir_model,
    ".",
    markersize=0.5,
    label="pysolar",
)
plt.plot(
    pir_sza,
    pir["sw_radiation"] / pir_model_Q,
    ".",
    markersize=0.5,
    label="Qcos(theta)",
)
plt.plot(
    pir_sza[iclear],
    pir["sw_radiation"][iclear] / pir_model_clear,
    "o",
    fillstyle="none",
    label="~clear (pysolar)",
)
plt.plot(
    pir_sza[iclear],
    pir["sw_radiation"][iclear] / pir_model_clear_Q,
    "o",
    fillstyle="none",
    label="~clear, (Qcos(theta))",
)
plt.plot(x, y, "o", label="selected points, pysolar")
plt.plot(x, y_Q, "o", label="selected points, Qcos(theta)")
plt.plot(
    sza_sort[ipos],
    np.polyval(pfac, sza_sort[ipos]),
    "r",
    label="polyval, pysolar",
)
plt.plot(
    sza_sort[ipos],
    np.polyval(pfac_Q, sza_sort[ipos]),
    "m",
    label="polyval, (Qcostheta))",
)
plt.axis([0, 60, 0, 1])
plt.legend()


# Plot the radiance
fig = plt.figure(2)
plt.cla()
plt.plot(pir["date"], pir["sw_radiation"])
plt.plot(rad_dates, diffuse, label="diffuse")
plt.plot(rad_dates, diffuse * np.polyval(pfac, sza), "r")
# plt.plot(rad_dates, .89*diffuse, label='diffuse')
# plt.plot(rad_dates, .8*diffuse, label='diffuse')
# plt.plot(rad_dates, .7*diffuse, label='diffuse')
plt.ylabel("SW flux (W/m${^2}$)")
fig.autofmt_xdate()

# 2017/02/06 - fairly smooth for part of day
# 2017/12/27 - fairly smooth


# .. Save pfac
f = open(outfile_sw_factor, "w")
f.write(
    str(pfac[0])
    + "    "
    + str(pfac[1])
    + "    "
    + str(pfac[2])
    + "    "
    + str(pfac[3])
    + "    "
    + str(pfac[4])
)
f.close()


"""
start_date = datetime(2018,11,30)
end_date = datetime(2018,12,4)

start_date = datetime(2017,11,5)
end_date = datetime(2017,12,31)

# .. Get the clearest day
# Clear day: 2017/12/27
iclear = np.intersect1d(np.where(pir['date'] >= pir['date'][154093]),
                        np.where(pir['date'] <= pir['date'][156670]))



# Get direct/diffuse radiation from pysolar
ipir = np.intersect1d(np.where(pir['date'] >= pir['date'][100000]),
                      np.where(pir['date'] <= pir['date'][156670]))

rad_dates = pir['date'][ipir]
N_dates = len(ipir)
diffuse = np.zeros(N_dates)
#direct = 0 * np.zeros(N_dates)
#rad = np.zeros(N_dates); sglobal=0* np.zeros(N_dates)
for i in range(N_dates):
    thisdate = rad_dates[i] 
    
    altitude_deg = pysolar.solar.get_altitude(latitude_deg, longitude_deg, thisdate)
    # The radiation calculation when the sun is below the horizon is nonsense
    if altitude_deg>0:
        diffuse[i] = pysolar.util.diffuse_underclear(latitude_deg, longitude_deg, thisdate)
        #rad[i] = pysolar.solar.radiation.get_radiation_direct(thisdate, altitude_deg)
        #direct[i] = pysolar.util.direct_underclear(latitude_deg, longitude_deg, thisdate)
        #sglobal[i] = pysolar.util.global_irradiance_clear(latitude_deg, longitude_deg, thisdate)


# Plot it
fig = plt.figure(1); plt.cla()
plt.plot(pir['date'], pir['sw_radiation'])
plt.plot(pir['date'][iclear], pir['sw_radiation'][iclear])
plt.plot(pir['date'][ipir], diffuse, label='diffuse')
plt.ylabel('SW flux (W/m${^2}$)')
#plt.ylim([-5, 2200])
plt.xlim([start_date, end_date])
fig.autofmt_xdate()


# Compare SW to clear-sky from pysolar
startdate = datetime(2017,1,1); enddate = datetime(2018,1,5)


fig = plt.figure(3); plt.cla()
#plt.plot(rad_dates, rad, 'o-', fillstyle='none', label='direct')
plt.plot(pir['date'][ipir], diffuse, label='diffuse')
plt.plot(pir['date'], pir['sw_radiation'], label='pir')
plt.axis([pir['date'][0], pir['date'][-1], 0, 1500])
plt.ylabel('Shortwave radiation (W/m${^2}$)')
fig.autofmt_xdate(); plt.legend()

plt.figure(4)
plt.plot(pir['date'][ipir], diffuse-pir['sw_radiation'][ipir])



datecut = pir['date'][110657:]
radcut = pir['sw_radiation'][110657:]
radrepeat = np.hstack([pir['sw_radiation'][iclear], 
                       pir['sw_radiation'][iclear],
                       pir['sw_radiation'][iclear], 
                       pir['sw_radiation'][iclear]])

plt.figure(3); plt.cla()
plt.plot(datecut[:len(radrepeat)], radrepeat)
plt.plot(datecut, radcut)

plt.figure(4); plt.cla()
plt.plot([300, 334, 361, 365], [670, 825, 864, 863],'.-')

"""
