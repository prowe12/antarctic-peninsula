#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 11:22:25 2023

@author: prowe
"""

# Dependencies
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# My modules
from antarc import params
from antarc.escudero.all.analysis_rsrc import (
    get_era5_esc_utc,
    get_lwd_qc,
    get_swd_qc,
    get_files,
)
from antarc.escudero.all.results.get_era5_cbh import get_era5_cbh


from antarc.escudero.parameters import esc_params
from antarc.escudero.parameters import swd_params, lwd_params
from antarc.escudero.all.results.params import SAVEFIGS


def model_meas_matches(model_time, model_val, meas_time, meas_val):
    """
    Given a set of model values and times and measured values and times,
    return a new set of model values and measured values at the same times,
    where one measured value is selected closest in time to each model value,
    whereever measured values exist within 30 minutes of model values.
    Also return the new times.
    """
    inds0 = [i for i in range(len(model_time)) if model_time[i] > meas_time[0]]
    inds = [i for i in inds0 if model_time[i] < meas_time[-1]]

    ime = 0
    mod_match_val = []
    meas_match_val = []
    dtime_match = []
    for dtime, val in zip(model_time[inds], model_val[inds]):

        # If model does not have a cloud base height, skip
        if np.ma.is_masked(val) or np.isnan(val):
            continue

        # Skip the meas time ahead until it is just higher than the model time
        while meas_time[ime] < dtime:
            if ime >= len(meas_time) - 1:
                break
            ime += 1

        # Check meas time before and after, if either is within 30 minutes,
        # use closest
        tdiff_bef = (dtime - meas_time[ime - 1]).total_seconds()
        tdiff_aft = (meas_time[ime] - dtime).total_seconds()

        # If both are within 30 minutes, use weighted mean
        bef_real = np.isnan(meas_val[ime - 1]) == False
        aft_real = np.isnan(meas_val[ime]) == False
        if tdiff_bef < 1800 and tdiff_aft < 1800 and bef_real and aft_real:
            mod_match_val.append(val)
            tdiff = tdiff_bef + tdiff_aft
            wt_bef = tdiff_aft / tdiff
            wt_aft = tdiff_bef / tdiff
            if wt_bef < 0 or wt_aft < 0 or wt_bef + wt_aft != 1:
                raise ValueError("Bad weights!")
            val0 = wt_aft * meas_val[ime] + wt_bef * meas_val[ime - 1]
            meas_match_val.append(val0)
            dtime_match.append(dtime)
        elif tdiff_bef < 1800 and bef_real:
            mod_match_val.append(val)
            meas_match_val.append(meas_val[ime - 1])
            dtime_match.append(dtime)
        elif tdiff_aft < 1800 and aft_real:
            mod_match_val.append(val)
            meas_match_val.append(meas_val[ime])
            dtime_match.append(dtime)

        if tdiff_bef < 0 or tdiff_aft < 0 and ime != len(meas_time) - 1:
            raise ValueError("Times are not as expected")

    return dtime_match, meas_match_val, mod_match_val


def get_days_since_date(dtimes, start):
    """
    For an input list of datetimes, return a numpy array with the
    number of fractional days that have passed since a given day
    """
    return np.array(
        [(x - start).total_seconds() / 60 / 60 / 24 for x in dtimes]
    )


# def load_libradtran(dtimes):
#     """
#     Load pre-run libradtran results for clear sky downwelling shortwave
#     for the input dates indicates. The folders are by year, files by day
#     """
#     swd_dates = []
#     swd_clear = []

#     librad_dir = "/Users/prowe/Sync/measurements/Escudero/swd_clear/"

#     # Get all libradtran data from first day to last day
#     dtime = dtimes[0].replace(hour=0, minute=0, second=0)
#     while dtime < dtimes[-1]:
#         daystr = dtime.strftime("%Y%m%d")
#         clrfile = f"{librad_dir}{dtime.year}/libradtran_{daystr}.csv"
#         try:
#             libclear = pd.read_csv(clrfile, delimiter=",")
#         except:
#             print(f"{clrfile} missing")
#             return None, None
#         swd_dates += pd.to_datetime(libclear["Date"]).to_list()
#         swd_clear += libclear["global"].to_list()
#         dtime += dt.timedelta(days=1)

#     return swd_dates, swd_clear


def load_lidar_coltype(lidar_dir, fmt, date1, date2, verbose=False):
    """
    Get the lidar columntype for a specified date range

    Long_Name: Column Data Mask
    (
     0 = obscured column,
     1 = clear column,
     2 = unidentified cloud column,
     3 = liquid column,
     4 = ice column
     )
    """
    base_ht_subset = True
    nswitch = 0

    # If dates are offset-aware, convert to offset-naive
    date1 = date1.replace(tzinfo=None)
    date2 = date2.replace(tzinfo=None)

    # Get a list of the files for the selected time period
    files_daterange = get_files(lidar_dir, fmt, date1, date2)

    dtime = []
    coltype_temp = []
    coltype_alt = []
    alt = [400, 2000, 5000]
    reptemps_ice = []

    for fname in files_daterange:
        if verbose:
            print(fname)
        ymd = dt.datetime.strptime(fname, fmt).replace(hour=0, minute=0)
        with Dataset(lidar_dir + fname) as ncid:

            hour = ncid["DataTime"][:].data
            hour = hour[:].data
            coltype_file_alt = 10 * ncid["ColumnType"][:].data
            coltype_file_temp = 10 * ncid["ColumnType"][:].data

            cbalt = ncid["CloudBaseAltitude"][:].data
            cbtemp = ncid["CloudBaseTemperature"][:].data

            # Check for hours that are nan?
            # coltype_file_alt = 10 * ncid["ColumnType"][~np.isnan(hour)].data
            # coltype_file_temp = 10 * ncid["ColumnType"][~np.isnan(hour)].data
            # hour = hour[~np.isnan(hour)].data

            # Find cases where cloud was identified by altitude and
            # temperature are nan. These are cases where there was not
            # clear air under the clouds. For these, make the lowest
            # cloud base height the altitude and get the temperature
            icld = np.where(coltype_file_alt >= 20)[0]
            inan = np.where(np.isnan(cbtemp))[0]
            irep = np.intersect1d(icld, inan)

            # Get cloud alts and temperatures for the lowest clouds
            # which weren't previously classified
            if len(irep) > 0:
                ialt = [
                    np.where(ncid["CloudMask"][i, :] == 1)[0][0] for i in irep
                ]
                temp = ncid["Temperature"][:].data
                cbtemp[irep] = [
                    temp[irep[i], ialt[i]] for i in range(len(irep))
                ]
                cbalt[irep] = ncid["Range"][ialt]

                # coltype0 = coltype_file_temp[irep]
                # reptemps = cbtemp[irep]
                # reptemps_liq += list(reptemps[coltype0 == 30])
                # reptemps_ice += list(reptemps[coltype0 == 40])

            # Find cases where ice was above 276.15 K and reclassify as liquid
            iwarm = np.where(coltype_file_temp == 40)[0]
            iwarm = np.intersect1d(iwarm, np.where(cbtemp > 276.15)[0])
            if len(iwarm) > 0:
                coltype_file_temp[iwarm] = 30.5
                coltype_file_alt[iwarm] = 30.5
                nswitch += len(iwarm)
                reptemps_ice += list(cbtemp[iwarm])

            if any(np.isnan(hour)):
                raise ValueError("bad hour in data")

            # Indices for cloud base temperature ranges
            i273 = np.where(cbtemp >= 273.15)[0]
            i260 = np.where(cbtemp >= 260)[0]
            i260 = np.setxor1d(i273, i260)
            i240 = np.where(cbtemp >= 240)[0]
            i240 = np.setxor1d(i260, i240)
            i240 = np.setxor1d(i273, i240)

            # Indices for cloud base height ranges
            ilowest = np.where(cbalt <= alt[0])[0]
            ilo = np.where(cbalt <= alt[1])[0]
            ilo = np.setxor1d(ilowest, ilo)
            imed = np.where(cbalt <= alt[2])[0]
            imed = np.setxor1d(ilo, imed)
            imed = np.setxor1d(ilowest, imed)

            coltype_file_temp[i273] += 3
            coltype_file_temp[i260] += 2
            coltype_file_temp[i240] += 1

            coltype_file_alt[ilowest] += 3
            coltype_file_alt[ilo] += 2
            coltype_file_alt[imed] += 1

            if base_ht_subset:
                # Only use data with cloud heights
                ino_cloud_base = np.where(np.isnan(cbtemp))[0]
                ino_cloud_base = np.intersect1d(ino_cloud_base, icld)
                # coltype_file_temp[ino_cloud_base] += 5
                if np.any(ino_cloud_base):
                    raise ValueError("Got unexpected nan")

                ino_cloud_base = np.where(np.isnan(cbalt))[0]
                ino_cloud_base = np.intersect1d(ino_cloud_base, icld)
                if np.any(ino_cloud_base):
                    raise ValueError("Got unexpected nan")
                # coltype_file_alt[ino_cloud_base] += 5

            # Make sure same length
            if len(hour) != len(coltype_file_alt) or len(hour) != len(
                coltype_file_temp
            ):
                raise ValueError("Dates and column types should be same len")

            # Convert hour to datetime
            dtime_file = [ymd + dt.timedelta(hours=x) for x in hour]
            dtime += dtime_file
            coltype_temp += list(coltype_file_temp)
            coltype_alt += list(coltype_file_alt)

            # After converting:
            # -10: Data missing
            # 0: obscured
            # 10: clear
            # 20-23: unidentified cloud
            # 30-33: liq
            # 40-43: ice
    return dtime, coltype_temp, coltype_alt, alt, reptemps_ice


# Directories and variables needed from main param file
meas_dir = params.MEAS_DIR

# Directories and variables needed from other param files
lwd_base_dir = lwd_params.LWD_DIR
era_dir = f"{params.MEAS_DIR}Escudero/era5/broadband/"
lwd_fmt = lwd_params.LWD_FILEFORMAT
swd_fmt = swd_params.SWD_FILEFORMAT
lat = esc_params.LATITUDE
lon = esc_params.LONGITUDE
# lwd_clear_dir = meas_dir + "Escudero/lwd_clear/"
era_fmt = "era5_esc_broadband_*[0-9]*.nc"
lidar_dir = meas_dir + "Escudero/mpl/reprocessing_stillwell/v3/"
lidar_fmt = "KGI_MPLData_%Y%m%d.cdf"


# Flags
savedir = "antarc/escudero/all/figures/"

make_figs = False
plot_lwd_hists_for_clear_obsc_unk = False
plot_lwd_hists_for_alts = False
plot_lwd_hists_for_alts_inc_era5 = False
plot_lwforcing_hists_for_clear_liq_ice_obsc = True
plot_lwforcing_hists_for_temps = True
plot_lwforcing_hists_for_temps_inc_era5 = False
plot_lwforcing_hists_for_alts = True
plot_lwforcing_hists_for_alts_inc_era5 = False


# Main code


# Date range: do everything timezone-UNAWARE - all dates must be in UTC
date1 = dt.datetime(2017, 1, 1)
date2 = dt.datetime(2023, 12, 31)

# Load in the ERA5 broadband data, one year at a time
dt1 = dt.datetime(2017, 1, 1)
dt2 = dt.datetime(2017, 12, 31)
era5 = get_era5_esc_utc(era_dir, era_fmt, dt1, dt2, lat, lon)
for year in [2018, 2019, 2020, 2021, 2022, 2023]:
    dt1 = dt.datetime(year, 1, 1)
    dt2 = dt.datetime(year, 12, 31)
    era5b = get_era5_esc_utc(era_dir, era_fmt, dt1, dt2, lat, lon)
    for key in era5.keys():
        era5[key] = np.hstack([era5[key], era5b[key]])


# Load in the Escudero broadband measurements (slow!)
lwd_date, lwd = get_lwd_qc(lwd_base_dir, lwd_fmt, date1, date2)
swd_date, swd = get_swd_qc(swd_params.SWD_DIR, swd_fmt, date1, date2)


# ColumnType:
#   float32 ColumnType(Time)
#   Unit: Unitless
#   Long_Name: Column Data Mask (
#     0 = obscured column,
#     10 = clear column,
#     20 = unidentified cloud column,
#     30 = liquid column,
#     40 = ice column)
# Get the miniMPL cloud mask for the column for different cloud base heights
(
    mpl_date,
    col,
    mpl_col_alt,
    alt,
    reptemps_ice,
) = load_lidar_coltype(lidar_dir, lidar_fmt, date1, date2)
col = np.array(col)

col[col == 21] = 20
col[col == 22] = 20
col[col == 23] = 20
col_for_temps = col


col = np.array(mpl_col_alt)
col[col == 21] = 20
col[col == 22] = 20
col[col == 23] = 20
col_heights = col


bins = np.arange(260, 280, 1) + 0.15
plt.figure()
hist = plt.hist(reptemps_ice, bins)
plt.plot([273.15, 273.15], [0, 10000], "k:")
plt.plot([276.15, 276.15], [0, 10000], "k:")
plt.ylim([0, 10000])
plt.xlabel("Temperature (K)")

# Get dates as days past the first day of the given month
# Remember: ERA5 accumulates to end of hour: move back 1/2 hour
# Interpolate ERA5 clear onto lwd_date vector
day1 = dt.datetime(lwd_date[0].year, lwd_date[0].month, 1)
lwd_days = get_days_since_date(lwd_date, day1)
swd_days = get_days_since_date(swd_date, day1)
era_days = get_days_since_date(era5["date"], day1) - 0.5 / 24
mpl_days = get_days_since_date(mpl_date, day1)

# Interpolate era5 clear sky onto LWD grid
# ERA is hourly
era_lwd_clr = np.interp(lwd_days, era_days, era5["lwd_clr"])

# Calculate cloud forcing for longwave downwelling flux
lwd_force = lwd - era_lwd_clr

# Interpolate forcing onto mpl grid
# Forcing is by minute, based on measured DLW flux by minute and
# clear-sky DLW flux by hour from ERA5 interpolated to the minute
# MPL is every two minutes.
lwd_force_mpl = np.interp(
    mpl_days, lwd_days, lwd_force, left=np.nan, right=np.nan
)

col = np.floor(col_for_temps)
iliq = np.where(col >= 30)[0]
nliq = len(np.intersect1d(iliq, np.where(col < 40)[0]))
nice = len(np.where(col >= 40)[0])
nobsc = len(np.where(col == 0)[0])
nclr = len(np.where(col == 10)[0])

# Correct assuming same proportion for obscured
ntot = nliq + nice
nliq += nobsc * (nliq / ntot)
nice += nobsc * (nice / ntot)

ntot = nliq + nice + nclr

if plot_lwd_hists_for_clear_obsc_unk:
    frc = lwd_force_mpl
    col = col_heights

    # 3-panel histograms of LWD for clear, obscured, and unknown clouds
    bins = np.arange(0, np.nanmax(frc), 10)

    plt.figure(num=1, figsize=[6, 8], clear=True)
    plt.subplot(3, 1, 1)
    plt.hist(frc[col == 10], bins)
    plt.legend(["clear"])
    plt.ylabel("Number")

    plt.subplot(3, 1, 2)
    plt.hist(frc[col == 0], bins)
    plt.legend(["obscured"])
    plt.ylabel("Number")

    plt.subplot(3, 1, 3)
    plt.hist(frc[col == 20], bins, stacked=True)
    plt.legend(["unknown cloud"])
    plt.ylabel("Number")
    plt.xlabel("Downwelling longwave cloud forcing (W m$^{-2}$)")

    if SAVEFIGS:
        plt.savefig(savedir + "forcing_hists1.png", format="png")
        plt.savefig(savedir + "forcing_hists1.eps", format="eps")


if plot_lwd_hists_for_alts:
    frc = lwd_force_mpl
    col = col_heights
    # ilowest = np.where(cbalt <= alt[0])[0] => +3
    # ilo = np.where(cbalt <= alt[1])[0]
    # ilo = np.setxor1d(ilowest, ilo)  => +2
    # imed = np.where(cbalt <= alt[2])[0]
    # imed = np.setxor1d(ilo, imed)
    # imed = np.setxor1d(ilowest, imed) => +1
    # So if alts = [400, 2000, 5000], we have
    # 30 => liq < 400 m
    # 31 => 400 -2000 m
    # 32 => 2000-5000 m
    # 33 => > 5000 m

    bins = np.arange(0, np.nanmax(frc), 10)

    histcolorsl = ["cyan", "blue", "orange", "green"]
    histcolorsi = ["cyan", "blue", "orange"]

    # Since we start with 30 and combine 31 and 32 below
    legl = [
        f"base > {alt[2]}m",
        f"{alt[0]}m < base < {alt[2]}m",
        f"base < {alt[0]}m",
        f"base < {alt[0]}m; reclassified ice",
    ]
    legi = [
        f"base > {alt[0]}m",
        f"{alt[0]}m < base < {alt[1]}m",
        f"base < {alt[2]}m",
    ]
    # fmt: off
    fig, [ax1, ax2] = plt.subplots(2, 1, num=20, clear=True) # , figsize=[8, 8]
    ax1.set_position([0.14, 0.57, 0.83, 0.4])
    ax2.set_position([0.14, 0.1, 0.83, 0.4])
    
    data = [frc[col == 30], np.hstack([frc[col == 31], frc[col == 32]]), frc[col == 33], frc[col==33.5]]
    ax1.hist(data, bins, stacked=True, color=histcolorsl)
    ax1.legend(legl)
    ax1.set_ylabel("Number")
    
    data = [frc[col == 40], np.hstack([frc[col == 41], frc[col == 42]]), frc[col == 43]]
    ax2.hist(data, bins, stacked=True, color=histcolorsi)
    ax2.legend(legi)
    
    ax2.set_ylabel("Number")
    ax2.set_xlabel("Downwelling longwave cloud forcing (W m$^{-2}$)")
    
    ax1.text(2, 25000, 'a) Liquid-containing')
    ax2.text(2, 2500, 'b) Ice-only')
    ax1.set_xlim([0, 100])
    ax2.set_xlim([0, 100])
    # fmt: on

    if SAVEFIGS:
        plt.savefig(savedir + "forcing_hists2.png", format="png", dpi=600)
        plt.savefig(savedir + "forcing_hists2.eps", format="eps")


if plot_lwforcing_hists_for_clear_liq_ice_obsc:
    frc = lwd_force_mpl
    col = col_heights

    # Get the total number of points excluding missing data within mpl files
    # (given that we are not accounting for missing data when there was no
    #  mpl file). That is, count all values but -15 (no data)
    tot_pts = len(frc) - len(frc[col < 0])

    wid = 0.8
    height = 0.2
    left = 0.15
    bot1 = 0.08
    bot2 = bot1 + height + 0.03
    bot3 = bot2 + height + 0.03
    bot4 = bot3 + height + 0.03

    bins = np.arange(0, np.nanmax(frc), 10)

    # Plot
    fig, axs = plt.subplots(4, 1, num=1, figsize=[6, 6], clear=True)
    axs = axs.flatten()

    axs[0].set_position([left, bot4, wid, height])
    axs[1].set_position([left, bot3, wid, height])
    axs[2].set_position([left, bot2, wid, height])
    axs[3].set_position([left, bot1, wid, height])

    # fmt: off
    ax = axs[0]
    pct = len(frc[col == 10]) / tot_pts * 100
    ax.hist(frc[col == 10], bins)
    ax.legend([f"clear {pct:1.0f}%"], loc="upper left", bbox_to_anchor=[0.15, 0.5, 0.5, 0.5])
    ax.set_ylabel("Number")

    ax = axs[1]
    data = np.hstack([frc[col == 30], frc[col == 31], frc[col == 32], frc[col == 33]])
    ax.hist(data, bins, stacked=True)
    ax.set_ylabel("Number")

    ax = axs[2]
    data = np.hstack(
        [frc[col == 40], frc[col == 41], frc[col == 42], frc[col == 43]]
    )
    pct = len(data) / tot_pts * 100
    pct_tot = (len(data) + len(frc[col == 45])) / tot_pts * 100
    ax.hist(data, bins)
    ax.set_ylabel("Number")
    ax.legend([f"Ice {pct:1.0f}% ({pct_tot:1.0f}%)"], loc="upper left")

    ax = axs[3]
    pct = len(frc[col == 0]) / tot_pts * 100
    ax.hist(frc[col == 0], bins)
    ax.legend([f"Obscured {pct:1.0f}%"], loc="upper left")
    ax.set_ylabel("Number")
    ax.set_xlabel("Downwelling longwave cloud forcing (W/m$^{2}$)")
    # fmt: on

    for ax in axs:
        ax.set_xlim([0, 100])
    for ax in axs[:-1]:
        ax.set_xticklabels([])

    # Same but without separating by height
    tot_pts = len(frc) - len(frc[col < 0])

    wid = 0.8
    height = 0.2
    left = 0.15
    bot1 = 0.08
    bot2 = bot1 + height + 0.03
    bot3 = bot2 + height + 0.03
    bot4 = bot3 + height + 0.03

    bins = np.arange(0, np.nanmax(frc), 10)

    # Plot for paper
    fig, axs = plt.subplots(4, 1, num=1, figsize=[6, 6], clear=True)
    axs = axs.flatten()

    axs[0].set_position([left, bot4, wid, height])
    axs[1].set_position([left, bot3, wid, height])
    axs[2].set_position([left, bot2, wid, height])
    axs[3].set_position([left, bot1, wid, height])

    # fmt: off
    ax = axs[0]
    pct = len(frc[col == 10]) / tot_pts * 100
    ax.hist(frc[col == 10], bins)
    ax.legend([f"clear {pct:1.0f}%"], loc="upper left", bbox_to_anchor=[0.15, 0.5, 0.5, 0.5])
    ax.set_ylabel("Number")
   
    ax = axs[1]
    data = np.hstack([frc[col == 30], frc[col == 31], frc[col == 32], frc[col == 33], frc[col == 35]])
    pct = len(data) / tot_pts * 100
    ax.hist(data, bins)
    ax.legend([f"Liquid-containing {pct:1.0f}%"])
    ax.set_ylabel("Number")
   
    ax = axs[2]
    data = np.hstack(
        [frc[col == 40], frc[col == 41], frc[col == 42], frc[col == 43], frc[col==45]]
    )
    pct = len(data) / tot_pts * 100
    pct_tot = (len(data) + len(frc[col == 45])) / tot_pts * 100
    ax.hist(data, bins)
    ax.set_ylabel("Number")
    ax.legend([f"Ice {pct:1.0f}%"], loc="upper left")
   
    ax = axs[3]
    pct = len(frc[col == 0]) / tot_pts * 100
    unk =  np.hstack([frc[col == 20], frc[col == 25]])
    pct_unk = len(unk) / tot_pts * 100
    ax.hist([frc[col == 0], unk], bins, stacked=True)
    ax.legend([f"Obscured {pct:1.0f}%", f"Unknown {pct_unk:1.0f}%"], loc="upper left")
    ax.set_ylabel("Number")
    ax.set_xlabel("Downwelling longwave cloud forcing (W m$^{-2}$)")
    # fmt: on

    for ax in axs:
        ax.set_xlim([0, 100])
    for ax in axs[:-1]:
        ax.set_xticklabels([])

    axs[0].text(95, 10000, "a)")
    axs[1].text(95, 120000, "b)")
    axs[2].text(95, 30000, "c)")
    axs[3].text(95, 30000, "d)")

    if SAVEFIGS:
        plt.savefig(savedir + "forcing_hists1b.png", format="png", dpi=600)
        plt.savefig(savedir + "forcing_hists1b.eps", format="eps")

if plot_lwforcing_hists_for_temps:

    col = np.floor(col_for_temps)
    frc = lwd_force_mpl

    bins = np.arange(0, np.nanmax(frc), 10)
    histcolorsl = ["cyan", "blue", "green", "orange", "red"]
    histcolorslflip = ["red", "orange", "green", "blue", "cyan"]

    histcolors = ["cyan", "blue", "green", "orange"]
    histcolorsflip = ["orange", "green", "blue", "cyan"]

    # "T $\geq$ 273 K: reclassified ice",
    legl = [
        "T $\geq$ 273 K",
        "260 K $\leq$ T < 273 K",
        "240 K $\leq$ T < 260 K",
        "T < 240 K",
    ]
    legi = [
        "273 K $\leq$ T < 276 K",
        "260 K $\leq$ T < 273 K",
        "240 K $\leq$ T < 260 K",
        "T < 240 K",
    ]

    # fmt: off
    fig, [ax1, ax2] = plt.subplots(2, 1, num=3, clear=True)
    ax1.set_position([0.15, 0.55, 0.75, 0.4])
    ax2.set_position([0.15, 0.08, 0.75, 0.4])

    # data = [frc[col == i] for i in [30, 31, 32, 33, 33.5]]
    # dataflip = [frc[col == i] for i in [33.5, 33, 32, 31, 30]]
    # ax1.hist(dataflip, bins, stacked=True, color=histcolorslflip)
    # ax1.legend(legl)
    # ax1.hist(data, bins, stacked=True, color=histcolorsl)
    data = [frc[col == i] for i in [30, 31, 32, 33]]
    dataflip = [frc[col == i] for i in [33, 32, 31, 30]]
    ax1.hist(dataflip, bins, stacked=True, color=histcolorsflip)
    ax1.legend(legl)
    ax1.hist(data, bins, stacked=True, color=histcolors)
    ax1.set_ylabel("Number")
    ax1.text(5, 30000,'a) Liquid-Containing')

    data = [frc[col == i] for i in [40, 41, 42, 43]]
    dataflip = [frc[col == i] for i in [43, 42, 41, 40]]
    hist2 = ax2.hist(dataflip, bins, stacked=True, color=histcolorsflip)
    ax2.legend(legi)
    ax2.hist(data, bins, stacked=True, color=histcolors)
    ax2.set_ylabel("Number")
    ax2.text(5, 6000, 'b) Ice') #, fontsize=12)
    ax2.set_xlabel("Downwelling longwave cloud forcing (W m$^{-2}$)")
    
    ax1.set_xlim([0,100])
    ax2.set_xlim([0,100])

    if SAVEFIGS:
        # plt.savefig(savedir + "forcing_hists_bytemp.png", format="png", dpi=600)
        # plt.savefig(savedir + "forcing_hists_bytemp.eps", format="eps")
        plt.savefig(savedir + "forcing_hists_bytemp2.png", format="png", dpi=600)
        plt.savefig(savedir + "forcing_hists_bytemp2.eps", format="eps")
    # fmt: on

    # Note: total number of points for ice is:
    np.sum(hist2[0][-1, :])

    # Get means and medians
    data = [frc[col == i] for i in [30, 31, 32, 33]]
    print(f"Mean for liq, T<240: {np.nanmean(data[0])}")
    print(f"Mean for liq, 240<T<260: {np.nanmean(data[1])}")
    print(f"Mean for liq, 260>T<273: {np.nanmean(data[2])}")
    print(f"Mean for liq, T>273: {np.nanmean(data[3])}")
    print()

    data = [frc[col == i] for i in [40, 41, 42, 43]]
    print(f"Mean for ice, T<240: {np.nanmean(data[0])}")
    print(f"Mean for ice, 240<T<260: {np.nanmean(data[1])}")
    print(f"Mean for ice, 260>T<273: {np.nanmean(data[2])}")
    print(f"Mean for ice, T>273: {np.nanmean(data[3])}")

    data = [frc[col == i] for i in [30, 31, 32, 33]]
    print(f"Median for liq, T<240: {np.nanmedian(data[0])}")
    print(f"Median for liq, 240<T<260: {np.nanmedian(data[1])}")
    print(f"Median for liq, 260>T<273: {np.nanmedian(data[2])}")
    print(f"Median for liq, T>273: {np.nanmedian(data[3])}")
    print()

    data = [frc[col == i] for i in [40, 41, 42, 43]]
    print(f"Median for ice, T<240: {np.nanmedian(data[0])}")
    print(f"Median for ice, 240<T<260: {np.nanmedian(data[1])}")
    print(f"Median for ice, 260>T<273: {np.nanmedian(data[2])}")
    print(f"Median for ice, T>273: {np.nanmedian(data[3])}")

    # Overall:
    data = np.hstack([frc[col == i] for i in [32, 33]])
    print(f"Overall mean liq, T>260: {np.nanmean(data)}")
    print(f"Overall median liq, T>260: {np.nanmedian(data)}")
    print()
    data = np.hstack([frc[col == i] for i in [42, 43]])
    print(f"Overall mean ice, T>260: {np.nanmean(data)}")
    print(f"Overall median ice, T>260: {np.nanmedian(data)}")

# Distributions for all temperatures/heights for liquid and ice
frc = lwd_force_mpl
col = col_heights


i40 = np.where(col >= 40)[0]
dice = frc[i40]
dliq = frc[np.setxor1d(np.where(col >= 30)[0], i40)]

bins = np.arange(0, np.nanmax(frc), 10)

fig, [ax1, ax2] = plt.subplots(2, 1, num=20, clear=True)  # , figsize=[8, 8]
ax1.set_position([0.14, 0.57, 0.83, 0.4])
ax2.set_position([0.14, 0.1, 0.83, 0.4])

ax1.hist(dice, bins)
ax1.legend(legl)
ax1.set_ylabel("Number")

ax2.hist(dliq, bins)
ax2.legend(legi)

ax2.set_ylabel("Number")
ax2.set_xlabel("Downwelling longwave cloud forcing (W m$^{-2}$)")

ax1.text(2, 25000, "a) Liquid-containing")
ax2.text(2, 2500, "b) Ice-only")
ax1.set_xlim([0, 100])
ax2.set_xlim([0, 100])

# Get medians and means
print(f"Mean, median for ice: {np.nanmean(dice)}, {np.nanmedian(dice)}")
print(f"Mean, median for liquid: {np.nanmean(dliq)}, {np.nanmedian(dliq)}")

if plot_lwd_hists_for_alts_inc_era5:
    raise ValueError("Altitude range labels on legend are incorrect!")

    # After converting:
    # -15: Data missing
    # 5: obscured
    # 15: clear
    # 20: unidentified cloud, above 1000 km
    # 23: unidentified cloud, below 400 m
    # 25: unidentified cloud, height unknown
    # 30: liq, above 1000 km
    # 31: liq,
    # 32: liq
    # 33: liq
    # 35: liq, height unknown
    # 40: ice, height unknown

    col = np.floor(col_heights)

    # Interpolate forcing onto mpl grid
    mlwd = np.interp(mpl_days, lwd_days, lwd, left=np.nan, right=np.nan)

    histcolors = ["blue", "orange", "cyan"]
    ehistcolors = ["blue", "orange", "cyan"]

    bins = np.arange(150, 400, 25)
    leg = [
        "base > 1 km",
        "400 m < base < 1 km",
        "base < 400 m",
    ]

    leg1 = ["ice, " + x for x in leg]
    leg2 = ["liq, " + x for x in leg]

    lwd_ice = [x for i, x in enumerate(mlwd) if col[i] in [40, 41, 42, 43]]
    lwd_ice = np.array(lwd_ice)
    lwd_ice = lwd_ice[np.isnan(lwd_ice) == False]

    lwd_liq = [x for i, x in enumerate(mlwd) if col[i] in [30, 31, 32, 33]]
    lwd_liq = np.array(lwd_liq)
    lwd_liq = lwd_liq[np.isnan(lwd_liq) == False]

    fig, [ax1, ax2, ax3] = plt.subplots(3, 1, num=16, clear=True)
    ax2.set_position([0.15, 0.7, 0.75, 0.25])
    ax1.set_position([0.15, 0.4, 0.75, 0.25])
    ax3.set_position([0.15, 0.1, 0.75, 0.25])

    ice_mean = 0
    hdata = ax1.hist(lwd_ice, bins, stacked=False, color="black")
    for i in range(len(bins) - 2):
        if hdata[0][i] > 0:
            tot = np.sum(hdata[0])
            pcbar = round(hdata[0][i] / tot * 100, 1)
            ice_mean += pcbar / 100 * np.mean([bins[i], bins[i + 1]])
            lab = f"{pcbar}%"
            ax1.text(hdata[1][i], hdata[0][i] + 100, lab)
    data = [
        mlwd[col == 40],
        np.hstack([mlwd[col == 41], mlwd[col == 42]]),
        mlwd[col == 43],
    ]
    ax1.hist(data, bins, stacked=True, color=histcolors, label=leg)
    ax1.set_ylabel("Number")
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles[::-1], labels[::-1])
    ax1.set_xlim([150, 350])
    ax1.set_ylim([0, 55000])
    ax1.text(155, 10000, "b) Measured, ice")
    # ax1.text(155, 6500, f"   |Flux| = {round(ice_mean)} W/m$^2$")

    data = [
        mlwd[col == 30],
        np.hstack([mlwd[col == 31], mlwd[col == 32]]),
        mlwd[col == 33],
    ]
    liq_mean = 0
    hdata = ax2.hist(lwd_liq, bins, stacked=False, color="black")
    for i in range(len(bins) - 1):
        if round(hdata[0][i] / tot * 100, 1) >= 0.01:
            tot = np.sum(hdata[0])
            pcbar = round(hdata[0][i] / tot * 100, 1)
            lab = f"{pcbar}%"
            ax2.text(hdata[1][i], hdata[0][i] + 1000, lab)
            liq_mean += pcbar / 100 * np.mean([bins[i], bins[i + 1]])

    ax2.hist(data, bins, stacked=True, color=histcolors)
    ax2.set_ylabel("Number")
    ax2.set_ylim([0, 220000])
    ax2.text(155, 180000, "a) Measured, liquid-containing")
    # ax2.text(155, 60000, f"  |Flux| = {round(liq_mean)} W/m$^2$")
    ax2.set_xlim([150, 350])

    # Make the plot for era_cbh
    era_cbh_dtime, era_cbh = get_era5_cbh("All seasons")

    # Match up the CBH and the LWD
    dtime_match, ecbh_match, elwd_match = model_meas_matches(
        era5["date"], era5["lwd"], era_cbh_dtime, era_cbh
    )
    elwd_match = np.array(elwd_match)
    ecbh_match = np.array(ecbh_match)

    i400 = np.where(ecbh_match <= 400)[0]
    i400_1000 = np.where(ecbh_match < 1000)[0]
    i1000 = np.where(ecbh_match >= 1000)[0]
    i400_1000 = np.setxor1d(i400_1000, i400)

    data = [elwd_match[i1000], elwd_match[i400_1000], elwd_match[i400]]
    datasum = np.array(elwd_match)
    datasum = datasum[np.isnan(datasum) == False]
    hdata = ax3.hist(datasum, bins, stacked=False, color="black")
    era5_mean = 0
    for i in range(len(bins) - 1):
        if hdata[0][i] > 0:
            tot = np.sum(hdata[0])
            pcbar = round(hdata[0][i] / tot * 100, 1)
            lab = f"{pcbar}%"
            ax3.text(hdata[1][i], hdata[0][i] + 100, lab)
            era5_mean += pcbar / 100 * np.mean([bins[i], bins[i + 1]])

    ax3.hist(data, bins, stacked=True, color=ehistcolors)
    ax3.set_ylabel("Number")
    ax3.set_xlabel("Downwelling longwave Flux (W m$^{-2}$)")
    ax3.text(155, 12000, "c) ERA5, all cloud")
    # ax3.text(155, 6000, f"  |Flux| = {round(era5_mean)} W/m$^2$")
    ax3.set_xlim([150, 350])
    ax3.set_ylim([0, 15000])

    if SAVEFIGS:
        figname2 = savedir + "lwd_down_hists_with_cbh_b.png"
        plt.savefig(figname2, format="png", dpi=600)


if plot_lwforcing_hists_for_alts:
    # Label ice-converted-to-liquid as ice
    # Should be a very small fraction of liquid
    col = np.floor(col_heights)

    frc = lwd_force_mpl
    bins = np.arange(0, np.nanmax(frc), 10)
    histcolors = ["cyan", "blue", "green", "orange"]
    histcolorsflip = ["orange", "green", "blue", "cyan"]

    leg = [
        f"base < {alt[0]}m",
        f"{alt[0]}m < base < {alt[1]}m",
        f"{alt[1]}m < base < {alt[2]}m",
        f"base > {alt[2]}m",
    ]

    fig, [ax1, ax2] = plt.subplots(2, 1, num=3, figsize=[6, 6], clear=True)
    ax1.set_position([0.15, 0.55, 0.75, 0.4])
    ax2.set_position([0.15, 0.08, 0.75, 0.4])

    # fmt: off
    data = [frc[col == i] for i in [30, 31, 32, 33]]
    dataflip = [frc[col == i] for i in [33, 32, 31, 30]]
    ax1.hist(dataflip, bins, stacked=True, color=histcolorsflip)
    ax1.legend(leg)
    ax1.set_ylabel("Number")
    ax1.hist(data, bins, stacked=True, color=histcolors)

    data = [frc[col == i] for i in [40, 41, 42, 43]]
    dataflip = [frc[col == i] for i in [43, 42, 41, 40]]
    ax2.hist(dataflip, bins, stacked=True, color=histcolorsflip)
    ax2.legend(leg)
    ax2.hist(data, bins, stacked=True, color=histcolors)
    ax2.set_ylabel("Number")
    ax2.set_xlabel("Downwelling longwave cloud forcing (W/m$^{2}$)")
    # fmt: on

    ax1.set_xlim([0, 100])
    ax2.set_xlim([0, 100])

    if SAVEFIGS:
        plt.savefig(savedir + "forcing_hists2b.png", format="png")
        plt.savefig(savedir + "forcing_hists2b.eps", format="eps")


if plot_lwforcing_hists_for_alts_inc_era5:
    raise ValueError("Altitude range labels on legend are incorrect!")
    frc = lwd_force_mpl

    # Label ice-converted-to-liquid as ice
    # Should be a very small fraction of liquid
    col = np.floor(col_heights)

    era_cbh_dtime, era_cbh = get_era5_cbh("All seasons")

    # Match up the ERA5 CBH, LWD, and LWD_clr
    bins = np.arange(0, np.nanmax(frc), 10)
    dtime_match, ecbh_match, elwd_match = model_meas_matches(
        era5["date"], era5["lwd"], era_cbh_dtime, era_cbh
    )
    elwd_match = np.array(elwd_match)
    ecbh_match = np.array(ecbh_match)

    dtm, ecbh_match_clr, elwd_match_clr = model_meas_matches(
        era5["date"], era5["lwd_clr"], era_cbh_dtime, era_cbh
    )
    elwd_match_clr = np.array(elwd_match_clr)
    ecbh_match_clr = np.array(ecbh_match_clr)

    # Make sure times are same
    if sum([(x - y).total_seconds() for x, y in zip(dtime_match, dtm)]) > 0:
        raise ValueError("Differences exist!")
    if sum([x - y for x, y in zip(ecbh_match, ecbh_match_clr)]) > 0:
        raise ValueError("Differences exist!")

    # histcolors = ["blue", "orange", "green", "cyan"]
    # leg = [
    #     "base > 1000m",
    #     "500m < base < 1000m",
    #     "400m < base < 500m",
    #     "base<400m",
    # ]
    ehistcolors = ["cyan", "blue", "orange"]
    histcolors = ["cyan", "blue", "orange"]
    histcolorsflip = ["orange", "blue", "cyan"]
    leg = [
        "base > 1000m",
        "400m < base < 1000m",
        "base<400m",
    ]

    i400 = np.where(ecbh_match <= 400)[0]
    i400_1000 = np.where(ecbh_match < 1000)[0]
    i1000 = np.where(ecbh_match >= 1000)[0]
    i400_1000 = np.setxor1d(i400_1000, i400)

    era_data = [
        elwd_match[i1000] - elwd_match_clr[i1000],
        elwd_match[i400_1000] - elwd_match_clr[i400_1000],
        elwd_match[i400] - elwd_match_clr[i400],
    ]
    # datasum = np.array(elwd_match)
    # datasum = datasum[np.isnan(datasum) == False]
    # hdata = ax3.hist(datasum, bins, stacked=False, color="black")
    # era5_mean = 0
    # for i in range(len(bins) - 1):
    #     if hdata[0][i] > 0:
    #         tot = np.sum(hdata[0])
    #         pcbar = round(hdata[0][i] / tot * 100, 1)
    #         lab = f"{pcbar}%"
    #         ax3.text(hdata[1][i], hdata[0][i] + 100, lab)
    #         era5_mean += pcbar / 100 * np.mean([bins[i], bins[i + 1]])

    # fmt: off
    fig, [ax1, ax2, ax3] = plt.subplots(3, 1, num=20, clear=True) # , figsize=[8, 8]
    ax1.set_position([0.14, 0.7, 0.83, 0.25])
    ax2.set_position([0.14, 0.4, 0.83, 0.25])
    ax3.set_position([0.14, 0.1, 0.83, 0.25])
    
    leg = [
        "base < 400m",
        "400m < base < 1000m",
        "base > 1000m",
    ]
    
    dum = np.array([-.3])
    for ax in [ax1, ax2, ax3]:
        ax.hist([dum, dum, dum], stacked=True, color=histcolorsflip)
    ax2.legend(leg)
    
    data = [frc[col == 30], np.hstack([frc[col == 31], frc[col == 32]]), frc[col == 33]]
    ax1.hist(data, bins, stacked=True, color=histcolors)
    ax1.set_ylabel("Number")
    
    data = [frc[col == 40], np.hstack([frc[col == 41], frc[col == 42]]), frc[col == 43]]
    ax2.hist(data, bins, stacked=True, color=histcolors)
    ax2.set_ylabel("Number")
    
    ax3.hist(era_data, bins, stacked=True, color=ehistcolors)
    ax3.set_ylabel("Number")
    ax3.set_xlabel("Downwelling longwave cloud forcing (W/m$^{2}$)")
    
    ax1.text(2, 130000, 'a) Measured, Liquid-containing')
    ax2.text(2, 2700, 'b) Measured, Ice-only')
    ax3.text(2, 13000, 'c) ERA5, All')
    
    ax1.set_xlim([0, 100])
    ax2.set_xlim([0, 100])
    ax3.set_xlim([0, 100])

    if SAVEFIGS:
        plt.savefig(savedir + "forcing_hists_w_era5.png", format="png", dpi=600)
        plt.savefig(savedir + "forcing_hists_w_era5.eps", format="eps")
