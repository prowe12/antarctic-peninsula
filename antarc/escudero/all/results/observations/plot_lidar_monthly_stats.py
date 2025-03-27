#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 11:32:46 2023

@author: prowe
"""

# Built-in modules
from matplotlib import pyplot as plt
import numpy as np
from netCDF4 import Dataset
import datetime as dt
import calendar


from antarc.escudero.all.results.get_lidar_stats_rsrc import (
    get_files,
    get_cloud_phase_stats,
    organize_for_hist,
)
from antarc.params import MEAS_DIR
from antarc.escudero.all.results.params import SAVEFIGS, SAVEDIR

# Directories
meas_dir = MEAS_DIR + "/Escudero/mpl/"


# Params
lidar_dir = meas_dir + "reprocessing_stillwell/v3/"
datestring = "20220516"
fname = f"KGI_MPLData_{datestring}.cdf"
titlestr = "King George Island MPL Data Processed " + datestring
fmt = "KGI_MPLData_%Y%m%d.cdf"
nmin_day = 24 * 60


# Get all the data

# For convenience, keys to ncid variables
tempkey = "CloudBaseTemperature"
altkey = "CloudBaseAltitude"

min_coverage = 1  # Percent uptime requirement to plot
years = [2017, 2018, 2019, 2020, 2021, 2022, 2023]
months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

nyears = len(years)
nmonths = len(months)


# Numpy array of year by month
liq = np.zeros([nyears, nmonths])
ice = np.zeros([nyears, nmonths])
cloud = np.zeros([nyears, nmonths])
obsc = np.zeros([nyears, nmonths])
clr = np.zeros([nyears, nmonths])
nodata = np.zeros([nyears, nmonths])
ndays = np.zeros([nyears, nmonths])  # number of days when data exists
ndays_mo = np.zeros([nyears, nmonths])  # number of days in the month

npts_month = np.zeros([nyears, nmonths])
npts_poss_month = np.zeros([nyears, nmonths])

ndays_alt = np.zeros([nyears, nmonths])  # number of days cloud alts exist
ndays_temp = np.zeros([nyears, nmonths])  # number of days cloud temps exist
coverage = np.nan * np.ones([nyears, nmonths])  # Percent of month with data
min_tstep = 99 * np.ones([nyears, nmonths])  # smallest time step all month
max_tstep = np.zeros([nyears, nmonths])  # largest time step all month
alt = np.zeros([nyears, nmonths])
temp = np.zeros([nyears, nmonths])
nalt = np.zeros([nyears, nmonths])
ntemp = np.zeros([nyears, nmonths])
ncloud = np.zeros([nyears, nmonths])
nunk_alt = np.zeros([nyears, nmonths])

ave_obsc = np.zeros(nmonths)
ave_clr = np.zeros(nmonths)
ave_cloud = np.zeros(nmonths)
ave_liq = np.zeros(nmonths)
ave_ice = np.zeros(nmonths)
ave_nodata = np.zeros(nmonths)

liqalt = [[[] for i in range(len(years))] for i in range(len(months))]
liqtemp = [[[] for i in range(len(years))] for i in range(len(months))]
icealt = [[[] for i in range(len(years))] for i in range(len(months))]
icetemp = [[[] for i in range(len(years))] for i in range(len(months))]
cldalt = [[[] for i in range(len(years))] for i in range(len(months))]
cldtemp = [[[] for i in range(len(years))] for i in range(len(months))]

other = []
spacings = []
basediff = []
for imo, month in enumerate(months):
    for iyr, year in enumerate(years):

        # Get starting and ending date for this year and month
        # (beginning insclusive; end non-inclusive)
        date1 = dt.datetime(year, month, 1)
        if month < 12:
            date2 = dt.datetime(year, month + 1, 1)
        elif month == 12:
            date2 = dt.datetime(year + 1, 1, 1)
        else:
            raise ValueError(f"Month must be 1 - 12 but is {month}")

        # Get list of files for the year and month
        files_daterange = get_files(lidar_dir, fmt, date1, date2)

        # Loop over files (1 per day)
        for fname in files_daterange:
            with Dataset(lidar_dir + fname) as ncid:

                # ColumnType:
                #   float32 ColumnType(Time)
                #   Unit: Unitless
                #   Long_Name: Column Data Mask (
                #     0 = obscured column,
                #     1 = clear column,
                #     2 = unidentified cloud column,
                #     3 = liquid column,
                #     4 = ice column)

                # Times with no measurement are given as -2 in ColumnType
                # Thus we can get the number of measurements on the day
                # Also, the number of possible measurements is len(ColumnType)
                npts = sum(ncid["ColumnType"][:].data >= 0)
                npts_poss = len(ncid["ColumnType"])

                # QC of dates
                if np.any(np.isnan(ncid["DataTime"])) or np.any(
                    not np.allclose(
                        np.diff(ncid["DataTime"]), 0.03333282, atol=1e-5
                    )
                ):
                    raise ValueError("Bad date time")

                # Check for nans in the data
                if np.any(np.isnan(ncid["ColumnType"][:].data)):
                    raise ValueError("Nan found in column type!")

                # Indices to times where there is each type of cloud column
                iobsc = np.where(ncid["ColumnType"][:].data == 0)[0]
                icld = np.where(ncid["ColumnType"][:].data == 2)[0]
                iliq = np.where(ncid["ColumnType"][:].data == 3)[0]
                iice = np.where(ncid["ColumnType"][:].data == 4)[0]
                inodata = np.where(ncid["ColumnType"][:].data == -2)[0]
                icloud = np.union1d(np.union1d(iliq, iice), icld)
                ncloud_day = len(icloud) + len(iobsc)

                # Make sure the number of points totals as expected
                nclr = sum(ncid["ColumnType"][:].data == 1)
                if len(iobsc) + len(icloud) + len(inodata) + nclr != npts_poss:
                    raise ValueError("Points do not sum to 100%")

                if len(iobsc) + len(icloud) + nclr != npts:
                    raise ValueError("Points do not sum to 100%")

                cbtemp = ncid["CloudBaseTemperature"][:].data
                alt_km = ncid["CloudBaseAltitude"][:].data

                # Get cloud base altitudes and temperatures
                # for cases where clouds are low
                inan = np.where(np.isnan(ncid[tempkey]))[0]
                irep = np.intersect1d(icloud, inan)

                # Get cloud alts and temperatures for the lowest clouds
                # which weren't previously classified
                if len(irep) > 0:
                    ialt = [
                        np.where(ncid["CloudMask"][i, :] == 1)[0][0]
                        for i in irep
                    ]
                    temper = ncid["Temperature"][:].data
                    cbtemp[irep] = [
                        temper[irep[i], ialt[i]] for i in range(len(irep))
                    ]
                    alt_km[irep] = ncid["Range"][ialt]

                # Indices to times where there is cloud base information
                ialt = np.where(np.isnan(alt_km) == False)[0]
                itemp = np.where(np.isnan(cbtemp) == False)[0]

                # Indices to times where there is base info for each type
                ialt_liq = np.intersect1d(iliq, ialt)
                ialt_ice = np.intersect1d(iice, ialt)
                ialt_cld = np.intersect1d(np.union1d(iobsc, icld), ialt)
                itemp_liq = np.intersect1d(iliq, itemp)
                itemp_ice = np.intersect1d(iice, itemp)
                itemp_cld = np.intersect1d(np.union1d(iobsc, icld), itemp)

                if (
                    len(ialt_liq) != len(iliq)
                    or len(ialt_ice) != len(iice)
                    or len(itemp_liq) != len(iliq)
                    or len(itemp_ice) != len(iice)
                ):
                    raise ValueError("Missing cloud heights?")
                alt_km = alt_km / 1000.0

                # Lists of cloud base alts and temps for each type
                liqalt[imo][iyr] += list(alt_km[ialt_liq])
                icealt[imo][iyr] += list(alt_km[ialt_ice])
                cldalt[imo][iyr] += list(alt_km[ialt_cld])
                liqtemp[imo][iyr] += list(cbtemp[itemp_liq].data)
                icetemp[imo][iyr] += list(cbtemp[itemp_ice].data)
                cldtemp[imo][iyr] += list(cbtemp[itemp_cld].data)
                ncloud[iyr, imo] += ncloud_day

                # Count this day for the monthly averages
                npts_month[iyr, imo] += npts
                npts_poss_month[iyr, imo] += npts_poss
                ndays[iyr, imo] += 1

                # Get the coverage stats for this day and get the overall
                # max and min for the month
                min_tstep0 = min(np.diff(ncid["DataTime"]))  # hour
                max_tstep0 = max(np.diff(ncid["DataTime"]))  # hour
                min_tstep[iyr, imo] = min(min_tstep[iyr, imo], min_tstep0)
                max_tstep[iyr, imo] = max(max_tstep[iyr, imo], max_tstep0)

                # Get the fraction of times with data for each mask type
                # and the fraction of all time for which there is no data
                nfobsc = len(iobsc) / npts
                nfclr = nclr / npts
                nfcld = len(icld) / npts
                nfliq = len(iliq) / npts
                nfice = len(iice) / npts
                n_nodata = len(inodata) / npts_poss

                # Add this fraction to monthly sum. Thus will need to divide
                # by days per month later
                obsc[iyr, imo] += nfobsc
                clr[iyr, imo] += nfclr
                cloud[iyr, imo] += nfcld
                liq[iyr, imo] += nfliq
                ice[iyr, imo] += nfice
                nodata[iyr, imo] += n_nodata

                # Get the average cloud base and temperature for the day
                # add them to monthly sum (will average over days/month later)
                # Some of the altitudes and temperatures are nan, so we
                # need to get the average for times when they were real
                # First get the number of real alts and temps
                # If there are any clouds on this day, increment the days
                # for cloud alts and temps and the sums
                nalt = len(alt_km)  # sum(~np.isnan(ncid["CloudBaseAltitude"]))
                if nalt > 0:
                    ndays_alt[iyr, imo] += 1
                    alt[iyr, imo] += np.nansum(alt_km) / nalt
                    # np.nansum(ncid["CloudBaseAltitude"]) / nalt

                ntemp = len(cbtemp)  # sum(~np.isnan(ncid[tempkey]))
                if ntemp > 0:
                    ndays_temp[iyr, imo] += 1
                    temp[iyr, imo] += np.nansum(cbtemp) / ntemp
                    # temp[iyr, imo] += np.nansum(ncid[tempkey]) / ntemp

        # Get coverage for the month as number of days with measurements
        # divided by total number of days
        dayrange = calendar.monthrange(year, month)
        ndays_mo[iyr, imo] = dayrange[1]
        coverage[iyr, imo] = 100 * ndays[iyr, imo] / ndays_mo[iyr, imo]

        # Get monthly percentages for each type by dividing by the number of
        # days measurements were made in the month.
        # For cloud alts and temps, this may be fewer
        # than the number of days with data, so use corresponding ndays
        # Todo: Note that ndays_alt and ndays_temp should be same but deal
        # with that later. Also could have an entire month without clouds
        # if there were very few days of measurements in the month
        if coverage[iyr, imo] > min_coverage:

            # Get the sums over all years
            ave_obsc[imo] += obsc[iyr, imo]
            ave_clr[imo] += clr[iyr, imo]
            ave_cloud[imo] += cloud[iyr, imo]
            ave_liq[imo] += liq[iyr, imo]
            ave_ice[imo] += ice[iyr, imo]
            ave_nodata[imo] += nodata[iyr, imo]

            # obsc[iyr, imo] = 100 * obsc[iyr, imo] / ndays[iyr, imo]
            # clr[iyr, imo] = 100 * clr[iyr, imo] / ndays[iyr, imo]
            # cloud[iyr, imo] = 100 * cloud[iyr, imo] / ndays[iyr, imo]
            # liq[iyr, imo] = 100 * liq[iyr, imo] / ndays[iyr, imo]
            # ice[iyr, imo] = 100 * ice[iyr, imo] / ndays[iyr, imo]
            # nodata[iyr, imo] = 100 * nodata[iyr, imo] / ndays[iyr, imo]
            alt[iyr, imo] = alt[iyr, imo] / ndays_alt[iyr, imo]
            temp[iyr, imo] = temp[iyr, imo] / ndays_temp[iyr, imo]

        else:
            obsc[iyr, imo] = np.nan
            clr[iyr, imo] = np.nan
            cloud[iyr, imo] = np.nan
            liq[iyr, imo] = np.nan
            ice[iyr, imo] = np.nan
            nodata[iyr, imo] = np.nan
            alt[iyr, imo] = np.nan
            temp[iyr, imo] = np.nan


# convert sums to aves
tot = ave_obsc + ave_clr + ave_cloud + ave_liq + ave_ice
ave_obsc = ave_obsc / tot
ave_clr = ave_clr / tot
ave_cloud = ave_cloud / tot
ave_liq = ave_liq / tot
ave_ice = ave_ice / tot

# Print some stats
# Fraction of the time there liquid, of time when phase is known
fliq_of_known = ave_liq / (ave_liq + ave_ice)
fice_of_known = ave_ice / (ave_liq + ave_ice)
# Corrected fration of lidar-operating time when it is liquid
liq_time = ave_liq + fliq_of_known * (ave_obsc + ave_cloud)
ice_time = ave_ice + fice_of_known * (ave_obsc + ave_cloud)

print("Stats for paper")
print(f"Liq/ice cloud that is liquid: {100*np.mean(fliq_of_known)}%")
print(f"Liq/ice cloud that is ice: {100*np.mean(fice_of_known)}%")

print(f"Fraction of lidar op time when clear: {100*np.mean(ave_clr)}%")
print(f"Fraction of lidar op time when liq-containg: {100*np.mean(liq_time)}%")
print(f"Fraction of lidar op time when ice-only: {100*np.mean(ice_time)}%")

clr_seas = []
liq_seas = []
ice_seas = []
for inds in [[11, 0, 1], [2, 3, 4], [5, 6, 7], [8, 9, 10]]:
    clr_seas += [np.mean(ave_clr[inds])]
    liq_seas += [np.mean(liq_time[inds])]
    ice_seas += [np.mean(ice_time[inds])]
print(
    "Fraction of time when clear, by season (summer, fall, winter, spring): "
)
print(f"{[100*x for x in clr_seas]}%")
print("Fraction of time when liq-containg by seas: ")
print(f"{[100*x for x in liq_seas]}%")
print("Fraction of time when ice-only, by season: ")
print(f"{[100*x for x in ice_seas]}%")

coverage_pts = npts_month / (ndays_mo * 720) * 100

# QC: make sure the totals sum to 100% (or nan) as expected
tots = (obsc + clr + cloud + liq + ice) / ndays
if not np.allclose(tots[np.isnan(tots) == False], 1):
    raise ValueError("Cloud type percentages must sum to 100%")


[totice, totcld] = get_cloud_phase_stats(icetemp, ncloud)
[totliq, totcld] = get_cloud_phase_stats(liqtemp, ncloud)
[totunkcld, totcld] = get_cloud_phase_stats(cldtemp, ncloud)

print("Total number of ice, liquid-containing, unknown, and all cloud")
print(totice)
print(totliq)
print(totunkcld)
print(totcld)
print()

# print("\nPercentages")
# print(f"Ice-only: {totice / totcld * 100}")
# print(f"Liq&ice:  {totliq / totcld * 100}")
# print(f"Unknown:  {totunkcld / totcld * 100}")
# print(f"total:    {(totice + totliq + totunkcld) / totcld * 100}\n")


# # # # #     Plot cloud mask as percentage totaling up to 100%      # # # # #
# Get the percentage of time each cloud mask applies for year and month
fclr = 100 * clr / ndays
fice = 100 * ice / ndays
fliq = 100 * liq / ndays
fobsc = 100 * obsc / ndays
funk = 100 * cloud / ndays

# Get the stacked means over years
sclr = np.nanmean(fclr, axis=0)
sice = np.nanmean(fclr + fice, axis=0)
sliq = np.nanmean(fclr + fice + fliq, axis=0)
sobsc = np.nanmean(fclr + fice + fliq + fobsc, axis=0)
sunk = np.nanmean(fclr + fice + fliq + fobsc + funk, axis=0)

sclr = 100 * ave_clr
sice = 100 * (ave_clr + ave_ice)
sliq = 100 * (ave_clr + ave_ice + ave_liq)
sobsc = 100 * (ave_clr + ave_ice + ave_liq + ave_obsc)
sunk = 100 * (ave_clr + ave_ice + ave_liq + ave_obsc + ave_cloud)

mon = [m for m in months]
mon[0] = 0.9
mon[-1] = 12.1
unk_cloud = obsc + cloud

# Make stacked histogram plot of lidar classifications
fig, ax = plt.subplots(1, 1, num=1, clear=True)
ax.set_position([0.1, 0.1, 0.75, 0.85])

ax.plot(months, clr.T, "+")
ax.set_prop_cycle(None)
ax.fill_between(mon, sunk, color=[0.8, 0.8, 0.8])
ax.fill_between(mon, sobsc, color=[0.6, 0.6, 0.6])
ax.fill_between(mon, sliq, color=[0.6, 0.9, 1])
ax.fill_between(mon, sice, color="white")
ax.fill_between(mon, sclr, color="cyan")

ax.set_prop_cycle(None)
ax.plot(months, fclr.T, "+")
ax.set_prop_cycle(None)
ax.plot(months, (fclr + fice).T, "x")
ax.set_prop_cycle(None)
ax.plot(months, (fclr + fice + fliq).T, "o", fillstyle="none")

ax.text(7.8, 60, "Liquid-Containing")
ax.text(8.1, 15, "Ice")
ax.text(10.3, 94, "Obscured")
ax.text(9.6, 102, "Unknown Cloud")
ax.text(8.1, 3, "Clear")
ax.set_ylabel("Cloud Mask (%)")
ax.set_xlabel("Month")
ax.set_xlim([0.9, 12.1])
ax.set_ylim([0, 100])
ax.legend(years, loc="upper left", bbox_to_anchor=(1.0, 0.92))

if SAVEFIGS:
    plt.savefig(SAVEDIR + "lidar_class_by_month_stack.png", format="png")
    plt.savefig(SAVEDIR + "lidar_class_by_month_stack.eps", format="eps")


# When was the covereage greater than 80%?
fclr_80 = np.nan + np.zeros(np.shape(coverage))
for iyr in range(len(coverage)):
    iuse = np.where(coverage[iyr] >= 80)[0]
    fclr_80[iyr, iuse] = fclr[iyr, iuse]
np.nanstd(fclr, axis=0)
np.nanstd(fclr_80, axis=0)


# Seasonal:
tot_liq = ave_liq + ave_obsc
print(f"Summer liquid fraction: {(np.sum(tot_liq[:2]) + tot_liq[12-1])/3}")
print(f"Fall liquid fraction: {np.mean(tot_liq[2:5])}")
print(f"Winter liquid fraction: {np.mean(tot_liq[5:8])}")
print(f"Spring liquid fraction: {np.mean(tot_liq[8:11])}\n")

print(f"Summer ice fraction: {(np.sum(ave_ice[:2]) + ave_ice[12-1])/3}")
print(f"Fall ice fraction: {np.mean(ave_ice[2:5])}")
print(f"Winter ice fraction: {np.mean(ave_ice[5:8])}")
print(f"Spring ice fraction: {np.mean(ave_ice[8:11])}\n")

print(f"Clear fraction summer: {(np.sum(ave_clr[:2]) + ave_clr[12-1])/3}")
print(f"Clear fraction fall: {np.mean(ave_clr[2:5])}")
print(f"Clear fraction winter: {np.mean(ave_clr[5:8])}")
print(f"Clear fraction sprint: {np.mean(ave_clr[8:11])}\n")

print(f"Cloudy fraction summer: {1-(np.sum(ave_clr[:2]) + ave_clr[12-1])/3}")
print(f"Cloudy fraction fall: {1-np.mean(ave_clr[2:5])}")
print(f"Cloudy fraction winter: {1-np.mean(ave_clr[5:8])}")
print(f"Cloudy fraction spring: {1-np.mean(ave_clr[8:11])}\n")

# # # # # #     Plot instrument coverage    # # # # # # # # # # # # # # # #
coverage[coverage == 0] = np.nan
leg = [str(x) for x in years] + ["cutoff"]
symbols = ["s", "o", "^", "+", ".", "x", "v"]

fig, ax = plt.subplots(figsize=[5, 2.6], num=2, clear=True)
ax.set_position([0.13, 0.2, 0.84, 0.78])

for i in range(3):
    ax.plot(
        months,
        coverage[i, :],
        marker=symbols[i],
        linestyle="none",
        fillstyle="none",
    )
for i in range(3, len(symbols)):
    ax.plot(months, coverage[i, :], marker=symbols[i], linestyle="none")
ax.plot([0, 13], [min_coverage, min_coverage], "k:", label="cutoff")
ax.set_ylabel("Instrument uptime (%)")
ax.set_ylim([0, 102])
ax.set_xlim([0.9, 12.1])
ax.set_xlabel("Month")
ax.legend(leg, loc="center", bbox_to_anchor=(0.76, 0.45))
ax.set_ylim([0, 105])
if SAVEFIGS:
    plt.savefig(SAVEDIR + "lidar_uptime.png", format="png")
    plt.savefig(SAVEDIR + "lidar_uptime.eps", format="eps")


# # # # # #     Print averages for seasons    # # # # # # # # # # # # # # # #
cloudiness = 100 - np.nanmean(fclr, axis=0)
cloudpc_djf = (cloudiness[11] + cloudiness[0] + cloudiness[1]) / 3
cloudpc_mam = (cloudiness[2] + cloudiness[3] + cloudiness[4]) / 3
cloudpc_jja = (cloudiness[5] + cloudiness[6] + cloudiness[7]) / 3
cloudpc_son = (cloudiness[8] + cloudiness[9] + cloudiness[10]) / 3

print(f"Overall cloudiness (%) {np.mean(cloudiness)}")
print(f"DJF cloudiness (%) {cloudpc_djf}")
print(f"MAM cloudiness (%) {cloudpc_mam}")
print(f"JJA cloudiness (%) {cloudpc_jja}")
print(f"SON cloudiness (%) {cloudpc_son}")


# # # # # #    Plot probability density with altitude by season   # # # # # #
liq_alt = organize_for_hist(liqalt)
ice_alt = organize_for_hist(icealt)
liq_temp = organize_for_hist(liqtemp)
ice_temp = organize_for_hist(icetemp)


bot = 0.15
wid = 0.39
height = 0.8
col = ["tab:red", [1, 0.6, 0], "blue", [1, 0.8, 0]]


liq_alt_tot = 0
ice_alt_tot = 0
for lalts, ialts in zip(liq_alt, ice_alt):
    liq_alt_tot += len(lalts)
    ice_alt_tot += len(ialts)
tot_alts = liq_alt_tot + ice_alt_tot

bins = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
flip = "horizontal"
season = ["DJF", "MAM", "JJA", "SON"]
season_liq = []
season_ice = []
for i, alts in enumerate(liq_alt):
    lpc = round(100 * len(alts) / tot_alts)
    ipc = round(100 * len(ice_alt[i]) / tot_alts)
    season_liq.append(f"{season[i]} ({lpc}%)")
    season_ice.append(f"{season[i]} ({ipc}%)")

plt.figure(num=3, clear=True, figsize=(4.5, 3))
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
ax1.set_position([0.115, bot, wid, height])
ax2.set_position([0.575, bot, wid, height])

liq_hist = ax1.hist(
    liq_alt,
    orientation=flip,
    density=True,
    bins=bins,
    label=season_liq,
    color=col,
)
ax1.set_ylim([0, 10])
ax1.legend(loc="upper left", bbox_to_anchor=[0.06, 0.9])
ax1.set_xlabel("Probability Density")
ax1.set_ylabel("Altitude (km)")
ax1.text(0.05, 9.3, "a) Liquid-containing")

ice_hist = ax2.hist(
    ice_alt,
    orientation=flip,
    density=True,
    bins=bins,
    label=season_ice,
    color=col,
)
ax2.set_ylim([0, 10])
ax2.legend(bbox_to_anchor=[0.06, 0.9])
ax2.set_xlabel("Probability Density")
ax2.text(0.05, 9.3, "b) Ice-only")

if SAVEFIGS:
    plt.savefig(
        SAVEDIR + "fig4_cloud_base_height_hist.png", format="png", dpi=300
    )
    plt.savefig(SAVEDIR + "fig4_cloud_base_height_hist.eps", format="eps")


lcb = []
for x in liq_alt:
    lcb += x
lcb = np.array(lcb)
print(f"Percentage of liq clouds > 4km: {len(lcb[lcb>4])/len(lcb)*100}")

lcb = []
for x in ice_alt:
    lcb += x
lcb = np.array(lcb)
print(f"Percentage of ice clouds > 4km: {len(lcb[lcb>4])/len(lcb)*100}")


# Altitude stats
lcb = []
for x in liq_alt:
    lcb += x
for x in ice_alt:
    lcb += x
lcb = np.array(lcb)

print(f"Percentage of clouds <= 500 m: {len(lcb[lcb<=.5])/len(lcb)*100}")
print(f"Percentage of clouds <= 1 km: {len(lcb[lcb<=1])/len(lcb)*100}")
print(f"Percentage of clouds <= 3 km: {len(lcb[lcb<=3])/len(lcb)*100}")

# Assuming 9.7% of cloudy cases are obscured, 90.3% not obscured
fobsc = 9.7
fcbh = 100 - 9.7
for x in [0.5, 1, 3]:
    fract = len(lcb[lcb <= x]) / len(lcb) * fcbh
    low = fract
    high = fract + fobsc
    print(f"500 m uncertainty range: {low} - {high}")


lcb = []
for x in liq_alt:
    lcb += x
lcb = np.array(lcb)
print(f"Percentage of liquid clouds > 3 km: {len(lcb[lcb>3])/len(lcb)*100}")

lcb = []
for x in ice_alt:
    lcb += x
lcb = np.array(lcb)
print(f"Percentage of ice clouds > 3 km: {len(lcb[lcb>3])/len(lcb)*100}")


# # # # # #    Plot cloud base temperature # # # # # # # # # # # # # #
nice = 0
nliq = 0
liq_tempp = [[] for i in range(4)]
ice_tempp = [[] for i in range(4)]
for i in range(4):
    liq_tempp[i] = [x for x in liq_temp[i]]
    liq_tempp[i] += [x for x in ice_temp[i] if x > 276.15]
    ice_tempp[i] = [x for x in ice_temp[i] if x <= 276.15]
    nice += len(ice_tempp[i])
    nliq += len(liq_tempp[i])

pc_baseid = (nice + nliq) / np.sum(np.sum(ncloud)) * 100
print(f"Percentage of clouds with bases identified: {round(pc_baseid)}%")

bins = np.array([200, 213, 223, 233, 243, 253, 263, 273, 280]) + 3.15
left = 0.11
wid = 0.87
height = 0.38
col = ["tab:red", [1, 0.6, 0], "blue", [1, 0.8, 0]]

plt.figure(num=4, clear=True)
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
ax1.set_position([left, 0.58, wid, height])
ax2.set_position([left, 0.1, wid, height])

ax1.hist(liq_tempp, density=True, bins=bins, label=season, color=col)
for x in bins:
    ax1.plot([x, x], [0, 1], ":", color=[0.8, 0.8, 0.8])
ax1.set_ylim([0, 0.09])
ax1.legend()
ax1.set_ylabel("Probability Density")
ax1.text(201, 0.008, f"a) Liquid-containing ({nliq} observations)")

ax2.hist(ice_tempp, density=True, bins=bins, label=season, color=col)
for x in bins:
    ax2.plot([x, x], [0, 1], ":", color=[0.8, 0.8, 0.8])
ax2.set_ylim([0, 0.09])
ax2.legend()
ax2.set_ylabel("Probability Density")
ax2.set_xlabel("Temperature (K)")
ax2.text(201, 0.008, f"b) Ice-only ({nice} observations)")

if SAVEFIGS:
    plt.savefig(
        SAVEDIR + "cloud_base_temperature_hist.png", format="png", dpi=600
    )
    plt.savefig(SAVEDIR + "cloud_base_temperature_hist.eps", format="eps")


# Get fraction supercooled
bins = np.array([200, 213, 223, 233, 243, 253, 263, 273, 280, 290]) + 0.15
liq_temp_all = [x for xs in liq_temp for x in xs]
ice_temp_all = [x for xs in ice_temp for x in xs]
temp_all = liq_temp_all + ice_temp_all

plt.figure()
temp_hist = plt.hist(liq_temp_all, bins=bins)
fract = np.sum(temp_hist[0][:7]) / np.sum(temp_hist[0])
print(f"fraction of liquid cloud bases above 0C: {fract}")

# zoom in on warm clouds
bins = np.arange(273, 280) + 0.15
left = 0.11
wid = 0.87
height = 0.38
col = ["tab:red", [1, 0.6, 0], "blue", [1, 0.8, 0]]

plt.figure(num=5, clear=True)
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
ax1.set_position([left, 0.58, wid, height])
ax2.set_position([left, 0.1, wid, height])

ax1.hist(liq_temp, density=False, bins=bins, label=season, color=col)
# ax1.set_ylim([0, 10])
ax1.legend()
ax1.set_ylabel("Probability Density")
ax1.text(277, 5000, f"a) Liquid-containing (of {nliq} cases)")

ax2.hist(ice_temp, density=False, bins=bins, label=season, color=col)
# ax2.set_ylim([0, 10])
ax2.legend()
ax2.set_ylabel("Probability Density")
ax2.set_xlabel("Temperature (K)")
ax2.text(277, 500, f"b) Ice-only (of {nice} cases)")


# # # # # #    Plot mean cloud base temperature and altitude    # # # # # #
fix, axs = plt.subplots(2, 1, figsize=[8, 8], num=6, clear=True)
axs[0].set_position([0.1, 0.55, 0.75, 0.4])
axs[1].set_position([0.1, 0.08, 0.75, 0.4])

axs[0].plot(months, np.nanmean(alt, 0))
axs[0].plot(months, alt.T, "o")
axs[0].set_ylabel("Cloud base altitude (m)")
axs[0].legend(["mean"] + [str(x) for x in years], bbox_to_anchor=(1.01, 1.05))

axs[1].plot(months, np.nanmean(temp, 0))
axs[1].plot(months, temp.T, "o")
axs[1].set_xlabel("Month")
axs[1].set_ylabel("Cloud base temperature (K)")

if SAVEFIGS:
    plt.savefig(SAVEDIR + "cloud_base.png", format="png")
    plt.savefig(SAVEDIR + "cloud_base.eps", format="eps")
