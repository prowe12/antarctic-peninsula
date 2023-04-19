#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 08:56:46 2022

@author: prowe

Purpose: Plot broadband radiation at Marambio on 2022/02
"""

# # Dependencies
# from matplotlib import pyplot as plt
# import numpy as np
# import datetime as dt
# import matplotlib.dates as mdates

# # My modules
# from antarc.marambio.case202202.get_swd_flux import get_swd
# from antarc.marambio.case202202.get_lwd_flux import get_lwd, get_lwd_clear
# from antarc.marambio.parameters import rad_flux_params
# from antarc.marambio.case202202.load_rad_flux import load_rad_flux

# # from antarc.marambio.case202202.get_rad_flux import get_rad_flux

# bbrad_dir = rad_flux_params.RAD_FLUX_DIR
# bbrad_fmt = rad_flux_params.RAD_FLUX_FILEFORMAT

# out_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/Feb_2022/figures/"
# fig1name = out_dir + "mbio_lwd_swd.png"
# fig2name = out_dir + "mbio_lwd_swd_zoom.png"

# utc = dt.timezone.utc
# savefigs = False

# # A longer date range
# date1 = dt.datetime(2022, 2, 1, tzinfo=utc)
# date2 = dt.datetime(2022, 2, 15, tzinfo=utc)

# # Get the radiative fluxes
# dtimes, swdir, sw_spn1 = load_rad_flux(bbrad_dir, bbrad_fmt, date1, date2)
# # start_dt, end_dt)


# # SWD
# swd_date, swd, swd_clear_lib, swd_clear, diffuse, diffuse_fac = get_swd()
# pyr_dir += str(2022) + "/"
# swd_date_long, swd_long = load_rad_flux(pyr_dir, pyr_fmt, date1, date2)

# # LWD
# lwd_clear_date, lwd_clear = get_lwd_clear()
# lwd_date, lwd = get_lwd(date1, date2)
# lwd_tstamp = [x.timestamp() for x in lwd_date]
# lwd_clear_tstamp = [x.timestamp() for x in lwd_clear_date]
# lwd_clear = np.array(lwd_clear)

# # Get LWD forcing
# lwd_clear_interp = np.interp(lwd_tstamp, lwd_clear_tstamp, lwd_clear)

# # Index to where data is missing
# isw = np.where(swd_date < dt.datetime(2022, 2, 9, tzinfo=utc))[0][-1] + 1
# ilw = np.where(lwd_date < dt.datetime(2022, 2, 9, tzinfo=utc))[0][-1] + 1
# iswl = np.where(swd_date_long < dt.datetime(2022, 2, 9, tzinfo=utc))[0][-1] + 1

# # Cumulative forcing for LWD
# start_date = dt.datetime(2022, 2, 6, tzinfo=utc)  # swd_date[0]
# final_date = swd_date[isw - 1]
# # dt.datetime(fdate.year, fdate.month, fdate.day, fdate.hour + 1, 0)
# i1 = np.where(np.array(lwd_tstamp) >= start_date.timestamp())[0][0]
# i2 = np.where(np.array(lwd_tstamp) <= final_date.timestamp())[0][-1] + 1
# lwd_force = lwd - lwd_clear_interp
# ntimes = i2 - i1
# lwf_tot = np.zeros([ntimes])
# lwf_tot_mwh = np.zeros([ntimes])
# count = 0
# for i in range(i1, i2):
#     lwf = lwd_force[i] * (lwd_tstamp[i] - lwd_tstamp[i - 1])
#     lwf_mwh = lwd_force[i] * (lwd_tstamp[i] - lwd_tstamp[i - 1]) / 3600 / 1000
#     lwf_tot[count] = lwf_tot[count - 1] + lwf
#     lwf_tot_mwh[count] = lwf_tot_mwh[count - 1] + lwf_mwh
#     count += 1
# il1 = i1
# il2 = i2

# # Cumulative forcing for SWD
# swd_tstamp = [x.timestamp() for x in swd_date]
# i1 = np.where(np.array(swd_tstamp) >= start_date.timestamp())[0][0]
# i2 = np.where(np.array(swd_tstamp) <= final_date.timestamp())[0][-1]
# swd_force = swd - swd_clear
# ntimes = i2 - i1
# swf_tot = np.zeros([ntimes])
# swf_tot_mwh = np.zeros([ntimes])
# count = 0
# for i in range(i1, i2):
#     swf = swd_force[i] * (swd_tstamp[i] - swd_tstamp[i - 1])
#     swf_mwh = swd_force[i] * (swd_tstamp[i] - swd_tstamp[i - 1]) / 3600 / 1000
#     swf_tot[count] = swf_tot[count - 1] + swf
#     swf_tot_mwh[count] = swf_tot_mwh[count - 1] + swf_mwh
#     count += 1
# is1 = i1
# is2 = i2


# # lwd_tot = np.interp(lwd_tstamp, lwd_clear_interp)
# # swd_tot = np.interp(lwd_tstamp, lwd_clear_interp)
# lwf_tot_sw_date = np.interp(swd_tstamp[is1:is2], lwd_tstamp[il1:il2], lwf_tot)
# lwd_force_sw_date = np.interp(swd_tstamp, lwd_tstamp, lwd_force)

# rot = 0
# datefmt = "%d"
# plt.figure(num=1)
# ax1 = plt.subplot(211)
# ax3 = plt.subplot(212, sharex=ax1)
# wid = 0.85
# hit = 0.4
# ax1.set_position([0.12, 0.57, wid, hit])
# ax3.set_position([0.12, 0.10, wid, hit])

# ax1.plot(swd_date_long[:iswl], swd_long[:iswl], color="red", label="Measured")
# ax1.plot(swd_date_long[iswl:], swd_long[iswl:], color="red")
# ax1.plot(swd_date[:isw], swd_clear[:isw], "--", color="red", label="Clear")
# ax1.plot(swd_date[isw:], swd_clear[isw:], "--", color="red")
# ind = range(0, len(swd_date), 50)  # , 200)
# ax1.plot(swd_date[ind], swd_clear_lib[:, 1], "k--")  # global down
# # ax1.plot(swd_date[ind], swd_clear_lib[:, 0], "ks")  # direct
# # ax1.plot(swd_date[ind], swd_clear_lib[:, 2], "c.")  # diffuse down
# # ax1.plot(swd_date[ind], swd_clear_lib[:, 3], "m+")  # diffuse up
# ax1.legend()
# ax1.set_xlim(
#     dt.datetime(2022, 2, 5, tzinfo=utc), dt.datetime(2022, 2, 11, tzinfo=utc)
# )
# ax1.set_ylim([0, 1000])
# ax1.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# ax1.tick_params(axis="x", rotation=rot)
# ax1.set_ylabel("SWD (W/m$^{2}$)")
# # ax1.legend(loc=(0.1,.3))

# ax3.plot(lwd_date[:ilw], lwd[:ilw], color="blue", label="Measured")
# ax3.plot(lwd_date[ilw:], lwd[ilw:], color="blue")
# ax3.plot(lwd_clear_date, lwd_clear, "o:", color="blue", label="Clear")
# ax3.legend()
# ax3.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# ax3.tick_params(axis="x", rotation=rot)
# ax3.set_xlabel("Day of 2022/02")
# ax3.set_ylabel("LWD (W/m$^{2}$)")
# ax3.set_ylim([0, 1000])
# ax3.legend(loc="upper left")

# if savefigs:
#     plt.savefig(fig1name)


# ax1 = plt.subplot(221)
# ax2 = plt.subplot(222, sharex=ax1)
# ax3 = plt.subplot(223, sharex=ax1)
# ax4 = plt.subplot(224, sharex=ax1)

# ax2.plot(lwd_date[:ilw], lwd_force[:ilw], color="blue", label="LWD Forcing")
# ax2.plot(swd_date[:isw], swd_force[:isw], color="red", label="SWD Forcing")
# ax2.legend()
# ax2.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# ax2.tick_params(axis="x", rotation=rot)

# ax4.plot(lwd_date[:ilw], lwf_tot[:ilw] / 1e6, color="blue", label="LWD")
# ax4.plot(swd_date[:isw], swf_tot[:isw] / 1e6, color="red", label="SWD")
# ax4.plot(
#     swd_date[:isw],
#     (swf_tot[:isw] + lwf_tot_sw_date[:isw]) / 1e6,
#     "k",
#     label="Total",
# )
# ax4.legend()
# ax4.set_ylabel("Cumulative Forcing (MJ/m$^{2}$)")
# ax4.xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
# ax4.set_xlabel("Day of 2022/02")
# ax4.tick_params(axis="x", rotation=rot)

# totf = swf_tot + lwf_tot_sw_date
# totfstr = str(round(totf[-1] / 1e6))

# rot = 40
# datefmt = "%d %H"
# fdate = swd_date[isw - 1]
# # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
# plt.close(2)
# plt.figure(num=2)
# ax1 = plt.subplot(221)
# ax2 = plt.subplot(222, sharex=ax1)
# ax3 = plt.subplot(223, sharex=ax1)
# ax4 = plt.subplot(224, sharex=ax1)
# wid = 0.37
# hit = 0.36
# ax1.set_position([0.1, 0.62, wid, hit])
# ax2.set_position([0.6, 0.62, wid, hit])
# ax3.set_position([0.1, 0.15, wid, hit])
# ax4.set_position([0.6, 0.15, wid, hit])

# ax1.plot(swd_date[:isw], swd[:isw], color="red", label="Meas")
# ax1.plot(swd_date[:isw], swd_clear[:isw], "--", color="red", label="Clear")
# ax1.legend()
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
# ax1.tick_params(axis="x", rotation=rot)
# ax1.set_xlim(start_date, final_date)
# ax1.set_ylabel("SWD (W/m$^{2}$)")
# ax1.legend(loc=[0.4, 0.7])

# ax3.plot(lwd_date[:ilw], lwd[:ilw], color="blue", label="Measured")
# ax3.plot(lwd_clear_date, lwd_clear, "o:", color="blue", label="Clear")
# ax3.legend()
# ax3.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
# ax3.tick_params(axis="x", rotation=rot)
# ax3.set_xlabel("Day and Hour in 2022")
# ax3.set_ylabel("LWD (W/m$^{2}$)")

# tot_force = swd_force[:isw] + lwd_force_sw_date[:isw]
# ax2.plot(lwd_date[:ilw], lwd_force[:ilw], color="blue", label="LWD")
# ax2.plot(swd_date[:isw], swd_force[:isw], color="red", label="SWD")
# ax2.plot(swd_date[:isw], tot_force, "k", label="Total")
# ax2.plot(lwd_date[:ilw], np.zeros(ilw), "k:")
# ax2.legend()
# ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d %H"))
# ax2.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
# ax2.tick_params(axis="x", rotation=rot)
# ax2.set_ylabel("Forcing (W/m$^{2}$)")
# ax2.legend(loc=[0.4, 0.02])

# ax4.plot(lwd_date[il1:il2], lwf_tot / 1e6, color="blue", label="LWD")
# ax4.plot(swd_date[is1:is2], swf_tot / 1e6, color="red", label="SWD")
# ax4.plot(swd_date[is1:is2], totf / 1e6, "k", label="Total")
# ax4.plot(swd_date[is2], totf[-1] / 1e6, "ko")
# ax4.plot(lwd_date[il1:il2], np.zeros(il2 - il1), "k:")
# ax4.legend()
# ax4.set_ylabel("Cumulative Forcing (MJ/m$^{2}$)")
# ax4.xaxis.set_major_formatter(mdates.DateFormatter("%d %H"))
# ax4.tick_params(axis="x", rotation=rot)
# ax4.set_xlabel("Day and Hour in 2022")
# ax4.text(
#     dt.datetime(2022, 2, 7, 15, 30, tzinfo=utc), -15, totfstr, fontsize=14
# )

# if savefigs:
#     plt.savefig(fig2name)
