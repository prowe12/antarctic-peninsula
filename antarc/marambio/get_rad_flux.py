#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 16:22:16 2019

@author: prowe

From Claudia:
    Quality controlled data from Marambio:
    singles files because for me it is easier to use them this way, 
    but if you need I can divide it. Not knowing which format you prefer 
    I also uploaded them as csv and plain text (.dat).
    ERL stands for the data within the Extremely Rare Limits, whereas 
    PPL for those within the Physically Possible Limits. 
    Such limits are from the recommended controls by the BSRN 
    (see the attachment), but these data are not checked for the 
    comparisons tests because it is not possible to apply such control here.
    LW.c means that the values for LW are already corrected by the 
    temperature of the instrument, so they are good to go.

"""

# Dependencies
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt


# import pysolar

# # My modules
# from antarc.load_rad_flux import load_rad_flux
# from antarc.run_libradtran import run_libradtran

# Parameter modules
from antarc.marambio.parameters import rad_flux_params


def string_to_float(pdseries) -> list[float]:
    """
    Replace commas with periods in numbers represented as strings
    and convert from string to float
    @param pdseries: numbers given as strings, with some nans
    @return list of resulting floats
    """
    dum = list(pdseries)
    return [
        float(x.replace(",", ".")) if isinstance(x, str) else x for x in dum
    ]


dirname = rad_flux_params.RAD_FLUX_ALL_DIR
filename = rad_flux_params.RAD_FLUX_ALL_FILE
outdir = "/Users/prowe/Sync/measurements/Marambio/broadbandRadiation/stand/"

# Columns:
# "";"date";"T.pt100";"T.datalogger";"SWDir";"SWD.spn1";SWD.nr01";"SWU";
# "LWD.c";"LWU.c"
#
# The ones we will use are:
# "date";"SWDir";"SWD.spn1";SWD.nr01";"SWU"; "LWD.c";"LWU.c"
# "date": date as yyyy-mm-dd hh:mm:ss;
# "SWDir": Direct SW down?
# "SWD.spn1": Global SW Down?
# "SWD.nr01": Direct Normal SW
# "SWU": SW up
# "LWD.c": corrected LW down
# "LWU.c": corrected LW up
data = pd.read_csv(
    dirname + filename,
    delimiter=";",
    usecols=("date", "SWD.spn1", "SWU", "LWD.c", "LWU.c"),
)
#    usecols=("date", "SWDir", "SWD.spn1", "SWD.nr01", "SWU", "LWD.c", "LWU.c"),

# Write out the dates within a given date range to a file
start_date = "2022-02-01"
end_date = "2022-02-15"

# Remove the dates before start_date or after end_date
irem = np.union1d(
    np.where(data["date"] < start_date)[0],
    np.where(data["date"] > end_date)[0],
)
data = data.drop(index=irem)


# swdir = string_to_float(data["SWDir"])
sw_spn1 = string_to_float(data["SWD.spn1"])
# sw_nr01 = string_to_float(data["SWD.nr01"])
swu = string_to_float(data["SWU"])
lwd = string_to_float(data["LWD.c"])
lwu = string_to_float(data["LWU.c"])


# Plot the data to check it out
plt.figure()
# plt.plot(pd.to_datetime(data["date"]), swdir, label="SWDir")
plt.plot(pd.to_datetime(data["date"]), sw_spn1, label="SWD spn1")
# plt.plot(pd.to_datetime(data["date"]), sw_nr01, label="SWD nr01")
plt.plot(pd.to_datetime(data["date"]), swu, label="SWU")
plt.legend()

plt.figure()
plt.plot(pd.to_datetime(data["date"]), swu, label="SWDir")

plt.figure()
plt.plot(pd.to_datetime(data["date"]), lwd, label="SWD spn1")
# plt.plot(pd.to_datetime(data["date"]), lwu, label="SWD nr01")
plt.legend()

# Loop over days, collect data from that day, and
# write to a file
datefmt = "%Y-%m-%d %H:%M:%S"
# dtime = dt.datetime.strptime(data['date'][0], datefmt)
dtime_list = list(data["date"])
dtime = [dt.datetime.strptime(x, datefmt) for x in dtime_list]

# dtime0 = dtime[0]
# day = dtime0.day()
# thisday = dtime0.day()
# while day == thisday:
#     np.write(dtime_list[i], swdir, sw_spn1, sw_nr01, swu, lwd, lwu)
#     day = dtime0.day()

file_prefix = "mbio_broadband_"
file_ext = ".csv"
header = "date, swdir, sw_spn1, sw_nr01, swu, lwd, lwu\n"


i = 0
dtime0 = dtime[i]
fname = dtime[i].strftime(file_prefix + "%Y%m%d" + file_ext)
thisday = dtime[i].day
fid = open(outdir + fname, "w")
fid.write(header)
for i, dtime0 in enumerate(dtime):
    if dtime0.day > thisday:
        # close the old file, open the new
        fid.close()
        fname = dtime0.strftime(file_prefix + "%Y%m%d" + file_ext)
        thisday = dtime0.day
        fid = open(outdir + fname, "w")
        fid.write(header)
    elif dtime0.day < thisday:
        raise ValueError("Days should be increasing")
    # Write the data to the file
    # content = f"{dtime_list[i]}, {swdir[i]}, {sw_spn1[i]}, {sw_nr01[i]}, {swu[i]}, {lwd[i]}, {lwu[i]}\n"
    content = f"{dtime_list[i]}, {sw_spn1[i]}, {swu[i]}, {lwd[i]}, {lwu[i]}\n"
    fid.write(content)
fid.close()
