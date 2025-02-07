#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  26 09:33:02 2023

@author: prowe

ssrd
-----
This parameter is the amount of solar radiation (shortwave radiation) that 
reaches a horizontal plane at the surface of the Earth. This parameter 
comprises both direct and diffuse solar radiation.

Radiation from the Sun (solar, or shortwave, radiation) is partly reflected 
back to space by clouds and particles in the atmosphere (aerosols) and some 
of it is absorbed. The rest is incident on the Earth's surface (represented 
by this parameter). See further documentation.

To a reasonably good approximation, this parameter is the model equivalent 
of what would be measured by a pyranometer (an instrument used for measuring 
solar radiation) at the surface. However, care should be taken when comparing 
model parameters with observations, because observations are often local to 
a particular point in space and time, rather than representing averages over 
a model grid box.

This parameter is accumulated over a particular time period which depends 
on the data extracted. The units are joules per square metre (J m-2). 
To convert to watts per square metre (W m-2), the accumulated values should 
be divided by the accumulation period expressed in seconds. The ECMWF 
convention for vertical fluxes is positive downwards.
"""

from netCDF4 import Dataset, num2date
import numpy as np
import datetime as dt
import glob

# from typing import Optional


# def add_to_era5_data(
#     erafile: str, lat: Optional[float] = None, lon: Optional[float] = None
# ):
#     if lat and np.any(lat != nci["latitude"][:]):
#         raise ValueError("Latitudes changed in era5 file")
#     if lon and np.any(lon != nci["longitude"][:]):
#         raise ValueError("Longitudes changed in era5 file")
#     return


def get_and_pack_era5_data(erafile):
    (
        date,
        lat,
        lon,
        swd,
        lwd,
        swd_clr,
        lwd_clr,
        sw_net,
        lw_net,
        sw_net_clr,
        lw_net_clr,
    ) = get_era5_data(erafile)

    era5 = {}
    era5["date"] = date
    era5["lat"] = lat
    era5["lon"] = lon
    era5["swd"] = swd
    era5["lwd"] = lwd
    era5["swd_clr"] = swd_clr
    era5["lwd_clr"] = lwd_clr
    era5["sw_net"] = sw_net
    era5["lw_net"] = lw_net
    era5["sw_net_clr"] = sw_net_clr
    era5["lw_net_clr"] = lw_net_clr

    return era5


def get_era5_broadband(erafile: str):
    """
    Get SWD, LWD, SWD_clr and LWD_clr from ERA5 files,
    where
    # ssr: Surface net short-wave (solar) radiation
    # ssrc: Surface net short-wave (solar) radiation, clear sky
    # ssrd: Surface short-wave (solar) radiation downwards
    # ssrdc: Surface solar radiation downward clear-sky
    # str: Surface net long-wave (thermal) radiation
    # strc: Surface net long-wave (thermal) radiation, clear sky
    # strd: Surface long-wave (thermal) radiation downwards
    # strdc: Surface thermal radiation downward clear-sky
    """

    with Dataset(erafile, "r") as nci:

        date = num2date(
            nci["time"][:].data,
            nci["time"].units,
            only_use_python_datetimes=True,
            only_use_cftime_datetimes=False,
        )

        if "expver" in nci.variables.keys():
            # According to https://confluence.ecmwf.int/pages/
            # viewpage.action?pageId=173385064:
            # The ERA5 hourly and monthly data are made available with a
            # 3 month delay. This means that after a month has passed,
            # another month's worth of ERA5 data is written to the dataset.
            # ERA5T (near real time) preliminary data are used to fill
            # the gap between the end of the ERA5 data and 5 days before
            # the present date. The oldest month of these is overwritten
            # each month as new ERA5 data become available.
            # Both expver dimensions use the full time extent of time
            # coordinate but the expver 1 data only covers the first 7
            # timesteps, the remaining timesteps are 'padded' with empty
            # fields. For the expver 5 data, the first 7 timesteps are
            # padded with empty fields, with the remaining timesteps
            # coming from the ERA5T data.
            # When the last ERA5 data are released, they will overwrite the
            # ERA5T data for the entire month and for accumulated variables
            # for 00-06 in next month. This process will be repeated each
            # month.
            raise ValueError("Fix this!")
            # swd = np.vstack(
            #     [nci["ssrd"][:7, 0, :, :], nci["ssrd"][7:, 1, :, :]]
            # )
            # lwd = np.vstack(
            #     [nci["strd"][:7, 0, :, :], nci["strd"][7:, 1, :, :]]
            # )
            # swd_clr = np.vstack(
            #     [nci["ssrdc"][:7, 0, :, :], nci["ssrdc"][7:, 1, :, :]]
            # )
            # lwd_clr = np.vstack(
            #     [nci["strdc"][:7, 0, :, :], nci["strdc"][7:, 1, :, :]]
            # )

        else:
            # ssr: Surface net short-wave (solar) radiation
            # ssrc: Surface net short-wave (solar) radiation, clear sky
            # ssrd: Surface short-wave (solar) radiation downwards
            # ssrdc: Surface solar radiation downward clear-sky
            # str: Surface net long-wave (thermal) radiation
            # strc: Surface net long-wave (thermal) radiation, clear sky
            # strd: Surface long-wave (thermal) radiation downwards
            # strdc: Surface thermal radiation downward clear-sky

            era5 = {
                "swd": nci["ssrd"][:] / 3600,
                "lwd": nci["strd"][:] / 3600,
                "swd_clr": nci["ssrdc"][:] / 3600,
                "lwd_clr": nci["strdc"][:] / 3600,
                "swnet": nci["ssr"][:] / 3600,
                "lwnet": nci["str"][:] / 3600,
                "swnet_clr": nci["ssrc"][:] / 3600,
                "lwnet_clr": nci["strc"][:] / 3600,
            }

        # Radiation is accumulated over the hour and has units of J/m2
        # Therefore it must be divided by 3600 seconds to get W/m2
        era5["date"] = date
        era5["lat"] = nci["latitude"][:]
        era5["lon"] = nci["longitude"][:]

        return era5


def read_era5_broadband(
    direc: str, fmt: str, date1: dt.datetime, date2: dt.datetime
) -> tuple[np.ndarray, np.ndarray, np.ndarray,]:
    """
    Read ERA5 downwelling-to-surface broadband fluxes from files
    @param direc  Directory with data
    @param fmt  Format of data
    @param date1  First date in range to get
    @param date2  Last date in range to get
    Reminders:
    ncid.variables.keys()
    """
    # Get all filenames in between date1 and date2
    fmt = direc + "era5_esc_broadband_*[0-9]*.nc"
    files = glob.glob(fmt)

    if len(files) == 0:
        return None

    islash = files[0].rfind("/")
    files = [f[islash + 1 :] for f in files]

    iyear = files[0].find(str(date1.year))
    idot = files[0].find(".nc")
    dfmt = "%Y%m%d"

    if iyear == -1:
        raise ValueError(f"First file does not include {date1.year}")

    erafiles = []
    for fname in np.sort(files):
        thisdate = dt.datetime.strptime(fname[iyear:idot], dfmt)
        # thisdate = thisdate.replace(tzinfo=dt.timezone.utc)
        if date1 <= thisdate <= date2:
            erafiles.append(fname)

    if len(erafiles) == 0:
        return None

    # Set it up
    erafile = direc + erafiles[0]

    era5 = get_era5_broadband(erafile)

    # Note that the 0th file is done above, so start with 1 here
    for fname in erafiles[1:]:
        erafile = direc + fname
        era = get_era5_broadband(erafile)

        if np.any(era5["lat"] != era["lat"]):
            raise ValueError("Latitudes changed in era5 file")
        if np.any(era5["lon"] != era["lon"]):
            raise ValueError("Longitudes changed in era5 file")

        era5["date"] = np.hstack([era5["date"], era["date"]])
        era5["swd"] = np.vstack([era5["swd"], era["swd"]])
        era5["lwd"] = np.vstack([era5["lwd"], era["lwd"]])
        era5["swd_clr"] = np.vstack([era5["swd_clr"], era["swd_clr"]])
        era5["lwd_clr"] = np.vstack([era5["lwd_clr"], era["lwd_clr"]])

        era5["swnet"] = np.vstack([era5["swnet"], era["swnet"]])
        era5["lwnet"] = np.vstack([era5["lwnet"], era["lwnet"]])
        era5["swnet_clr"] = np.vstack([era5["swnet_clr"], era["swnet_clr"]])
        era5["lwnet_clr"] = np.vstack([era5["lwnet_clr"], era["lwnet_clr"]])

    return era5


def get_era5_data(erafile: str):
    """
    Get SWD, LWD, SWD_clr and LWD_clr from ERA5 files,
    where
    # ssr: Surface net short-wave (solar) radiation
    # ssrc: Surface net short-wave (solar) radiation, clear sky
    # ssrdc: Surface solar radiation downward clear-sky
    # ssrd: Surface short-wave (solar) radiation downwards
    # str: Surface net long-wave (thermal) radiation
    # strc: Surface net long-wave (thermal) radiation, clear sky
    # strd: Surface long-wave (thermal) radiation downwards
    # strdc: Surface thermal radiation downward clear-sky
    """

    with Dataset(erafile, "r") as nci:

        date = num2date(
            nci["time"][:].data,
            nci["time"].units,
            only_use_python_datetimes=True,
            only_use_cftime_datetimes=False,
        )

        if "expver" in nci.variables.keys():
            # According to https://confluence.ecmwf.int/pages/
            # viewpage.action?pageId=173385064:
            # The ERA5 hourly and monthly data are made available with a
            # 3 month delay. This means that after a month has passed,
            # another month's worth of ERA5 data is written to the dataset.
            # ERA5T (near real time) preliminary data are used to fill
            # the gap between the end of the ERA5 data and 5 days before
            # the present date. The oldest month of these is overwritten
            # each month as new ERA5 data become available.
            # Both expver dimensions use the full time extent of time
            # coordinate but the expver 1 data only covers the first 7
            # timesteps, the remaining timesteps are 'padded' with empty
            # fields. For the expver 5 data, the first 7 timesteps are
            # padded with empty fields, with the remaining timesteps
            # coming from the ERA5T data.
            # When the last ERA5 data are released, they will overwrite the
            # ERA5T data for the entire month and for accumulated variables
            # for 00-06 in next month. This process will be repeated each
            # month.

            def get_expver(val):
                return np.vstack([val[:7, 0, :, :], val[7:, 1, :, :]])

            swd = get_expver(nci["ssrd"])
            lwd = get_expver(nci["strd"])
            swd_clr = get_expver(nci["ssrdc"])
            lwd_clr = get_expver(nci["strdc"])

            swnet = get_expver(nci["ssr"][:])
            lwnet = get_expver(nci["str"][:])
            swnet_clr = get_expver(nci["ssrc"][:])
            lwnet_clr = get_expver(nci["strc"][:])

        else:
            swd = nci["ssrd"][:]
            lwd = nci["strd"][:]
            swd_clr = nci["ssrdc"][:]
            lwd_clr = nci["strdc"][:]

            swnet = nci["ssr"][:]
            lwnet = nci["str"][:]
            swnet_clr = nci["ssrc"][:]
            lwnet_clr = nci["strc"][:]

        # Radiation is accumulated over the hour and has units of J/m2
        # Therefore it must be divided by 3600 seconds to get W/m2
        return (
            date,
            nci["latitude"][:],
            nci["longitude"][:],
            swd / 3600,
            lwd / 3600,
            swd_clr / 3600,
            lwd_clr / 3600,
            swnet / 3600,
            lwnet / 3600,
            swnet_clr / 3600,
            lwnet_clr / 3600,
        )


def read_era5_broadband_down(
    direc: str, fmt: str, date1: dt.datetime, date2: dt.datetime
) -> tuple[np.ndarray, np.ndarray, np.ndarray,]:
    """
    Read ERA5 downwelling-to-surface broadband fluxes from files
    @param direc  Directory with data
    @param fmt  Format of data
    @param date1  First date in range to get
    @param date2  Last date in range to get
    Reminders:
    ncid.variables.keys()
    """
    # Get all filenames in between date1 and date2
    fmt = direc + "era5_esc_broadband_*[0-9]*.nc"
    files = glob.glob(fmt)

    if len(files) == 0:
        return None

    islash = files[0].rfind("/")
    files = [f[islash + 1 :] for f in files]

    iyear = files[0].find(str(date1.year))
    idot = files[0].find(".nc")
    dfmt = "%Y%m%d"

    if iyear == -1:
        raise ValueError(f"First file does not include {date1.year}")

    erafiles = []
    for fname in np.sort(files):
        thisdate = dt.datetime.strptime(fname[iyear:idot], dfmt)
        # thisdate = thisdate.replace(tzinfo=dt.timezone.utc)
        if date1 <= thisdate <= date2:
            erafiles.append(fname)

    if len(erafiles) == 0:
        return None

    # Set it up
    keys = ["swd", "lwd", "swd_clr", "lwd_clr"]
    keys += ["sw_net", "lw_net", "sw_net_clr", "lw_net_clr"]

    erafile = direc + erafiles[0]
    era5 = get_and_pack_era5_data(erafile)

    # Note that the 0th file is done above, so start with 1 here
    for fname in erafiles[1:]:
        if fname == "era5_esc_broadband_20231020.nc":
            print()
        erafile = direc + fname

        era5b = get_and_pack_era5_data(erafile)
        # date, lat, lon, swd, lwd, swd_clr, lwd_clr = get_era5_data(erafile)

        if np.any(era5["lat"] != era5b["lat"]):
            raise ValueError("Latitudes changed in era5 file")
        if np.any(era5["lon"] != era5b["lon"]):
            raise ValueError("Longitudes changed in era5 file")

        era5["date"] = np.hstack([era5["date"], era5b["date"]])
        for key in keys:
            era5[key] = np.vstack([era5[key], era5b[key]])
        # era5["swd"] = np.vstack([era5["swd"], swd])
        # era5["lwd"] = np.vstack([era5["lwd"], lwd])
        # era5["swd_clr"] = np.vstack([era5["swd_clr"], swd_clr])
        # era5["lwd_clr"] = np.vstack([era5["lwd_clr"], lwd_clr])

    return era5


def read_era5_broadband_down_utc(
    direc: str, fmt: str, date1: dt.datetime, date2: dt.datetime
) -> tuple[np.ndarray, np.ndarray, np.ndarray,]:
    """
    Read ERA5 downwelling-to-surface broadband fluxes from files
    assuming all times are in UTC, and thus dates are  offset-unaware
    @param direc  Directory with data
    @param fmt  Format of data
    @param date1  First date in range to get
    @param date2  Last date in range to get
    Reminders:
    ncid.variables.keys()
    """
    # Get all filenames in between date1 and date2
    fmt = direc + "era5_esc_broadband_*[0-9]*.nc"
    files = glob.glob(fmt)
    islash = files[0].rfind("/")
    files = [f[islash + 1 :] for f in files]

    iyear = files[0].find(str(date1.year))
    idot = files[0].find(".nc")
    dfmt = "%Y%m%d"

    if iyear == -1:
        raise ValueError(f"First file does not include {date1.year}")

    erafiles = []
    for fname in np.sort(files):
        thisdate = dt.datetime.strptime(fname[iyear:idot], dfmt)
        if date1 <= thisdate <= date2:
            erafiles.append(fname)

    if len(erafiles) == 0:
        return None

    era5 = {}
    # Set it up
    erafile = direc + erafiles[0]

    date, lat, lon, swd, lwd, swd_clr, lwd_clr, _, _, _, _ = get_era5_data(
        erafile
    )
    era5["date"] = date
    era5["lat"] = lat
    era5["lon"] = lon
    era5["swd"] = swd
    era5["lwd"] = lwd
    era5["swd_clr"] = swd_clr
    era5["lwd_clr"] = lwd_clr

    # Note that the 0th file is done above, so start with 1 here
    for fname in erafiles[1:]:
        erafile = direc + fname
        date, lat, lon, swd, lwd, swd_clr, lwd_clr, _, _, _, _ = get_era5_data(
            erafile
        )

        if np.any(era5["lat"] != lat):
            raise ValueError("Latitudes changed in era5 file")
        if np.any(era5["lon"] != lon):
            raise ValueError("Longitudes changed in era5 file")

        era5["date"] = np.hstack([era5["date"], date])
        era5["swd"] = np.vstack([era5["swd"], swd])
        era5["lwd"] = np.vstack([era5["lwd"], lwd])
        era5["swd_clr"] = np.vstack([era5["swd_clr"], swd_clr])
        era5["lwd_clr"] = np.vstack([era5["lwd_clr"], lwd_clr])

    return era5
