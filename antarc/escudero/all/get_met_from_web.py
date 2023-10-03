#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:22:11 2022

@author: prowe

"""

from antarc.met_data_get_from_web import get_met_data_for_year
import json
import requests


def get_met_15_minutes(station: float, dtime, outdir: str):
    """
    Get met data from metochile on 15 minute intervals
    @param station  The station number
    @param dtime  The desired date (just year and month)
    """

    daykeys = {
        "temperatura": "temperature (C),",
        "puntoDeRocio": "dewpoint (C),",
        "humedadRelativa": "RH (%),",
        "presionEstacion": "Pressure (hPa),",
        "presionNivelDelMar": "Sea level pressure (hPa),",
        "presionNivelEstandar": "Standard level pressure (hPa),",
        "aguaCaidaDelMinuto": "Precipitation by minute (mm),",
        "aguaCaida6Horas": "Precipitation over 6 hours (mm),",
        "aguaCaida24Horas": "Precipitation over 24 hours (mm),",
        "direccionDelViento": "Wind direction (degree),",
        "fuerzaDelViento": "Wind speed (knot),",
    }

    station = str(station)
    datestr = dtime.strftime("/%Y/%m?")

    "https://climatologia.meteochile.gob.cl/application/servicios/getDatosRecientesEma/950001/2019/11?usuario=prowe@harbornet.com&token=92508d823a380131dc33a3b3"
    # response = json.loads(requests.get("https://climatologia.meteochile.gob.cl/application/geoservicios/getEmaResumenDiarioGeo?usuario=correo@correo.cl&token=apiKey_personal").text)
    url = (
        "https://climatologia.meteochile.gob.cl/application/servicios/"
        + "getDatosRecientesEma/"
        + str(station)
        + datestr
    )
    username = "prowe@harbornet.com"
    apikey = "92508d823a380131dc33a3b3"
    try:
        response = json.loads(
            requests.get(url + "usuario=" + username + "&token=" + apikey).text
        )
    except:
        # If there is no response then presumably there is no data
        # for this month, so quit
        return

    fname = outdir + "esc_met_" + dtime.strftime("%Y%m.txt")
    with open(fname, "w") as fid:
        # fid.write("Year,Month,Day,Hour,Minute,Second,")
        fid.write("Date,")
        for key, value in daykeys.items():
            fid.write(value)
        fid.write("\n")

        for resp in response["datosEstaciones"]["datos"]:
            # date_n_time = resp["momento"].split()
            # ymd = date_n_time[0].split("-")
            # hms = date_n_time[1].split(":")
            # fid.write(ymd[0] + "," + ymd[1] + "," + ymd[2] + ",")
            # fid.write(hms[0] + "," + hms[1] + "," + hms[2] + ",")
            fid.write(resp["momento"] + ",")
            for key in daykeys:
                if resp[key] is None:
                    fid.write("nan,")
                else:
                    fid.write(resp[key].split()[0] + ",")
            fid.write("\n")


# # #     BY HOUR     # # # #
# Columns: Spanish to English
# old_header_to_new = {
#     "Momento UTC ": "Date (UTC)",
#     "Temperatura (°c)": "Temperature (C)",
#     "Temperatura a 3 Metros (°c)": "Temperature at 3 m (C)",
#     "Punto de Rocio (°c)": "Dewpoint (C) ",
#     "Humedad (%)": "Humidity (%)",
#     "Razón de Mezcla (Gr/Kg)": "Mixing Ratio (g/kg)",
#     "Presión de la Estación (hPa)": "Station Pressure (hPa)",
#     "Presión a nivel del Mar (hPa)": "Sea Level Pressure (hPa)",
#     "Viento (°/Kph)": "Wind (deg/Kph)",
#     "Agua Caida Acumulada Diaria (mm)": "Accumulated Daily Precipitation (mm)",
# }
# year = 2016
# direc = esc_met_params.DIREC_BY_HOUR + str(year) + "/"
# get_met_data_for_year(
#     "by_hour",
#     esc_met_params.STATION,
#     old_header_to_new,
#     year,
#     direc,
#     esc_met_params.PFIX,
#     esc_met_params.EXT,
# )


def get_met(year: int, direc: str, station, prefix, ext):
    """
    Get Escudero met data from the web
    @param year  Desired year
    @param outdir  Desired location of output data
    Inputs, e.g.
    year = 2022
    direc = met_params.DIREC_BY_MINUTE + str(year) + "/"
    """
    # # # #     BY MINUTE     # # # #
    # Columns: Spanish to English
    old_header_to_new = {
        "HORA": "Time",
        "VIENTO (Grados/Nudos) DDInst": "Wind DDInst (deg)",
        "VIENTO (Grados/Nudos) FFInst": "Wind FFInst (knots)",
        "VIENTO (Grados/Nudos) DD2Min": "Wind DD2Min (deg)",
        "VIENTO (Grados/Nudos) FF2Min": "Wind FF2Min (knots)",
        "VIENTO (Grados/Nudos) DD10Min": "Wind DD10Min (deg)",
        "VIENTO (Grados/Nudos) FF10Min": "Wind FF10Min (knots)",
        "VIENTO (Grados/Nudos) FFMáximo": "Wind FFMax (knots)",
        "TEMPERATURA (°) Seco": "Dry Temperature (C)",
        "TEMPERATURA (°) Rocio": "Dewpoint Temperature (C)",
        "TEMPERATURA (°) TMin12Horas": "Temperature TMin12Hours (C)",
        "TEMPERATURA (°) TMax12Horas": "Temperature TMax12Hours (C)",
        "HUMEDAD (%)": "Humidity (%)",
        "PRESION (hPa) QFE": "Pressure QFE (hPa)",
        "PRESION (hPa) QFF": "Pressure QFF (hPa)",
        "PRESION (hPa) QNH": "Pressure QNH (hPa)",
        "AGUA CAÍDA (mm) 1 Min": "Precipitation 1 Min (mm)",
        "AGUA CAÍDA (mm) 6 Horas": "Precipitation 6 Hours (mm)",
        "AGUA CAÍDA (mm) 24 Horas": "Precipitation 24 Hours (mm)",
    }

    # get_met_data_for_year(
    #     "by_minute",
    #     met_params.STATION,
    #     old_header_to_new,
    #     year,
    #     direc,
    #     met_params.PFIX,
    #     met_params.EXT,
    # )

    get_met_data_for_year(
        "by_minute",
        station,
        old_header_to_new,
        year,
        direc,
        prefix,
        ext,
    )
