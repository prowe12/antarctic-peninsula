#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:22:11 2022

@author: prowe

"""

from antarc.met_data_get_from_web import get_met_data_for_year
from antarc.escudero.parameters import met_params


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
year = 2022
direc = met_params.DIREC_BY_MINUTE + str(year) + "/"

get_met_data_for_year(
    "by_minute",
    met_params.STATION,
    old_header_to_new,
    year,
    direc,
    met_params.PFIX,
    met_params.EXT,
)
