#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:22:11 2022

@author: prowe

Based on code from:
    https://www.geeksforgeeks.org/convert-html-table-into-csv-file-in-python/
"""

# Importing the required modules
import requests
from datetime import datetime, timedelta
import pandas as pd
from bs4 import BeautifulSoup
import re


def gen_days(year: int) -> list:
    """
    Get a list of all the datetimes in the given year
    @param year
    @return A list of datetimes of all dates in the given year
    """
    start_date = datetime(year, 1, 1)
    end_date = datetime(year, 12, 31)
    curr_date = start_date
    dates = [start_date]
    while curr_date < end_date:
        curr_date += timedelta(days=1)
        dates.append(curr_date)
    return dates


def get_met_data_by_hour(
    url: str, old_header_to_new: dict, filename: str
) -> pd.DataFrame():
    """
    Get the metorological data from the web as HTML data
    and convert to pandas dataframe
    @param url  The URL for the met data
    @param old_header_to_new  dictionary of old and new headers
    @param filename  The filename to save to
    @return The met data as a pandas dataframe
    """
    # Get the URL and read in the content
    resp = requests.get(url)

    # Get header from HTML file
    data = []
    list_header = []
    soup = BeautifulSoup(resp.content, "html.parser")
    header = soup.find_all("table")[0].find("tr")

    for items in header:
        try:
            list_header.append(items.get_text())
        except:
            continue

    # Check if the header is as expected
    for colname in old_header_to_new:
        if colname not in list_header:
            print(f"Bad or missing data for {url}")
            return None

    html_table = soup.find_all("table")[0].find_all("tr")[1:]

    for element in html_table:
        sub_data = []
        for sub_element in element:
            try:
                sub_data.append(sub_element.get_text())
            except:
                continue
        data.append(sub_data)

    # Put the data into Pandas DataFrame
    met = pd.DataFrame(data=data, columns=list_header)

    # Remove empty columns
    met.drop("\n", axis=1, inplace=True)

    # Rename columns from Spanish to English
    met.rename(columns=old_header_to_new, inplace=True)

    # Converting Pandas DataFrame
    # into CSV file
    met.to_csv(filename)

    return met


def get_met_data_by_minute(
    url: str, old_header_to_new: dict, filename: str
) -> pd.DataFrame():
    """
    Get the metorological data from the web as HTML data
    and convert to pandas dataframe
    @param url  The URL for the met data
    @param old_header_to_new  dictionary of old and new headers
    @param filename  The filename to save to
    @return The met data as a pandas dataframe
    """
    # Get the URL and read in the content
    resp = requests.get(url)

    # Get header from HTML file
    data = []
    list_header = []

    soup = BeautifulSoup(resp.content, "html.parser")
    try:
        header1 = soup.find_all("tr", {"class": re.compile("bg*")})[0]
        header2 = soup.find_all("tr", {"class": re.compile("bg*")})[1]
    except:
        print(f"Bad or missing data for {url}")
        return None

    list_header = []
    rows = []
    for items in header1:
        if items.get_text() == "\n":
            continue
        # Handle multiple rows
        rowspan = (
            int(items.attrs["rowspan"]) if "rowspan" in items.attrs else 1
        )
        # Handle multiple columns
        colspan = (
            int(items.attrs["colspan"]) if "colspan" in items.attrs else 1
        )

        for icol in range(colspan):
            list_header.append(items.get_text())
            rows.append(rowspan)

    icol = 0
    for items in header2:
        if items.get_text() == "\n":
            continue
        if rows[icol] > 2:
            raise ValueError("Bad value for rowspan")
        if rows[icol] == 2:
            # This column spans two rows, so move on to the next column before
            # appending text
            icol += 1
        list_header[icol] += " " + items.get_text()
        icol += 1

    # Check if the header is as expected
    for colname in old_header_to_new:
        if colname not in list_header:
            print(f"Bad or missing data for {url}")
            return None

    html_table = soup.find_all("tr", {"class": re.compile("bg*")})[2:]

    for element in html_table:
        sub_data = []
        for sub_element in element:
            if sub_element.get_text() == "\n":
                continue
            sub_data.append(sub_element.get_text())
        data.append(sub_data)

    # Put the data into Pandas DataFrame
    met = pd.DataFrame(data=data, columns=list_header)

    # Rename columns from Spanish to English
    met.rename(columns=old_header_to_new, inplace=True)

    # Converting Pandas DataFrame
    # into CSV file
    met.to_csv(filename)

    return met


def get_met_data_for_year(
    meas_frequency: str,
    station: int,
    header: dict,
    year: int,
    direc: str,
    pfix: str,
    ext: str,
):
    """
    Get the metorological data from the web as HTML data
    and convert to pandas dataframe for an entire year, by looping over day
    @param station  Station ID number
    @param year  Desired year of data
    @param direc  The directory to save to
    @param pfix  The output file prefix
    @param ext  The output file extension
    """

    url = {
        "by_hour": "https://climatologia.meteochile.gob.cl/application/diario/datosEmaRecienteHorario/",
        "by_minute": "https://climatologia.meteochile.gob.cl/application/diario/datosDiariosEstacionAutomatica/",
    }
    runner = {
        "by_hour": get_met_data_by_hour,
        "by_minute": get_met_data_by_minute,
    }

    get_met_data = runner[meas_frequency]
    datelist = gen_days(year)

    for dtime in datelist:
        # Output file
        datestr_file = dtime.strftime("%Y%m%d")
        filename = direc + pfix + datestr_file + ext

        datestr_url = dtime.strftime("%Y/%m/%d")
        fullurl = url[meas_frequency] + str(station) + "/" + datestr_url
        get_met_data(fullurl, header, filename)
