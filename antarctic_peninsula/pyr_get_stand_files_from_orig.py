#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:18:25 2022

@author: prowe
"""

# Dependencies
import datetime as dt
import pandas as pd
import numpy as np


class StandFilesFromOrig:
    """
    Save the pyranometer data to standard format
    """

    def __init__(self, direc, fname, fname_fmt):
        self.direc = direc
        self.filename = fname
        self.filename_fmt = fname_fmt

        # self.qc_filedate()
        self.dataframe = self.read_data()

    def get_skiprows(self, firstcol: str = ".data") -> int:
        """
        Get the number of unused lines at the top
        @param filename  The name of the pyranometer file
        @param firstcol  Value found in the first column of the data
        @return  Number of lines to skip, including header
        @raises ValueError if firstcol value not found
        """
        nlines = 0
        with open(self.direc + self.filename, encoding="utf8") as file:
            try:
                for line in file:
                    thisline = line.split(";")
                    if thisline[0] == firstcol:
                        return nlines
                    nlines += 1
                return nlines
            except Exception as epn:
                raise ValueError(
                    f"While attempting to read {self.filename}: {epn}"
                ) from epn

    def get_headerline(self, firstcol: str = ".data.hdr"):
        """
        Get the number of unused lines at the top
        @param filename  The name of the pyranometer file
        @param firstcol  Value found in the first column of the data
        @return  Number of header lines
        @raises ValueError if firstcol value not found
        """
        headerfirstline = 0
        with open(self.direc + self.filename, encoding="utf8") as file:
            for line in file:
                thisline = line.split(";")
                if thisline[0] == firstcol:
                    return headerfirstline
                headerfirstline += 1
            raise ValueError(f"No data found in {self.filename}")

    def read_data(self):
        """
        Read in the dataframe and put into standard form
        """

        # We need to skip all rows until we get to the row that starts with
        # firstcol (e.g. ".data"), excluding the header line (see below)
        # Here we find the number of lines to skip
        nskip = self.get_skiprows(".data")

        # There are three header lines, but there may be any number of lines
        # before them. We want the 2nd header line. Here we get the first,
        # then increment it by one
        headerline = self.get_headerline(".data.hdr") + 1

        # Now we want to skip all rows except for the header row
        skiprows = [x for x in range(nskip) if x != headerline]
        return pd.read_csv(
            self.direc + self.filename,
            skiprows=skiprows,
            sep=";",
            usecols=[1, 2, 3, 4, 5, 6, 7],
        )

    def get_filedtime(self) -> dt.datetime:
        """
        Return the datetime from the filename
        @return the date
        """
        try:
            return dt.datetime.strptime(self.filename, self.filename_fmt)
        except Exception as epn:
            raise ValueError(
                f"Failed reading date from {self.filename}, check format."
            ) from epn

    def isempty(self):
        """
        Return True if the file has no data, else false
        @returns True if the file has no data, else false
        """
        if (
            self.dataframe.empty
            or np.all(np.isnan(self.dataframe["#Samples"]))
            or np.all(self.dataframe["#Samples"] == 0)
        ):
            return True
        return False

    def get_start_date(self, date_fmt: str, time_fmt: str) -> tuple:
        """
        Return the first time from the dataframe
        @return the first time
        @raises ValueError if the dataframe is empty
        """
        if self.isempty():
            raise ValueError(
                f"Cannot get first time: no data found in {self.direc + self.filename}"
            )
        datestr = self.dataframe["Date"][0] + " " + self.dataframe["Time"][0]
        return dt.datetime.strptime(datestr, date_fmt + " " + time_fmt)

    def get_output_filename(
        self, date_fmt: str, time_fmt: str, outfile_fmt: str
    ) -> str:
        """
        Get the output file name in the desired format
        @param outfile_fmt Desired output format
        """
        # The hours in the filename are not typically UTC, but the times in
        # the file are, so replace time with hours and minutes from the file
        dtime = self.get_start_date(date_fmt, time_fmt)

        # Create the outputfile from the corrected time
        return dt.datetime.strftime(dtime, outfile_fmt)

    def save_dataframe(self, outputfile: str):
        """
        Save the output data
        @param outputfile  The output file path
        """
        if self.isempty():
            print(f"Not saving: no data found in {self.direc + self.filename}")
            return
        self.dataframe.to_csv(outputfile)


class StandFilesFromOrigV1(StandFilesFromOrig):
    """
    Save the pyranometer data to standard format
    for the version 1 format for dates before 2022-02-07-1906.csv
    """

    def read_data(self):
        """
        Read in the dataframe and put into standard form
        """

        # We need to skip all rows until we get to the row that starts with
        # firstcol (e.g. ".data"), excluding the header line (see below)
        # Here we find the number of lines to skip
        nskip = super().get_skiprows(".data")

        # There are three header lines, but there may be any number of lines
        # before them. We want the 2nd header line. Here we get the first,
        # then increment it by one
        headerline = self.get_headerline(".data.hdr") + 1

        # Now we want to skip all rows except for the header row
        skiprows = [x for x in range(nskip) if x != headerline]

        return pd.read_csv(
            self.direc + self.filename,
            skiprows=skiprows,
            sep=";",
            usecols=[1, 2, 3, 4, 5, 6, 7],
        )


class StandFilesFromOrigV2(StandFilesFromOrig):
    """
    Load broadband data from the original format and save to standard format
    for the version 2 format from 2022-02-07-1906.csv on
    .text;Info;<>;#Samples;Radiation;<>;Body Temperature;Power;Status;#Samples;
    """

    def read_data(self):
        nskiprows = super().get_skiprows(".data")

        return pd.read_csv(
            self.direc + self.filename,
            header=None,
            sep=";",
            skiprows=nskiprows,
            usecols=[1, 2, 3, 4, 6, 7, 8],
            names=[
                "Date",
                "Time",
                "#Samples",
                "Radiation",
                "Temperature",
                "Power",
                "Status",
            ],
        )


class LwdFilesFromOrigV1(StandFilesFromOrigV1):
    """
    Load pyrgeometer data from the original format, convert to downward
    longwave radiation, and save to standard format
    """

    def __init__(self, direc, fname, fname_fmt):
        super().init(direc, fname, fname_fmt)

    # def qc_filename(self):
    #     """Make sure the filename has the correct format"""
    #     fname = self.filename
    #     fnamelen = len(self.filename)
    #     samplen = len(params.SAMPLEFNAME)
    #     if fnamelen != samplen:
    #         raise ValueError(
    #             f"Filename {fname} is length {fnamelen} but \
    #                            should be length {samplen}"
    #         )

    # def qc_filedate(self):
    #     """
    #     Make sure the date in the filename is within the allowed range
    #     """
    #     if dateno + timeno / 1000 < self.start_date:
    #         raise ValueError("File date")
