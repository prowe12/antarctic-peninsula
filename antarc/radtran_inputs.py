#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 10:58:00 2022

@author: prowe
"""

import os
import datetime as dt


class Inputs:
    __slots__ = [
        "lblrtm_dir",
        "run_lblrtm",
        "prof_dir",
        "od_dir",
        "pmomfiles_liq",
        "pmomfiles_ice",
        "prof_sample_file",
        "prof_files",
        "prof_datenums",
        "iprofbeg",
        "iprofdate",
        "iprofhour",
        "iprofend",
        "solarSourceFunFile",
        "surfEmissDataFile",
        "lat",
        "lon",
        "alt",
        "v1",
        "v2",
        "dnu",
        "NSTR",
        "NLYR",
        "NTAU",
        "UTAU",
        "NPHI",
        "PHI",
        "IBCND",
        "UMU0",
        "PHI0",
        "ALBEDO",
        "SSALB",
        "FBEAM",
        "FISOT",
        "DTAUC",
        "WVNMLO",
        "WVNMHI",
        "USRTAU",
        "NUMU",
        "UMU",
        "USRANG",
        "LAMBER",
        "TEMIS",
        "PLANK",
        "ONLYFL",
        "TEMPER",
        "TTEMP",
        "BTEMP",
        "OUTFILE" "HEADER",
        "CLDLYR",
        "NMOM",
        "NCLDLYR",
        "PMOM",
    ]

    def __init__(self, radtran_params):
        ip = radtran_params.__dict__
        self.lblrtm_dir = ip["lblrtm_dir"]
        self.run_lblrtm = ip["run_lblrtm"]
        self.prof_dir = ip["prof_dir"]
        self.od_dir = ip["od_dir"]
        self.pmomfiles_liq = ip["pmomfiles_liq"]
        self.pmomfiles_ice = ip["pmomfiles_ice"]
        self.prof_sample_file = ip["prof_sample_file"]
        self.prof_datenums = []
        self.iprofbeg = (ip["iprofbeg"],)
        self.iprofdate = (ip["iprofdate"],)
        self.iprofhour = (ip["iprofhour"],)
        self.iprofend = (ip["iprofend"],)
        self.solarSourceFunFile = ip["solarSourceFunFile"]
        self.surfEmissDataFile = ip["surfEmissDataFile"]
        self.lat = ip["lat"]
        self.lon = ip["lon"]
        self.alt = ip["alt"]
        self.v1 = ip["v1"]
        self.v2 = ip["v2"]
        self.dnu = ip["dnu"]
        self.NSTR = ip["NSTR"]
        self.NLYR = ip["NLYR"]
        self.NTAU = ip["NTAU"]
        self.UTAU = ip["UTAU"]
        self.NPHI = ip["NPHI"]
        self.PHI = ip["PHI"]
        self.IBCND = ip["IBCND"]
        self.UMU0 = ip["UMU0"]
        self.PHI0 = ip["PHI0"]
        self.ALBEDO = ip["ALBEDO"]
        self.SSALB = ip["SSALB"]
        self.FBEAM = ip["FBEAM"]
        self.FISOT = ip["FISOT"]
        self.DTAUC = ip["DTAUC"]
        self.WVNMLO = ip["WVNMLO"]
        self.WVNMHI = ip["WVNMHI"]
        self.USRTAU = ip["USRTAU"]
        self.NUMU = ip["NUMU"]
        self.UMU = ip["UMU"]
        self.USRANG = ip["USRANG"]
        self.LAMBER = ip["LAMBER"]
        self.TEMIS = ip["TEMIS"]
        self.PLANK = ip["PLANK"]
        self.ONLYFL = ip["ONLYFL"]
        # self.TEMPER = ip['']
        # self.TTEMP = ip['']
        # self.BTEMP = ip['']
        # self.OUTFILE
        # self.HEADER = ip['']
        self.CLDLYR = ip["CLDLYR"]
        self.NMOM = ip["NMOM"]
        self.NCLDLYR = ip["NCLDLYR"]
        self.PMOM = ip["PMOM"]

    def prepInputs(self):

        # # # # # # # # # # # # # #
        # .. Atmospheric Profile Info
        #
        self.getFilesInDir(
            "prof_files",
            "prof_dir",
            "prof_sample_file",
            "iprofbeg",
            "iprofend",
        )
        self.getDno("prof_datenums", "prof_files", "iprofdate", "iprofhour")
        # self.getFieldVal('maxAtmosphereHt', 60)  #% km
        # # # # # # # # # # # # # #

    def getFilesInDir(self, filestr, dirstr, sampleFnameStr, ibegV, iendV):
        """
        If the field self.(filestr) does not exist, cd to dirstr and
        get all the files/folders that have names similar to sampleFnameStr,
        e.g.
            - name is same length as sampleFnameStr
            - name has same beginning (match to indices ibeg)
            - name has same ending (match to indices iend)


        Return the list as self.(filestr)
        """
        fnamelist = []
        ibeg = getattr(self, ibegV)[0]
        iend = getattr(self, iendV)[0]
        samplefname = getattr(self, sampleFnameStr)
        setattr(self, filestr, [])
        dirname = getattr(self, dirstr)

        if hasattr(self, filestr) == False or len(getattr(self, filestr)) <= 0:
            dirstuff = sorted(os.listdir(dirname))
            count = 0
            for i in range(len(dirstuff)):
                fname = dirstuff[i]
                # print(fname)
                if (
                    len(fname) == len(samplefname)
                    and fname[ibeg[0] : ibeg[1]]
                    == samplefname[ibeg[0] : ibeg[1]]
                    and fname[iend[0] : iend[1]]
                    == samplefname[iend[0] : iend[1]]
                ):
                    count += 1
                    fnamelist.append(fname)
            setattr(self, filestr, fnamelist)

        else:
            raise NameError("Variable " + filestr + " has size ")

        if len(fnamelist) == 0:
            msg = (
                "No files like "
                + samplefname
                + " found in "
                + "directory "
                + getattr(self, dirstr)
            )
            print("Inputs.py: Warning: " + msg)

    def getDno(self, datenumstr, namestr, iDate, iTime):
        fnames = getattr(self, namestr)
        if hasattr(self, datenumstr) and len(getattr(self, datenumstr)) > 0:
            # If the field exists, make sure it has same length as the names it
            # is supposed to correspond to (namestr)
            if len(getattr(self, datenumstr)) != len(fnames):
                raise NameError("Number of dates and file names not same.")

        else:
            idate = getattr(self, iDate)[0]
            itime = getattr(self, iTime)[0]
            if itime[1] - itime[0] == 6:
                format_str = "%Y%m%d%H%M%S"
            elif itime[1] - itime[0] == 4:
                format_str = "%Y%m%d%H%M"
            elif itime[1] - itime[0] == 2:
                format_str = "%Y%m%d%H"
            else:
                raise NameError(
                    "index to time in filename must have length",
                    "2, 4, or 6, for HH, HHMM or HHMMSS",
                )

            d = [
                dt.datetime.strptime(
                    x[idate[0] : idate[1]] + x[itime[0] : itime[1]], format_str
                )
                for x in fnames
            ]
            setattr(self, datenumstr, d)
