#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:21:08 2022

@author: prowe
"""

import datetime as dt
import pytz

DTIME = dt.datetime(2022, 2, 8, 8, 13, 0, tzinfo=pytz.utc)
DATE1 = dt.datetime(2022, 2, 6, 0, 0, 0, tzinfo=pytz.utc)
DATE2 = dt.datetime(2022, 2, 11, 0, 0, 0, tzinfo=pytz.utc)
DT_FOEHN1 = dt.datetime(2022, 2, 7, 0, 0, 0, tzinfo=pytz.utc)
DT_FOEHN2 = dt.datetime(2022, 2, 9, 0, 0, 0, tzinfo=pytz.utc)
DT_FOEHNMAX = dt.datetime(2022, 2, 8, 3, 0, 0, tzinfo=pytz.utc)


# dtfoehn1 = dt.datetime(2022, 2, 7, 21, 0, 0, tzinfo = pytz.utc)
# dtfoehn2 = dt.datetime(2022, 2, 8, 6, 0, 0, tzinfo = pytz.utc)
