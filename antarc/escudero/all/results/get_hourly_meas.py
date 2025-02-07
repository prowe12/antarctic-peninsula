#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:05:04 2024

@author: prowe
"""

from antarc.escudero.all.analysis_rsrc import ave_over_every_hour


def get_hourly_meas_fill(lwd_date, lwd, swd_date, swd, date1, date2):
    """
    Average measurements over an hour. All hours between date1 and date2
    should be represented in the result, with nans where data does not exist
    """
    nlwd, lwd_date_a, lwd_a = ave_over_every_hour(lwd_date, lwd, date1, date2)
    nswd, swd_date_a, swd_a = ave_over_every_hour(swd_date, swd, date1, date2)

    # QC
    if len(lwd_date_a) != len(swd_date_a):
        raise ValueError("LW and SW dates are not the same")

    if lwd_date_a != swd_date_a:
        raise ValueError("LW and SW dates are not the same")

    return {
        "date": swd_date_a,
        "lwd": lwd_a,
        "swd": swd_a,
        "lwd_pts": nlwd,
        "swd_pts": nswd,
    }


def get_hourly_meas(lwd_date, lwd, swd_date, swd, date1, date2):
    """
    Average measurements to an hour
    """
    nlwd, lwd_date_a, lwd_a = ave_over_every_hour(lwd_date, lwd, date1, date2)
    nswd, swd_date_a, swd_a = ave_over_every_hour(swd_date, swd, date1, date2)

    if len(lwd_date_a) == 0:
        return {
            "date": swd_date_a,
            "lwd": lwd_a,
            "swd": swd_a,
            "lwd_date": lwd_date_a,
            "swd_date": swd_date_a,
            "lwd_pts": nlwd,
            "swd_pts": nswd,
        }
    elif len(lwd_date_a) != len(swd_date_a):
        raise ValueError("Problem")
    else:
        return {
            "date": swd_date_a,
            "lwd": lwd_a,
            "swd": swd_a,
            "lwd_date": lwd_date_a,
            "swd_date": swd_date_a,
            "lwd_pts": nlwd,
            "swd_pts": nswd,
        }

    return None
