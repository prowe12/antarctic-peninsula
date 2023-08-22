#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 10:48:09 2023

@author: prowe
"""

# Modules
from antarc.escudero.case202202.get_atm_profs import get_atm_profs

# Parameter modules
from antarc.escudero.parameters import esc_params, atm_prof_params
from antarc.escudero.parameters import esc202202 as esc_case


get_atm_profs(atm_prof_params, esc_params, esc_case)
