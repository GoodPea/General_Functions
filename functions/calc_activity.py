#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 20:00:33 2024
Function to calculate activity 
@author: ethanschneider
"""

import math
from datetime import datetime 

def calc_activity(date, source):
    # returns the new activity of a source following a series of decays. 
    # dates for check sources beside cf252 to verify.
    # does not distinguish between weak or strong check sources 
    if source == "cf252": 
        initialActivity = 5.63*10**-6
        halflife = 2.647 # in years 
        recieveDate = "08-05-2015"
    if source == "cs137":
        halflife = 30.5
        recieveDate = "08-05-2008"
        initialActivity = 1*10**-6
    if source == "na22": 
        halflife = 2.6
        initialActivity = 1*10**-6
        recieveDate = "08-05-2008"
    if source == "co60": 
        halflife = 5.3
        recieveDate = "08-05-2008"
        initialActivity = 1*10**-6
    recieveDate = datetime.strptime(recieveDate, "%d-%m-%Y")
    final = datetime.strptime(date, "%d-%m-%Y")
    timeDif = final - recieveDate
    t = timeDif.days/365 
    finalActivity = initialActivity*math.exp(-math.log(2)/halflife*t)
    finalActivity = finalActivity*37*10**9
    return finalActivity

