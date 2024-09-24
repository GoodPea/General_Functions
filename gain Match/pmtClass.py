

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed July  11 08:26:26 2024
The purpose of this file is to define a class to read in processed files into
a class for conveient plotting and sorting of data.
@author: ethanschneider
"""

'''
Components of name in order:
DataR_CH0@DT5730S_14167 -> get channel from here
date
source
bar
location
notes


Components of OGS name  in order:
DataR_CH0@DT5730S_14167_7-10-2024_bkg_bar1_pmt394.csv

...while in directory with root file...
root <filename.root>
TBrowser b
'''
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from functions import *
from scipy.signal import savgol_filterls

from scipy.signal import find_peaks
from scipy.interpolate import interp1d

###############################################################################
class histClass:
    ''' A class to store information of filenames and processed pulses'''
    def __init__(self, fullFile, name_scheme):
        if isinstance(fullFile, list):
            objList = fullFile
            self.aging = objList[0].aging
            self.coating = objList[0].coating
            self.hist = objList[0].hist
            self.error  = objList[0].error**2
            self.date = objList[0].date
            objNum = 0
            for obj in objList:
                if objNum == 0 :
                    objNum += 1
                    continue
                self.hist = self.hist + obj.hist
                self.error = self.error + obj.error**2
            self.hist = self.hist/len(objList)
            self.error = np.sqrt(self.error)
            return

        COMP_POS = 1
        fileName = os.path.basename(fullFile)
        self.file = fileName
        nameList = fileName.split('_')
        compassName = nameList[COMP_POS]
        self.channel = int(compassName[2])

        if name_scheme == "pmtTest":
            DATE_POS = 3
            SOURCE_POS = 4
            COMP_POS = 1
            BAR_POS = 5
            pmtNumPos = 6
            if (nameList[5])[:7] == "control":
                self.bar = None
            elif nameList[5] != "control":
                self.bar = int((nameList[5])[3])

            try:
                self.date = nameList[DATE_POS]
            except:
                self.date = 0
            try:
                self.source = nameList[SOURCE_POS]
            except:
                self.source = 0
            try:
                pmtNum = int((nameList[pmtNumPos])[3:6])
                self.pmtNum = pmtNum
            except:
                self.pmtNum = "control"
            try:
                self.barNum = int((nameList[BAR_POS])[3])
            except:
                self.barNum = None


        data = np.loadtxt(fullFile, delimiter = ";")
       #  self.tails = data[:,0]
        self.integral = data
        # self.amp = data[:,2]
        self.psd = None
        self.notes = os.path.splitext(nameList[-1])[0]
        self.hist = None
        self.hist_err = None
        self.bins = None
        self.gaussParm = None
        self.gaussError = None
        self.MeasTime = None
        self.energyRes = None
        self.comptEdge = None
        self.height = None
        self.depth = None
        self.width = None
        self.comptEdge = None
        self.comptEdgeHeight = None
        self.comptEdgeChange = None
        self.FWHM = None
        return





    def clear_data(self):
        self.tails = 0
        self.integral = 0
        self.amp = 0
    def print_data(self):
        print("Date: " + str(self.date))
        print("Souce: " + str(self.source))
        print("Notes: " + str(self.notes))
        print("Channel: " + str(self.channel))
        print("Position: " + str(self.position))
        print("Bar Number: " + str(self.barNum))

    def print_cube(self):
        print("date " + str(self.date))
        print("Gain Match: " + str(self.gainMatch))
        print("Souce: " + str(self.source))
        print("Aged: " +str(self.aging))

    def copy_aging_coating(self, source):
    # Ensure the source is an instance of the same class
        if not isinstance(source, histClass):
            raise ValueError("Source must be an instance of MyClass")

        if hasattr(source, 'coating'):
            self.coating = source.coating
        if hasattr(source, 'aging'):
            self.aging = source.aging



    def fit_gauss(self, lower, upper):
        chunk = self.hist[(lower):upper]
        bins = self.bins[(lower-(upper-lower)):upper]
        peak = np.append(np.flip(chunk), chunk)
        guess = [np.max(peak),
                 np.argmax(peak),
                 np.std(peak)]
        plt.plot(bins, peak)
        plt.plot(self.bins[:-1], self.hist)
        plt.show()
        self.gaussParm, self.ErrorMatrix = curve_fit(gauss, bins, peak, p0 = guess)
        plt.plot(bins, gauss(bins, self.gaussParm[0],
                                  self.gaussParm[1], self.gaussParm[2]) )
        plt.plot(self.bins[:-1], self.hist)

        return

    def FWHM_calc(self, a, b, c):

        bins = self.bins[0,:]

        lower = a
        upper = b
        chunk = self.hist[(lower):upper]
        bins1 = bins[lower:upper]
        guess = [np.max(chunk),
                 np.argmax(chunk),
                 np.std(chunk)]

        peak1, _  = curve_fit(gauss, bins1, chunk, p0 = guess)


        lower = b
        upper = c
        chunk = self.hist[(lower):upper]
        bins2 = bins[lower:upper]

        guess = [np.max(chunk),
                 np.argmax(chunk),
                 np.std(chunk)]
        peak2, _  = curve_fit(gauss, bins2, chunk, p0 = guess)

        self.FWHM = round(np.abs(peak1[1] - peak2[1])/(2.35*peak1[2] + 2.35 * peak2[2]),2)
        plt.plot(bins1, gauss(bins1, peak1[0], peak1[1], peak1[2]) )
        plt.plot(bins2, gauss(bins2, peak2[0], peak2[1], peak2[2]) )
        plt.plot(bins[:-1], self.hist)
        plt.title("Gauss Fit " + str(self.date) + str(self.FWHM))
        plt.show()

        return


    def set_hist(self, histData, bins):
        self.hist = histData
        self.bins = bins
        return
    
    def find_inflection(self): 
        hist = self.hist
        bins = self.hist
        smoothed_hist = savgol_filter(hist, 4, 1)
        second_dir = savgol_filter(np.gradient(np.gradient(smoothed_hist)), 4, 1)
        # Now to find the bounds of this range: 
            
        peaks = find_peaks(second_dir)
        # for peak in peaks[0]: 
        #     plt.axvline(bins[peak])
        # plt.axvline(bins[4])
        # plt.plot(bins[:-1], second_dir)
        # plt.show()
        
        start = peaks[0][0] + np.argmin(second_dir[peaks[0][0]:])    # Removes first peak
        # plt.plot(bins[start:-1], second_dir[start:])
        start = start + np.argmin(second_dir[start:])
        end = start + np.argmax(second_dir[start:])
        
        # Now find the zeroes in this range: 
        # Might be best to modify this to find the slop along this line. 
        # The maximum could be used to calculate the window
        
        xmin = start + np.argmin(np.abs(second_dir[start:end]))
        slope = -1
        k = 1 
        while slope < 0: 
            # checks to see if the slope is positive or negative, if the slope is 
            # negative following an iteration, the range on which the code looks at
            # for interpolation increases 
            if np.abs(second_dir[xmin+k]) < np.abs(second_dir[xmin-k]):
                x = [bins[xmin], bins[xmin+k]]
                y = [second_dir[xmin], second_dir[xmin+k]]
            elif np.abs(second_dir[xmin+k]) > np.abs(second_dir[xmin-k]): 
                x = [bins[xmin-k], bins[xmin]]
                y = [second_dir[xmin-k], second_dir[xmin]]
            else: 
                # no need for interpolation 
                return bins[xmin]
            slope, b = line(x[0], y[0], x[1], y[1])
            k += 1 

            
        inflection = -b/slope 
        plt.plot(bins[:-1], hist, label = "Origional")
        plt.plot(bins[:-1], smoothed_hist, label = "Smoothed Hist")
        plt.plot(bins[:-1], second_dir, label = "Second Derivative")
        plt.axvline(inflection)
        plt.ylabel("Counts/second")
        plt.xlabel("Pulse Integral (V*ns)")
        plt.legend()
        plt.show()
        return inflection 
