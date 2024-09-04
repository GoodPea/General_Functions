#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 21:42:58 2024
Read a tally from an MCNP file 
@author: ethanschneider
"""


def read_mcnp(file, startline, endline): 
    # Read through these items and create an array with all of the necessary 
    # data 
    rfile = open(file, "r")
    lineNum = 1                 # line 
    i = 0 
    dataArray = np.zeros(((endline-startline+20),2))
    for line in rfile:
        if lineNum >= startline and lineNum < endline-5: 
            print(line)
            data = line.split(" " ) 
            data = list(filter(None, data))
            k=0
            while k < 2:
               dataArray[i,k] = float(data[k].lower())
               k+=1
            i+=1
        lineNum += 1 
    return dataArray
