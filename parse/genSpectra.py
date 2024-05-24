# -*- coding: utf-8 -*-
"""
Created on Sat May 25 00:04:39 2024

@author: phili
"""

import io
import tools.parse as parse
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os

def genSpectrum(voltage,src,offset):
    
    for ch in ["00","01","02","03","04","05","06","07","08","09","10"]:
        
        inputPathRaw = 'parseOutput/'+src+'_'+voltage+'_Offset'+offset+'/peak_ch'+ch+'.csv'
        inputPathBg = 'parseOutput/Untergrund_'+voltage+'_Offset'+offset+'/peak_ch'+ch+'.csv'

        rawBool = os.path.exists(inputPathRaw)
        bgBool = os.path.exists(inputPathBg)

        if not rawBool or not bgBool: #test if files needed for spectrum exist
            print("raw data file or background file does not exist")
            print("----------------------------\n")
            continue
        if os.path.isdir(ch+'/'+src) == False:
            os.mkdir(ch+'/'+src) #dynamically create output directory

        tMeasSrc =  #measurement time of source in minutes
        tMeasBg =   #measurement time of background in minutes
        #print(tMeasSrc)
        #print(tMeasBg)

        if tMeasBg == 0 or tMeasSrc == 0:
            print("source and or background measurement time is set to 0!")
            print("----------------------------\n")
            continue
    return