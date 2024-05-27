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

def genSpectrum(voltage,src,offset,tMeasSrc,tMeasBg):
    
    for ch in ["00","01","02","03","04","05","06","07","08","09","10"]:
        
        inputPathRaw = 'parseOutput/'+src+'_'+voltage+'_Offset'+offset+'/peak_ch'+ch+'.csv'
        inputPathBg = 'parseOutput/Untergrund_'+voltage+'_Offset'+offset+'/peak_ch'+ch+'.csv'
        
        outPath = "spectraOutput/" +ch+ '/' +src #+ '/' +voltage+ 'V_Offset'+offset

        rawBool = os.path.exists(inputPathRaw)
        bgBool = os.path.exists(inputPathBg)

        if not rawBool or not bgBool: #test if files needed for spectrum exist
            print("raw data file or background file does not exist")
            print("----------------------------\n")
            continue
        if os.path.isdir(outPath) == False:
            os.mkdir(outPath) #dynamically create output directory

        # tMeasSrc =  #measurement time of source in minutes
        # tMeasBg =   #measurement time of background in minutes
        #print(tMeasSrc)
        #print(tMeasBg)

        if tMeasBg == 0 or tMeasSrc == 0:
            print("source and/or background measurement time is set to 0!")
            print("----------------------------\n")
            continue
    
        #read in data
        specRead = np.genfromtxt(inputPathRaw, delimiter=',')
        bgRead = np.genfromtxt(inputPathBg, delimiter=',')
    
    
        #compute spectrum
        rawSpec = np.histogram(specRead, bins=int(np.ceil(np.amax(specRead))-np.amin(specRead)) )
        bg = np.histogram(bgRead, bins=int(np.ceil(np.amax(specRead))-np.amin(specRead)) )

        fullSpec = rawSpec[0]-(bg[0]*int(tMeasSrc/tMeasBg))

        # print(len(fullSpec))
        # print(len(rawSpec[1]))

        plt.figure()
        plt.title(src + ' gamma-spectrum, channel no. '+ch+', at '+voltage+'V')
        plt.ylabel('Num events')
        #plt.scatter(rawSpec[1][:-3], fullSpec[:-2], s=1)   #in case a spectrum with oveflow needs to be plotted,
                                                                    #overflow is truncated
        plt.scatter(rawSpec[1][:-1], fullSpec, s=1)
        plt.yscale('log')
        plt.savefig(outPath+ '/' +voltage+ 'V_Offset'+offset+'_spect.pdf', format="pdf", dpi=400)
        #plt.show()
        #vals = np.asarray([rawSpec[1][:-3], fullSpec[:-2]])
        vals = np.asarray([rawSpec[1][:-1], fullSpec])
        np.savetxt(outPath+ '/' +voltage+ 'V_Offset'+offset+'_spect.csv', vals, delimiter=",")
        print("\n----------------------------")
        
    return