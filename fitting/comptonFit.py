# -*- coding: utf-8 -*-
"""
Created on Mon May 27 15:57:16 2024

@author: phili
"""

import io
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import special
from scipy import optimize
from scipy.signal import savgol_filter, find_peaks
import os


def findEdge(spectrum,sampleRate,cutoff):
    scaling = 10
    delta = 30
    filty = savgol_filter(spectrum[1][:cutoff], 301, 3)
    der = savgol_filter(np.gradient(filty, sampleRate)*scaling, 401,1)
    peaks = find_peaks(-der, height=50, prominence=1, distance=10, width=100)
    peaksLoc = np.array([])
    peaksHgt = np.array([])
    muPeaks = np.array([])
    muLoc = np.array([])
    last = 0
    for i in range(len(peaks[0])):
        print("from " + str(last) +"to "+ str(peaks[0][i]))
        muPeak = find_peaks(filty[last:peaks[0][i]], height=10, prominence=1, distance=0.4*abs(peaks[0][i]-last), width=50)
        print(muPeak[0])
        print(len(muPeak[0]))
        if len(muPeak[0]) == 0:
            muPeaks = np.append(muPeaks, 0)
        else:
            muPeaks = np.append(muPeaks, muPeak[0][0]+last)
        last = peaks[0][i]
    for pI in peaks[0]:
        peaksLoc = np.append(peaksLoc, spectrum[0][pI])
        peaksHgt = np.append(peaksHgt, der[pI]/scaling)
    for mpI in muPeaks: 
        muLoc = np.append(muLoc, spectrum[0][int(mpI)])
        # for k in np.array(range(pI))[::-1]:
        #     if der[k] >= -delta and der[k] <= delta :
        #         muLoc = np.append(muLoc, spectrum[0][k])
        #         muInd = np.append(muInd, k)
        #         break
                    
    return np.stack((peaksLoc, peaks[0], peaksHgt, muLoc, muPeaks),axis=-1), der


def halfGauss(EArr,mu,sigma,a,d):               #this one works!!!!!!!
    
    results = np.array([])
    
    for E in EArr:
        if E < mu:
            results = np.append(results, 0.0)
        else:
            results = np.append(results, d+(a/(sigma*math.sqrt(2*math.pi)))*np.exp(-0.5*((E-mu)/sigma)**2))
    
    return results

mu0 = 500
sigma0 = 190000
a0 = 100000
d0 = 0

p0 = np.array([mu0,sigma0,a0])

def getCEdgeVals(voltage,src,offset):
    for ch in ["00","01","02","03","04","05","06","07","08","09","10"]:
        vals = np.array([])
        print("processing channel " + ch)
        
        dataFile = 'spectraOutput/'+ch+'/'+src+'/'+voltage+'V_Offset'+offset+'_spect.csv'
        
        dataBool = os.path.exists(dataFile)

        if not dataBool: #test if spectrum-files needed for fitting exist
            print("spectrum file does not exist, continuing")
            print("----------------------------\n")
            continue
        
        # if os.stat(dataFile).st_size/(1024*1024*1024) > 2:
        #     print("file bigger than 2GB, continuing")
        #     print("\n----------------------------")
        #     continue
        
        spect = np.genfromtxt(dataFile, delimiter=',') 
        EdgeGuesses = findEdge(spect, 0.01, 7000)
        for e in range(len(EdgeGuesses[0])):    
            print(EdgeGuesses[0][e])
        print(len(EdgeGuesses[0]))
        
        plt.figure()
        plt.title("channel no." + ch + ", " + src + ", " + voltage)
        plt.plot(spect[0][:7000],spect[1][:7000],color='C0')
        for Edge in EdgeGuesses[0]:     #cycle through all comptom edges and fit data to CE function
            
            print("EDGE: ")
            print(Edge)
            EdgeX = Edge[0]
            mu = Edge[3]
            print("MU GUESS:")
            print(mu)
            muInd = int(Edge[4])
            sigma = Edge[0] - Edge[3]
            print("SIGMA GUESS:")
            print(sigma)
            sigmaInd = int(Edge[1])
            print("SIGMA INDEX: ")
            print(sigmaInd)
            outer = sigmaInd - muInd
            
            print("OUTER: ")
            print(outer)
            fitStart = sigmaInd  - int(2*outer)
            fitEnd = sigmaInd  + int(3*outer)
            if fitStart < 0:
                print("/////////////FIT START LESS THAN ZERO/////////")
                continue
            print("FIT START: ")
            print(fitStart)
            print("FIT END: ")
            print(fitEnd)
            
            
            # plt.plot(spect[0][:7000], EdgeGuesses[1],color="black")   #plot derivative of spectrum
            # plt.ylim(-1000,1000)
            
            
            plt.yscale('log')
            
            
            plt.vlines(EdgeX, 0.0, max(spect[1]), color="purple", linestyles='dotted', label= "guesses,Edge")
            plt.vlines(mu, 0.0, max(spect[1]), color="magenta", linestyles='dotted', label= "guesses,mu")
            
            CE_params = optimize.curve_fit(halfGauss, spect[0][fitStart:fitEnd], spect[1][fitStart:fitEnd], [mu, sigma, a0, d0])
            
            print('FIT PARAMS:')
            print(CE_params[0])
            
            fitCEdge = CE_params[0][0]+CE_params[0][1]*1.177
            vals = np.append(vals, fitCEdge)
            
            CE = halfGauss(spect[0], *CE_params[0])
            # control = halfGauss(spect[0], *[mu, sigma, a0, CE_params[0][3]]) # d0])
            
            plt.hlines(1000, xmin=spect[0][fitStart:fitEnd][0], xmax=spect[0][fitStart:fitEnd][-1], color='red')

            plt.plot(spect[0][fitStart:fitEnd], CE[fitStart:fitEnd],color="orange")
            # plt.plot(spect[0][fitStart:fitEnd], control[fitStart:fitEnd],color="black")
            
            plt.vlines(fitCEdge, ymin=0, ymax=max(spect[1]), color='green', label= "fitted Compton edge")
            plt.hlines(1000, xmin=spect[0][fitStart:fitEnd][0], xmax=spect[0][fitStart:fitEnd][-1], color='red')
            
            plt.legend()
            
            
            
         # for i in range(len(EdgeGuesses[0])):     #cycle through all comptom edges and fit data to CE function
         #     outer = 1000
         #     # fitStart = int(list(spect[0]).index(CEs[0][i]) - outer*(1+i*.01))
         #     # fitEnd = int(list(spect[0]).index(CEs[0][i])  + outer*(1+i*.01))
         #     fitStart = int(EdgeGuesses[0][i])  - outer
         #     fitEnd = int(EdgeGuesses[0][i])  + outer
         #     if fitStart < 0:
         #         continue
         #     print(fitStart)
         #     print(fitEnd)
             
         #     plt.figure()
         #     plt.title("channel no." + ch + ", " + src + ", " + voltage)
         #     plt.plot(spect[0],spect[1],color='cyan')
         #     plt.yscale('log')
         #     plt.vlines(EdgeGuesses[0], 0.0, max(spect[1]), color="purple", linestyles='dotted', label= "guesses")
         #     # print(CEs[0][i])
         #     CE_params = optimize.curve_fit(halfGauss, spect[0][fitStart:fitEnd], spect[1][fitStart:fitEnd], np.append((EdgeGuesses[0][i]-p0[0]*1.177), p0))
             
         #     print(CE_params[0])
             
             
         #     CE = halfGauss(spect[0], *CE_params[0])
             
         #     plt.hlines(1000, xmin=spect[0][fitStart:fitEnd][0], xmax=spect[0][fitStart:fitEnd][-1], color='red')

         #     plt.plot(spect[0][fitStart:fitEnd], CE[fitStart:fitEnd],color="green")
         #     plt.vlines(CE_params[0][0]+CE_params[0][1]*1.177, ymin=0, ymax=max(spect[1]), color='orange', label= "fitted Compton edge")
             
             
         #     plt.legend()
        
        
        np.savetxt('fitOutput/'+ch+'_'+src+'_'+voltage+'V_Offset'+offset+'_CEdgeFitValues.csv', vals, delimiter=",")
        plt.savefig('fitOutput/'+ch+'_'+src+'_'+voltage+'V_Offset'+offset+'_fit.pdf', format="pdf", dpi=400)
        plt.show()