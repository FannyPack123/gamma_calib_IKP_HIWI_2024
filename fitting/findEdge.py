# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:28:14 2024

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


spect = np.genfromtxt('../spectraOutput/06/Bi207/1300V_Offset0000_spect.csv', delimiter=',')

def findEdge(spectrum,sampleRate,cutoff):

    # samplRate = 2000
    
    # print(int(len(spect[0])/samplRate))
    
    edgeX = np.array([spectrum[0][0]])
    edgeY = np.array([spectrum[1][0]])
    edges = np.array([])
    
    for i in range(int(len(spectrum[0])/sampleRate)):
        
        # print("hello")
        start = int(spectrum[0][0]) + int(i*sampleRate)
        end = int(spectrum[0][0]) + int((i+1)*sampleRate)
        avg = np.average(spectrum[1][start:end])
        
        # print(f"interval: [",start,", ",end,"]")
        # print(avg)
        
        edgeX = np.append(edgeX, start)
        edgeY = np.append(edgeY, avg)
        edgeX = np.append(edgeX, end)
        edgeY = np.append(edgeY, avg)
        
    for j in range(len(edgeX)-1):
        
        # print(math.log(edgeY[j],10) - math.log(edgeY[j+1],10))
        
        if math.log(edgeY[j],10) - math.log(edgeY[j+1],10) > 0.5 :
            # print("true")
            if edgeY[j] > cutoff:
                edges = np.append( edges, edgeX[j])
        # else:
            # print("false")
            
            
    return edges,edgeX,edgeY

    # print(edges)
    
    # edgeX = np.append(edgeX, [spect[0][-1]])
    # edgeY = np.append(edgeY, [spect[1][-1]])
    
    

def findEdge2(spectrum,sampleRate,cutoff):
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



cut = 7000

# CEs = findEdge(spect, 1000, 1)
CEs = findEdge2(spect, .01, cut)
# print(min(CEs[0][1]))
# yder = savgol_filter(np.exp(CEs[:cut]), 101, 3)#+abs(min(CEs[:cut]))
yder = CEs[1][:cut]
edgePoints = np.array([])
for edgeVal in CEs[0]:
    edgePoints = np.append(edgePoints,edgeVal[0])
    
muPoints = np.array([])
for muVal in CEs[0]:
    muPoints = np.append(muPoints,muVal[3])

# yder = np.exp(CEs[0][:cut]*np.log(10))      #10 hoch CEs
yhat = savgol_filter(spect[1][:cut], 301, 3)



plt.figure()
# plt.yscale('log')
plt.ylim(-1300, 1000)
plt.plot(spect[0][:cut],yhat)
plt.plot(spect[0][:cut],yder)
plt.vlines(edgePoints, min(yder), max(yder), color="red")
plt.vlines(muPoints, min(yder), max(yder), color="magenta")
# plt.plot(spect[0][:cut], yhat)
# plt.vlines(edgePoints,min(yder),max(yder),color='green')
# plt.plot(CEs[1],CEs[2])
# plt.vlines(CEs[0], 0.0, max(spect[1]), color="red")
plt.show()
