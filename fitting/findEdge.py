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
from scipy.signal import savgol_filter
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
    # print(np.diff(spectrum[1][:cutoff]))
    # yder = np.diff(spectrum[0][:cutoff])/np.diff(spectrum[1][:cutoff])
    # xder = np.array([])
    
    # for i in range(len(yder)-cutoff):
    #     print(len(yder)-cutoff)
    #     xtemp = (spectrum[0][i+1]+spectrum[0][i])/2
    #     xder = np.append(xder, xtemp)
    
    # return xder,yder
    # y = spectrum[1][:cutoff]
    filty = savgol_filter(spectrum[1][:cutoff], 301, 3)
    y = np.log10(filty)
    return np.gradient(filty, sampleRate)
    # return np.gradient(y, sampleRate)
    # return y


cut = 7000

# CEs = findEdge(spect, 1000, 1)
CEs = findEdge2(spect, 1, cut)
print(CEs)
# print(min(CEs[0][1]))
# yder = savgol_filter(np.exp(CEs[:cut]), 101, 3)#+abs(min(CEs[:cut]))
# yder = CEs[:cut]
yder = np.exp(CEs[:cut]*np.log(10))      #10 hoch CEs
yhat = savgol_filter(spect[1][:cut], 301, 3)
edgePoints = np.array([])
for i in range(len(yder)):
    if yder[i] < 0:
        edgePoints = np.append(edgePoints, spect[0][i])

plt.figure()
plt.yscale('log')
plt.plot(spect[0][:cut],yhat)
plt.plot(spect[0][:cut],yder)
# plt.plot(spect[0][:cut], yhat)
# plt.vlines(edgePoints,min(yder),max(yder),color='green')
# plt.plot(CEs[1],CEs[2])
# plt.vlines(CEs[0], 0.0, max(spect[1]), color="red")
plt.show()
