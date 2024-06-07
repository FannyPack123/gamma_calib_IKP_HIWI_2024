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
import os


spect = np.genfromtxt('../spectraOutput/06/Bi207/1300V_Offset0000_spect.csv', delimiter=',')

def findEdge(spectrum,sampleRate,cutoff):

    # samplRate = 2000
    
    # print(int(len(spect[0])/samplRate))
    
    edgeX = np.array([spect[0][0]])
    edgeY = np.array([spect[1][0]])
    edges = np.array([])
    
    for i in range(int(len(spect[0])/sampleRate)):
        
        # print("hello")
        start = int(spect[0][0]) + int(i*sampleRate)
        end = int(spect[0][0]) + int((i+1)*sampleRate)
        avg = np.average(spect[1][start:end])
        
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
    
    


CEs = findEdge(spect, 1000, 1)


plt.figure()
plt.yscale('log')
plt.plot(spect[0],spect[1])
plt.plot(CEs[1],CEs[2])
plt.vlines(CEs[0], 0.0, max(spect[1]), color="red")
plt.show()
