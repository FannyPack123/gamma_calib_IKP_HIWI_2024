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


def findEdge2(spectrum,sampleRate,cutoff):   #algorithm to find comptomn edges by looking for jumps in data
    
    edgeX = np.array([spectrum[0][0]])
    edgeY = np.array([spectrum[1][0]])
    edges = np.array([])
    
    for i in range(int(len(spectrum[0])/sampleRate)):
        
        start = int(spectrum[0][0]) + int(i*sampleRate)
        end = int(spectrum[0][0]) + int((i+1)*sampleRate)
        avg = np.average(spectrum[1][start:end])
        
        edgeX = np.append(edgeX, start)
        edgeY = np.append(edgeY, avg)
        edgeX = np.append(edgeX, end)
        edgeY = np.append(edgeY, avg)
        
    for j in range(len(edgeX)-1):
        
        if math.log(edgeY[j],10) - math.log(edgeY[j+1],10) > 0.5 and edgeY[j] > cutoff:
            
            edges = np.append(edges, (edgeX[j+2]+edgeX[j])/2)
            # edges = np.append( edges, edgeX[j])
            # edges = np.append( edges, j)
            
    return edges,edgeX,edgeY



def findEdge(spectrum,sampleRate,cutoff):
    scaling = 10
    filty = savgol_filter(spectrum[1][:cutoff], 301, 3)
    der = savgol_filter(np.gradient(filty, sampleRate)*scaling, 401,1)
    peaks = find_peaks(-der, height=10, prominence=1, distance=10, width=100)
    peaksLoc = np.array([])
    peaksHgt = np.array([])
    for i in peaks[0]:
        peaksLoc = np.append(peaksLoc, spectrum[0][i])
    for i in peaks[0]:
        peaksHgt = np.append(peaksHgt, der[i]/scaling)
    return np.stack((peaksLoc, peaksHgt),axis=-1), der



def halfGauss(EArr,mu,sigma,a,d):               #this one works!!!!!!!
    
    results = np.array([])
    
    for E in EArr:
        if E < mu:
            results = np.append(results, 0.0)
        else:
            results = np.append(results, (a/(sigma*math.sqrt(2*math.pi)))*np.exp(-0.5*((E-mu)/sigma)**2)+d)
    
    return results



def CEFunc(E,Ec,sigma,a,b,c,d):     #theoretical compton edge function
    
    # if a <=0:
    #     return 1000000000
    
    alpha = 0.5*(a*(E**2+sigma**2)+b*E+c)
    beta = ((-sigma)/(math.sqrt(2*math.pi))*a*(E+Ec)+b)
    exponent = (E-Ec)/(math.sqrt(2)*sigma)
    
    return alpha*special.erfc(exponent) + beta*np.exp(-(exponent**2)) + d



def logFunc(E,Ec,sigma,a,b,c,d):
    # print(CEFunc(E, Ec, sigma, a, b, c, d))
    return np.log(CEFunc(E, Ec, sigma, a, b, c, d))




# spect = np.genfromtxt('../spectraOutput/06/Bi207/1300V_Offset0000_spect.csv', delimiter=',')          #test
# EdgeGuesses = findEdge(spect, 0.01, 7000)
# print(EdgeGuesses[0])



# a = 5.0e-05
# b = 0.1
# c = 0.0
# Ec = 0.0
# sigma = 500.0
# d = 0.0

# pzero = np.array([Ec,sigma,a,b,c,d])

mu = 500
sigma = 190000
a = 30
d = 0



p0 = np.array([mu,sigma,a])


def getCEdgeVals(voltage,src,offset):
    for ch in ["00","01","02","03","04","05","06","07","08","09","10"]:
        
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
        print(EdgeGuesses[0])
        # print(EdgeGuesses[0])
        
        gamma = math.sqrt(a*10e3)*(2*math.pi*math.e)**-0.25
        
        for Edge in EdgeGuesses[0]:     #cycle through all comptom edges and fit data to CE function
            
            outer = gamma/math.sqrt(abs(Edge[1]))  #estimate of stddev of to be fitted gaussian, based on slope at Edge guess
            print(outer)
            # fitStart = int(list(spect[0]).index(CEs[0][i]) - outer*(1+i*.01))
            # fitEnd = int(list(spect[0]).index(CEs[0][i])  + outer*(1+i*.01))
            fitStart = int(Edge[0])  - int(1.5*outer)
            fitEnd = int(Edge[0])  + int(1.5*outer)
            if fitStart < 0:
                continue
            print(fitStart)
            print(fitEnd)
            
            plt.figure()
            plt.title("channel no." + ch + ", " + src + ", " + voltage)
            plt.plot(spect[0],spect[1],color='C0')
            plt.yscale('log')
            plt.vlines(Edge[0], 0.0, max(spect[1]), color="purple", linestyles='dotted', label= "guesses")
            
            CE_params = optimize.curve_fit(halfGauss, spect[0][fitStart:fitEnd], spect[1][fitStart:fitEnd], np.append((Edge[0]-p0[0]*1.177), p0))
            
            print(CE_params[0])
            
            
            CE = halfGauss(spect[0], *CE_params[0])
            
            plt.hlines(1000, xmin=spect[0][fitStart:fitEnd][0], xmax=spect[0][fitStart:fitEnd][-1], color='red')

            plt.plot(spect[0][fitStart:fitEnd], CE[fitStart:fitEnd],color="orange")
            plt.vlines(CE_params[0][0]+CE_params[0][1]*1.177, ymin=0, ymax=max(spect[1]), color='green', label= "fitted Compton edge")
            
            
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
        
        
        
        # plt.savefig("testFitting.pdf", format="pdf", dpi=400)
        plt.show()