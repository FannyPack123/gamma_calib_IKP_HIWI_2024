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
import os


a = 5.0e-05
b = 0.1
c = 0.0
Ec = 3400.0
sigma = 500.0
d = 0.0

pzero = np.array([Ec,sigma,a,b,c,d])

def findEdge(spectrum,sampleRate,cutoff):   #algorithm to find comptomn edges by looking for jumps in data
    
    edgeX = np.array([spect[0][0]])
    edgeY = np.array([spect[1][0]])
    edges = np.array([])
    
    for i in range(int(len(spect[0])/sampleRate)):
        
        start = int(spect[0][0]) + int(i*sampleRate)
        end = int(spect[0][0]) + int((i+1)*sampleRate)
        avg = np.average(spect[1][start:end])
        
        edgeX = np.append(edgeX, start)
        edgeY = np.append(edgeY, avg)
        edgeX = np.append(edgeX, end)
        edgeY = np.append(edgeY, avg)
        
    for j in range(len(edgeX)-1):
        
        if math.log(edgeY[j],10) - math.log(edgeY[j+1],10) > 0.5 and edgeY[j] > cutoff:
            
            edges = np.append(edges, (edgeX[j+2]+edgeX[j])/2)
            
    return edges,edgeX,edgeY



def CEFunc(E,Ec,sigma,a,b,c,d):     #theoretical compton edge function
    
    # if a <=0:
    #     return 1000000000
    
    alpha = 0.5*(a*(E**2+sigma**2)+b*E+c)
    beta = ((-sigma)/(math.sqrt(2*math.pi))*a*(E+Ec)+b)
    exponent = (E-Ec)/(math.sqrt(2)*sigma)
    
    return alpha*special.erfc(exponent) + beta*np.exp(-(exponent**2)) + d


# xVals = np.linspace(-b/(2*a)-Ec/2, Ec+10000, num = 1000)

# yVals = CEFunc(xVals, *pzero)


spect = np.genfromtxt('../spectraOutput/06/Bi207/1300V_Offset0000_spect.csv', delimiter=',')    #read in data

CEs = findEdge(spect, 1000, 1)


outer = 500

plt.figure()
plt.plot(spect[0],spect[1])
plt.yscale('log')

for i in range(len(CEs[0:-1])):     #cycle through all comptom edges and fit data to CE function
    
    fitStart = int(CEs[0][i]) - 2000
    fitEnd = int(CEs[0][i])  + 2000
    
    CE_params = optimize.curve_fit(CEFunc, spect[0][fitStart:fitEnd], spect[1][fitStart:fitEnd], np.append(CEs[0][0], pzero[1:]))
    
    print(CE_params[0])
    print(CE_params[0][0]-CEs[0][1])
    
    CE = CEFunc(spect[0], *CE_params[0])
    
    
    # plt.plot(spect[0][fitStart-outer:fitEnd+outer],spect[1][fitStart-outer:fitEnd+outer])
    # plt.vlines(CEs[0], 0.0, max(spect[1]), color="purple")
    plt.plot(spect[0][fitStart:fitEnd], CE[fitStart:fitEnd])
    plt.vlines(CE_params[0][0], ymin=0, ymax=max(spect[1]), color='orange', label= "fitted Compton edge")
    
    plt.vlines(CEs[0], 0.0, max(spect[1]), color="purple", linestyles='dotted', label= "guesses")
    # fit2Start = 3500
    # fit2End = 7000
    
    # CE2_params = optimize.curve_fit(CEFunc, spect[0][fit2Start:fit2End], spect[1][fit2Start:fit2End], np.append(CEs[0][1]+1000, pzero[1:]))
    # print(CE2_params[0])
    # print(CE2_params[0][0]-CEs[0][1])
    
    # print(params[0])
    
    # xVals = np.linspace(-b/(2*a)-Ec/2, Ec+1000, num = 1000)
    
    # CE2 = CEFunc(spect[0], *CE2_params[0])
    
    # yVals = CEFunc(spect[0], *pzero)
plt.legend()
plt.savefig("testFitting.pdf", format="pdf", dpi=400)
plt.show()


# print(xVals)
# print(yVals)

# print(spect)



# plt.figure()
# plt.yscale('log')
# plt.plot(spect[0][fit1Start-outer:fit1End+outer],spect[1][fit1Start-outer:fit1End+outer])
# plt.plot(spect[0][fit1Start:fit1End], CE1[fit1Start:fit1End])
# plt.vlines(CE1_params[0][0], ymin=0, ymax=max(spect[1]), color='orange', linestyles='dotted')
# plt.show()

# plt.figure()
# plt.yscale('log')
# plt.plot(spect[0][fit2Start-outer:fit2End+outer],spect[1][fit2Start-outer:fit2End+outer])
# plt.plot(spect[0][fit2Start:fit2End], CE2[fit2Start:fit2End])
# plt.vlines(CE2_params[0][0], ymin=0, ymax=max(spect[1]), color='orange', linestyles='dotted')
# plt.show()


# plt.figure()
# plt.yscale('log')
# plt.plot(spect[0],spect[1])


# plt.plot(spect[0][fit1Start:fit1End], CE1[fit1Start:fit1End])
# plt.vlines(CE1_params[0][0], ymin=0, ymax=max(spect[1]), color='red', linestyles='dotted')
# plt.plot(spect[0][fit2Start:fit2End], CE2[fit2Start:fit2End])
# plt.vlines(CE2_params[0][0], ymin=0, ymax=max(spect[1]), color='red', linestyles='dotted')


# plt.vlines(CEs[0], 0.0, max(spect[1]), color="purple")
# plt.plot(CEs[1],CEs[2])

# # plt.figure()
# # plt.plot(xVals, yVals)
# # plt.vlines(Ec, ymin=0, ymax=a*Ec**2+b*Ec+c, color='orange', linestyles='dotted')



plt.show()