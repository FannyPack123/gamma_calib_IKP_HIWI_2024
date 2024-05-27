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


a = 0.0001
b = 0.001
c = 0.0
Ec = 3400.0
sigma = 800
d = -100.0


pzero = np.array([Ec,sigma,a,b,c,d])




def CEFunc(E,*pzero):
    
    # Ec = pars[0]
    # sigma = pars[1]
    # a = pars[2]
    # b = pars[3]
    # c = pars[4]
    
    
    alpha = 0.5*(a*(E**2+sigma**2)+b*E+c)
    beta = -sigma/(math.sqrt(2*math.pi)*a*(E+Ec)+b)
    exponent = (E-Ec)/(math.sqrt(2)*sigma)
    
    return alpha*special.erfc(exponent) + beta*np.exp(-(exponent**2)) + d


# xVals = np.linspace(-b/(2*a)-Ec/2, Ec+10000, num = 1000)

# yVals = CEFunc(xVals, Ec, sigma, a, b, c)


spect = np.genfromtxt('../spectraOutput/06/Bi207/1300V_Offset0000_spect.csv', delimiter=',')

params = optimize.curve_fit(CEFunc, spect[0][1300:3500], spect[1][1300:3500], pzero)
# print(params[0])

xVals = np.linspace(-b/(2*a)-Ec/2, Ec+1000, num = 1000)

yVals = CEFunc(spect[0], *params[0])
# yVals = CEFunc(spect[0], *pzero)



# print(xVals)
# print(yVals)

# print(spect)

plt.figure()
# plt.yscale('log')
plt.plot(spect[0][1300:3500],spect[1][1300:3500])

# plt.figure()
plt.plot(spect[0][1300:3500], yVals[1300:3500])
# plt.vlines(Ec, ymin=0, ymax=a*Ec**2+b*Ec+c, color='orange', linestyles='dotted')

# plt.figure()
# plt.plot(xVals, yVals)
# plt.vlines(Ec, ymin=0, ymax=a*Ec**2+b*Ec+c, color='orange', linestyles='dotted')



plt.show()