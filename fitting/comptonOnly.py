# -*- coding: utf-8 -*-
"""
Created on Wed May 29 11:01:40 2024

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
b = 0.000001
c = 0.0
Ec = 1000.0
sigma = 100.0
d = 100.0


pzero = np.array([Ec,sigma,a,b,c,d])



def rFunc(E, Ec,sigma,a,b,c,d):
    
    result = np.array([])
    for energ in E:
        if energ <= Ec:
            result = np.append(result, (a*energ**2+b*energ+c+d))
            # print(result)
        else:
            result = np.append(result, 0.0+d)
            
    return result

def CEFunc(E,Ec,sigma,a,b,c,d):
    
    # Ec = pars[0]
    # sigma = pars[1]
    # a = pars[2]
    # b = pars[3]
    # c = pars[4]
    
    
    alpha = 0.5*(a*(E**2+sigma**2)+b*E+c)
    beta = ((-sigma)/(math.sqrt(2*math.pi))*a*(E+Ec)+b)
    exponent = (E-Ec)/(math.sqrt(2)*sigma)
    
    return alpha*special.erfc(exponent) + beta*np.exp(-(exponent**2)) + d

# [6.42255981e+03,  4.29010415e+02,  8.72830534e-06, -6.96115164e-02, 1.61983915e+02, 4.88246304e+00]

Ec = 6.42255981e+03
sigma = 4.29010415e+02
a = 8.72830534e-06
b = -6.96115164e-02
c = 1.61983915e+02
d = 4.88246304e+00

pzero = np.array([Ec,sigma,a,b,c,d])

# xVals = np.linspace(-b/(2*a)-Ec/2, Ec+4000, num = 1000)
xVals = np.linspace(4000, 9000, num = 1000)

# [ 6.42255981e+03,  4.29010415e+02,  8.72830534e-06, -6.96115164e-02, 1.61983915e+02, 4.88246304e+00]

CEVals = CEFunc(xVals, *pzero)
rVals = rFunc(xVals, *pzero)



print(rVals)

# print(xVals)
# print(yVals)

plt.figure()
plt.yscale('log')
plt.plot(xVals, CEVals)
# plt.plot(xVals, rVals)
plt.ylim(bottom = 0.0, top = 1000.0)
plt.vlines(Ec, ymin=0, ymax=max(rVals), color='orange', linestyles='dotted')










