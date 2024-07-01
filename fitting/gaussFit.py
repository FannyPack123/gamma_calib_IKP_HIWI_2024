# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 17:31:04 2024

@author: phili
"""

import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import special
from scipy import optimize

def gauss(x,sigma,mu):
    
    return (1/(sigma*math.sqrt(2*math.pi)))*np.exp(-0.5*((x-mu)/sigma)**2)


def halfGauss(xArr,sigma,mu,a,d):
    
    results = np.array([])
    
    for x in xArr:
        if x < mu:
            results = np.append(results, 0.0)
        else:
            results = np.append(results, (a/(sigma*math.sqrt(2*math.pi)))*np.exp(-0.5*((x-mu)/sigma)**2)+d)
    
    return results


dataFile = '../spectraOutput/06/Bi207/1300V_Offset0000_spect.csv'
spectrum = np.genfromtxt(dataFile, delimiter=',') 




#params

sigma = 10.0
mu = 0.0
d = 0.0
a = 1.0

p0 = np.array([sigma,mu])
pzero = np.append(p0*1.5, [a, d])

#test data

xVals = np.linspace(-50.0, 50.0, num = 41)
yVals = gauss(xVals, *p0)


#fit to data

popt, pcov = optimize.curve_fit(gauss, np.linspace(0.0, 50.0, num = 41), yVals, p0)
# print(*popt)

popt2, pcov2 = optimize.curve_fit(halfGauss, spectrum[0][2000:4500],spectrum[1][1500:4000], [3700,3250,50,50])
print(*popt2)

xFit = np.linspace(-50.0, 50.0, num = 201)
yFit = gauss(xFit,*popt)


yFit2 = halfGauss(spectrum[0][2000:4500],*popt2)

# print(xVals)
# print(yFit)


#plot

plt.figure()
# plt.vlines([sigma,sigma*1.177], min(yVals), [gauss(sigma,sigma,mu),max(yVals)], color=['orange','red'], linestyles='dashed')
plt.plot(xVals,yVals,'.')
plt.plot(xFit, yFit, color = 'green')
plt.vlines(popt[0]*1.177, min(yFit), max(yFit), color='red', linestyles='dashed')


plt.figure()
# plt.yscale('log')
plt.plot(spectrum[0][2000:4500],spectrum[1][1500:4000])
plt.plot(spectrum[0][2000:4500],yFit2)
plt.show()







