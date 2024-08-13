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


def halfGauss(xArr,sigma,mu,a,d):               #this one works!!!!!!!
    
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

popt, pcov = optimize.curve_fit(halfGauss, np.linspace(-50.0, 50.0, num = 41), yVals, pzero)
# print(*popt)

popt2, pcov2 = optimize.curve_fit(halfGauss, spectrum[0][2000:4500],spectrum[1][1500:4000], [500,3250,190000,30])
print(*popt2)

xFit = np.linspace(-50.0, 50.0, num = 201)
yFit = halfGauss(xFit,*popt)


yFit2 = halfGauss(spectrum[0][2000:4500],*popt2)
yFit2prime = halfGauss(spectrum[0][2000:4500],500,3250,190000,30)

# print(xVals)
# print(yFit)


#plot

plt.figure()
plt.title("fit1")
# plt.vlines([sigma,sigma*1.177], min(yVals), [gauss(sigma,sigma,mu),max(yVals)], color=['orange','red'], linestyles='dashed')
plt.plot(xVals,yVals,'.')
plt.plot(xFit, yFit, color = 'green')
plt.vlines(popt[0]*1.177, min(yFit), max(yFit), color='red', linestyles='dashed')


plt.figure()
plt.title("fit2")
# plt.yscale('log')
plt.plot(spectrum[0][2000:4500],spectrum[1][1500:4000])
plt.plot(spectrum[0][2000:4500],yFit2)
plt.plot(spectrum[0][2000:4500],yFit2prime)
plt.vlines(popt2[1]+popt2[0]*1.177,ymin=0.0, ymax=max(spectrum[1][1500:4000]),color="red",linestyles="dashed")
plt.show()







