# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 16:59:03 2024

@author: phili
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy import constants


def ComptEdge(E):
    return E*constants.e*(1-(1/(1+((2*E*constants.e)/(constants.m_e*(constants.c)**2)))))

def linear(x,m,b):
    return m*x+b


offset = 8000


calPts = np.array([[1.127362179693208782e+04-offset,ComptEdge(1063.656e03)/(constants.e*1e03)],
                   [1.016383211481090257e+04-offset,ComptEdge(661.657e03)/(constants.e*1e03)],
                   [1.247595848535043115e+04-offset,ComptEdge(1274.537e03)/(constants.e*1e03)],
                   [9.576123815302453295e+03-offset,ComptEdge(511.0e03)/(constants.e*1e03)],
                   [9.609066059023494745e+03-offset,ComptEdge(569.698e03)/(constants.e*1e03)]])


popt = optimize.curve_fit(linear, calPts[:,0], calPts[:,1], [1e-2,0])

xArr = np.linspace(0,6000,num=100)

plt.figure()
plt.title("calibration curve channel 03 @1000V")
plt.scatter(calPts[:,0],calPts[:,1], s=20, marker="x", color="red")
plt.plot(xArr, linear(xArr, *popt[0]))
plt.xlabel("bin number")
plt.ylabel("Energy/keV")
plt.hlines(0, xmin=0, xmax=6000, color="black")
plt.vlines(-(popt[0][1]/popt[0][0]), ymin=-500, ymax=500, color="orange", linestyles="dashed",label="root at: "+str(-(popt[0][1]/popt[0][0])))
plt.legend()
plt.show()

print("calibration channel 03: "+ str(popt[0][0]) + "*bin number + "+str(popt[0][1]))



calPts = np.array([[8.803423952858949633e+03-offset,ComptEdge(1063.656e03)/(constants.e*1e03)],
                   [9.286469546716758487e+03-offset,ComptEdge(1770.228e03)/(constants.e*1e03)],
                   [8.936097542362716922e+03-offset,ComptEdge(1274.537e03)/(constants.e*1e03)]])


popt = optimize.curve_fit(linear, calPts[:,0], calPts[:,1], [1e-2,0])

xArr = np.linspace(0,4000,num=100)

plt.figure()
plt.title("calibration curve channel 01 @1000V")
plt.scatter(calPts[:,0],calPts[:,1], s=20, marker="x", color="red")
plt.plot(xArr, linear(xArr, *popt[0]))
# plt.xlim(7000,11000)
plt.xlabel("bin number")
plt.ylabel("Energy/keV")
plt.hlines(0, xmin=0, xmax=4000, color="black")
plt.vlines(-(popt[0][1]/popt[0][0]), ymin=-2000, ymax=2000, color="orange", linestyles="dashed",label="root at: "+str(-(popt[0][1]/popt[0][0])))
plt.legend()
plt.show()

print("calibration channel 01: "+ str(popt[0][0]) + "*bin number + "+str(popt[0][1]))


calPts = np.array([[1.228762528498455140e+03,ComptEdge(569.698e3)/(constants.e*1e03)],
                   [2.123110844937242291e+03,ComptEdge(1063.656e03)/(constants.e*1e03)]])


popt = optimize.curve_fit(linear, calPts[:,0], calPts[:,1], [1e-2,0])

xArr = np.linspace(0,4000,num=100)

plt.figure()
plt.title("calibration curve channel 01 @1300V")
plt.scatter(calPts[:,0],calPts[:,1], s=20, marker="x", color="red")
plt.plot(xArr, linear(xArr, *popt[0]))
# plt.xlim(7000,11000)
plt.xlabel("bin number")
plt.ylabel("Energy/keV")
plt.hlines(0, xmin=0, xmax=4000, color="black")
plt.vlines(-(popt[0][1]/popt[0][0]), ymin=-2000, ymax=2000, color="orange", linestyles="dashed",label="root at: "+str(-(popt[0][1]/popt[0][0])))
plt.legend()
plt.show()

print("calibration channel 01 @1300V: "+ str(popt[0][0]) + "*bin number + "+str(popt[0][1]))




calPts = np.array([[1000.0,1.4184275014914989],
                    [1300.0,0.519189077685885]])


popt = optimize.curve_fit(linear, calPts[:,0], calPts[:,1], [1e-2,0])

xArr = np.linspace(0,2000,num=100)

plt.figure()
plt.title("voltage calibration curve channel 01")
plt.scatter(calPts[:,0],calPts[:,1], s=20, marker="x", color="red")
plt.plot(xArr, linear(xArr, *popt[0]))
# plt.xlim(7000,11000)
plt.ylabel("calibCoeff[keV/bin number]")
plt.xlabel("voltage[V]")
plt.hlines(0, xmin=min(xArr), xmax=max(xArr), color="black")
plt.vlines(-(popt[0][1]/popt[0][0]), ymin=-10, ymax=10, color="orange", linestyles="dashed",label="root at: "+str(-(popt[0][1]/popt[0][0])))
plt.legend()
plt.show()

print("calibration voltage channel 01: scale="+ str(popt[0][0]) + "*voltage + "+str(popt[0][1]))







# calPts = np.array([[1000.0,(9.286469546716758487e03-8000.0)],
#                    [1300.0,(2.123110844937242291e03)]])


# popt = optimize.curve_fit(linear, calPts[:,0], calPts[:,1], [1e-2,0])

# xArr = np.linspace(0,2000,num=100)

# plt.figure()
# plt.title("voltage calibration curve channel 01")
# plt.scatter(calPts[:,0],calPts[:,1], s=20, marker="x", color="red")
# plt.plot(xArr, linear(xArr, *popt[0]))
# # plt.xlim(7000,11000)
# plt.ylabel("scale")
# plt.xlabel("voltage/V")
# plt.hlines(0, xmin=min(xArr), xmax=max(xArr), color="black")
# plt.vlines(-(popt[0][1]/popt[0][0]), ymin=-2000, ymax=2000, color="orange", linestyles="dashed",label="root at: "+str(-(popt[0][1]/popt[0][0])))
# plt.legend()
# plt.show()

# print("calibration voltage channel 01: scale="+ str(popt[0][0]) + "*voltage + "+str(popt[0][1]))



# delta = 1063.656e03-569.698e3

# calPts = np.array([[1000.0,delta/(9.286469546716758487e03-8.803423952858949633e03)],
#                    [1300.0,delta/(2.123110844937242291e03-1.228762528498455140e03)]])


# popt = optimize.curve_fit(linear, calPts[:,0], calPts[:,1], [1e-2,0])

# xArr = np.linspace(0,2000,num=100)

# plt.figure()
# plt.title("voltage calibration curve channel 01")
# plt.scatter(calPts[:,0],calPts[:,1], s=20, marker="x", color="red")
# plt.plot(xArr, linear(xArr, *popt[0]))
# # plt.xlim(7000,11000)
# plt.ylabel("delta_E/delta_bin")
# plt.xlabel("voltage/V")
# plt.hlines(0, xmin=min(xArr), xmax=max(xArr), color="black")
# plt.vlines(-(popt[0][1]/popt[0][0]), ymin=-2000, ymax=2000, color="orange", linestyles="dashed",label="root at: "+str(-(popt[0][1]/popt[0][0])))
# plt.legend()
# plt.show()

# print("calibration voltage channel 01: delta_E/delta_bin="+ str(popt[0][0]) + "*voltage + "+str(popt[0][1]))
