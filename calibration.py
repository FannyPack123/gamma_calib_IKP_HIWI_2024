# -*- coding: utf-8 -*-
"""
Created on Fri May 24 23:27:21 2024

@author: phili
"""

import io
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os

from parse import quickParse
from parse import genSpectrum
from fitting import getCEdgeVals


    #gamma peak vals for sources, NUDAT
srcGammaEnergyDict = { 'Bi207'  : np.array([569.698e3,1063.656e3,1770.228e3]),  #328.10e3,511.0e3,897.77e3,1442.2e3 intensity too low
                       'Cs137'  : np.array([661.657e3]),                        #283.5e3
                       'Na22'   : np.array([511.0e3,1274.537e3])                        #511.0e3
                      }








measTimeDict= { '1000': {'Bi207':       30,
                         'Cs137':       30,
                         'Na22':        30,
                         'Untergrund':  214                                      
                         },
                '1100': {'Bi207':       10,
                         'Cs137':       10,
                         'Na22':        10,
                         'Untergrund':  302                                      
                         },
                '1200': {'Bi207':       10,
                         'Cs137':       10,
                         'Na22':        10,
                         'Untergrund':  1029                                      
                         },
                '1300': {'Bi207':       10,
                         'Cs137':       10,
                         'Na22':        10,
                         'Untergrund':  277                                      
                         },
                '1400': {'Bi207':       10,
                         'Cs137':       10,
                         'Na22':        10,
                         'Untergrund':  1276                                      
                         },
                '1500': {'Bi207':       10,
                         'Cs137':       10,
                         'Na22':        10,
                         'Untergrund':  1804                                      
                         }
               }


offset = "8000"
parseBool = True
genBool = True


voltages = measTimeDict.keys()

#parsing
# if parseBool:
#     for volt in voltages:
        
#         print("parsing for voltage: "+volt+"\n")
#         sources = measTimeDict[volt].keys()
        
#         for src in sources:
            
#             print("parsing for source: "+src+"\n")
#             quickParse(volt, src, offset)

#spectrum generation
# if genBool:
#     for volt in voltages:
        
#         print("generating spectra for voltage: "+volt+"\n")
#         sources = np.asarray(list(measTimeDict[volt].keys()))
        
#         for src in sources[:-1]:
            
#             print("generating spectra for source: "+src+"\n")
#             genSpectrum(volt, src, offset, measTimeDict[volt][src],measTimeDict[volt]['Untergrund'])
            

# if not parseBool and not genBool:
#     print("parsing and spectrum generation is set to False")
    
    
#compton edge fitting
for volt in ['1300']:
    
    print("getting compton edges for voltage: "+volt+"\n")
    sources = np.asarray(list(measTimeDict[volt].keys()))
    
    for src in ['Bi207']:
        
        print("source: "+src+"\n")
        getCEdgeVals(volt, src, offset)     

for volt in voltages:
    
    print("getting compton edges for voltage: "+volt+"\n")
    sources = np.asarray(list(measTimeDict[volt].keys()))
    
    for src in sources:
        
        print("source: "+src+"\n")
        getCEdgeVals(volt, src, offset)   

#calculate calibration coefficient
#TODO: 
#1. evaluate points (fitEdge, empiricalEnergy)
#2. fit linear function to find calib-coeff for each voltage
#3. use that info to find gain-curve

calPts = {}

for ch in ['00','01','02','03','04','05','06','07','08','09','10',]:
    calPts[ch] = {}
    
    for volt in ['1300','1000']:
        
        sources = np.asarray(list(measTimeDict[volt].keys()))
        
        for src in sources:
            dataFile = 'fitOutput/'+ch+'_'+src+'_'+volt+'V_Offset'+offset+'_CEdgeFitValues.csv'
            dataBool = os.path.exists(dataFile)

            if not dataBool:
                print("Edge file does not exist, continuing")
                print("----------------------------\n")
                continue
            CEdges = np.array([np.genfromtxt(dataFile, delimiter=',')])
            
            for i in range(len(CEdges)):
                calPts[ch][volt] = np.array(CEdges[i], srcGammaEnergyDict[src][i])

print(calPts)

















