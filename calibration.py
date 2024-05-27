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



measTimeDict= { '1000': {'Bi207':       10,
                         'Cs137':       10,
                         'Na22':        10,
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


offset = "0000"
parseBool = False
genBool = True


voltages = measTimeDict.keys()

#parsing
if parseBool:
    for volt in voltages:
        
        print("parsing for voltage: "+volt+"\n")
        sources = measTimeDict[volt].keys()
        
        for src in sources:
            
            print("parsing for source: "+src+"\n")
            quickParse(volt, src, offset)

#spectrum generation
if genBool:
    for volt in voltages:
        
        print("generating spectra for voltage: "+volt+"\n")
        sources = np.asarray(list(measTimeDict[volt].keys()))
        
        for src in sources[:-1]:
            
            print("generating spectra for source: "+src+"\n")
            genSpectrum(volt, src, offset, measTimeDict[volt][src],measTimeDict[volt]['Untergrund'])
            

if not parseBool and not genBool:
    print("parsing and spectrum generation is set to False")



















