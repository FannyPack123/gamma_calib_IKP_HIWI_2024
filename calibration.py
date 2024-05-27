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
#from parse import genSpectrum



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

voltages = measTimeDict.keys()
offset = "0000"

for volt in voltages:
    print("Voltage: "+volt+"\n")
    sources = measTimeDict[volt].keys()
    for src in sources:
        print("Source: "+src+"\n")
        quickParse(volt, src, offset)
    
    
    
    
    
    
    
    
    

# genSpectrum(voltage="1300", src="Bi207", offset="0000")







