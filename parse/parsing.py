# -*- coding: utf-8 -*-
"""
Created on Fri May 24 23:14:27 2024

@author: phili
"""

import io
import tools.parse as parse
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os

def quickParse(voltage,src,offset):
    
    # src = 'Bi207' #Bi207, Cs137, Na22, Untergrund
    # voltage = '1750' #0500, 0750, 1000, 1250, 1500, 1750
    outPath = 'parseOutput/'+src+'_'+voltage+'_Offset'+offset
    
    if os.path.isdir(outPath) == False:
        os.mkdir(outPath) #dynamically create output directory


    for ch in ["00","01","02","03","04","05","06","07","08","09","10"]: #cycle through all channels (detectors)
        
        print("processing channel " + ch)
        
        dataFile = 'data/Kalibrierung_2024/'+offset+'/'+src+'_'+voltage+'/raw-ch'+ ch +'.dat'
        #bgFile = 'data/Kalibrierung_2024/'+offset+'/Untergrund_'+voltage+'/raw-ch'+ ch +'.dat'

        dataBool = os.path.exists(dataFile)

        if not dataBool: #test if files needed for spectrum exist
            print("raw data file does not exist, continuing")
            print("----------------------------\n")
            continue
        
        
        
        if os.stat(dataFile).st_size/(1024*1024*1024) > 2:  #currently still needed, unoptimized code
            print("file bigger than 2GB, continuing")
            print("\n----------------------------")
            continue
        
        
        
        #parsing from Dougs code
        
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('-f','--file', type=str, default=dataFile, help='default=data/ch00.dat')
        # parser.add_argument('-v', '--verbose', action='store_true', help='Prints out many binary file components')
        # parser.add_argument('-es', '--eventSpacing', action='store_true', help='Plot histogram of event spacing')
        parser.add_argument('-N', '--N', type=int, default = 0, help='Plots first [N] events in the data file')
        parser.add_argument('-adc', '--adc', action='store_true', help='ADC histogram')
        args = parser.parse_args()
        
        args.adc = True
    
        infile = io.open(args.file, 'rb')
        p = parse.Parse(infile)
        # event_raw = []
        # event_maw = []
        # e_max = []
        timestamp = []
        true_time = []
        peak = []
        # maw_max = []
        totalEvents = 0
        HARDCODED_STOP = 1000 # Prevents over plotting of raw adc scope signals  
        infostr = ''
        
        for eventNum, event in enumerate(p):
            if args.N and eventNum >= args.N:
                break
            timestamp.append(event.ts)
            true_time.append( event.ts/(250_000_000) )
            
            totalEvents = eventNum
            
            # if hasattr(event,'e_max'):
                # e_max.append(event.e_max)
            if hasattr(event,'peak'):
                peak.append(event.peak)
            # if hasattr(event,'maw_max'):
                # maw_max.append(event.maw_max)
            # if args.verbose and eventNum<HARDCODED_STOP:
            #     if hasattr(event,'fmt'):
            #         print('format', bin(event.fmt))
            #     if hasattr(event,'peak'):
            #         print('peak ',event.peak)
            #     if hasattr(event,'e_max'):
            #         print('e_max ',event.e_max)
            #     if hasattr(event, 'acc7'):
            #         print('acc7 ', event.acc7)
            #     if hasattr(event, 'acc8'):
            #         print('acc8 ', event.acc8)
            #     if hasattr(event,'maw_max'):
            #         print('maw_max ',event.maw_max)
            # if eventNum<HARDCODED_STOP and hasattr(event,'raw'):
            #     event_raw.append( event.raw )
            # if eventNum<HARDCODED_STOP and hasattr(event,'maw'):
            #     event_maw.append( event.maw )
        
        print('total events: ', totalEvents+1)
        infostr = infostr + ('total events: ' + str(totalEvents+1) + '\n')
        if totalEvents == 0:
            print("no events or just one event in file, continuing to the next...")
            print("\n----------------------------")
            continue
    
    
        if peak and args.adc:
            
            plt.figure()
            plt.title(src + ' peak max, channel no. '+ch+', at '+voltage+'V')
            plt.ylabel('Num of events')
            
            # plt.hist(peak, range=[ np.amin(peak), int(np.ceil(np.amax(peak))) ],  bins=int(np.ceil(np.amax(peak))-np.amin(peak)) )
            plt.hist(peak,  bins=2**14)
            
            mostCommon = np.argmax(np.bincount(peak))
            print('most common value: ' + str(mostCommon))
            infostr = infostr + 'most common value: ' + str(mostCommon)+'\n'
            if mostCommon == (2**14)-1:
                print("-----------\nsuspected overflow!!!\n-----------")
                infostr = infostr + "-----------\nsuspected overflow!!!\n-----------"
                
            plt.yscale('log')
            
            plt.savefig(outPath+'/peak_ch'+ch+'.pdf', format="pdf", dpi=400)
            info = open(outPath+'/info_ch'+ch+'.txt', 'w')                      #verbesserbar
            info.write(infostr)
            info.close()
            vals = np.asarray(peak)
            np.savetxt(outPath+'/peak_ch'+ch+'.csv', vals, delimiter=",")
            
    
        print("\n----------------------------")
        
    return
