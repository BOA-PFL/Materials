# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 13:56:32 2023

Analyzing stretch testing results from Zwick Roell. File contains cycles of
stretch testing providing force and displacement of the data

"""

import pandas as pd
import os
from tkinter.filedialog import askopenfilename
import numpy as np
from matplotlib import pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal import butter,filtfilt
import scipy.signal as sig
from scipy import integrate
from tkinter import messagebox

#filename = askopenfilename()
fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Materials Testing\\Swatch Creation\\'
fileExt = r".csv"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]
check_data = 1
save_on = 1
outfileName = fPath + 'testResults.csv'
badFileList = []

for entry in entries:

    if entry.split(' ')[0].split('_')[-1] == 'Channels':
        
        dat = pd.read_csv(fPath + entry, sep=';', header = 0, skiprows=[1])
        dat.Force = dat.Force*1000
        
        fr = 500 #Hertz
        
        order = 2
        cutoff = 20
        nyq = 0.5 * fr
        
        def butter_lowpass_filter(data, cutoffVal, fs, order):
            normal_cutoff = cutoffVal / nyq
            # Get the filter coefficients 
            b, a = butter(order, normal_cutoff, btype='low', analog=False)
            y = filtfilt(b, a, data)
            return y
        
        normal_cutoff = cutoff/nyq
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        
        # extract force and displacement data alone
        ForceDat = dat.Force
        DispDat = dat.Displacement
        ForceDisp = ForceDat/DispDat
        
        # low pass filter both
        FilteredForceDat = butter_lowpass_filter(ForceDat,cutoff, fr, 2)
        filteredDatDisp = butter_lowpass_filter(DispDat,cutoff, fr, 2)
        FilteredForceDispDat = FilteredForceDat/filteredDatDisp
        
        # Check signal
        locs,_ = sig.find_peaks(-1 * FilteredForceDat, distance = 200)
        # find the local minima, filter out peaks below 10 and above 12 to find mean peak value
        # find an appropriate threshold value 
        pks = np.array(FilteredForceDat[locs])
        adaptiveThresh = pks[np.where((pks > 10) & (pks < 12))].mean() 
        padding = 0.5
        FilteredForceDat[FilteredForceDat < adaptiveThresh + padding] = adaptiveThresh
        #plt.plot(FilteredForceDat[5000:6000])
        
        minima = []
        stops = []
        for i, val in enumerate(FilteredForceDat):
            if i == (len(FilteredForceDat)-1):
                pass
            elif (val == adaptiveThresh) & (FilteredForceDat[i+1] > adaptiveThresh):
                minima.append(i)
            elif (val == adaptiveThresh) & (FilteredForceDat[i-1] > adaptiveThresh):
                stops.append(i)
                
        # get rid of first 5 indices of both since they seem to return unstable data
        if len(minima) == 0:
            # if no minima are detected, go back thru this file carefully
            print('No minima detected, adding to badFileList; check file '+entry)
            badFileList.append(entry)
        else:
            n = 70
            del minima[:n]
            del stops[:n]
            # if first element of stops is 'before' first of minima, delete first of stops
            if stops[0] < minima[0]:
                stops.pop(0)
            
            
            if check_data == 1:
                fig,ax = plt.subplots(1,1)    
                for i,val in enumerate(minima):
                    if i == len(minima)-1:
                        #ax.text(0.5,0,txt,transform=fig.transFigure)
                        ax.text(1.5,35,'Up', fontsize = 12,color = 'blue')
                        ax.text(1.5,30,'Down', fontsize = 12, color = 'orange')
                    else:
                        start = minima[i] 
                        mid = minima[i] + round((stops[i] - minima[i]) / 2)
                        end = stops[i]
                        
                        ax.plot(filteredDatDisp[start:mid],FilteredForceDat[start:mid], color = 'blue',label = 'up')
                        ax.plot(filteredDatDisp[mid:end],FilteredForceDat[mid:end], color = 'orange',label='down')
             
                
                answer = messagebox.askyesno("Question","Is data clean?")
                badFileList = []
                
                if answer == False:
                    plt.close('all')
                    print('Adding file to bad file list')
                    badFileList.append(entry)
                
                if answer == True:
                    plt.close('all')
                    print('Estimating point estimates')  
            else:
                answer = True
            
            
            EnergyLoss = []
            PercentReturn = []
            Stiffness = []
            
            
            for i, val in enumerate(minima):
                try:
                    tmpForce = FilteredForceDat[minima[i]:stops[i]]
                    tmpDispDat = filteredDatDisp[minima[i]:stops[i]]
                    tmpFDDat = FilteredForceDispDat[minima[i]:stops[i]]
                    tmpMax = np.argmax(tmpForce)
                    tmpMaxFP = np.argmax(tmpDispDat)
                    tmp25 = round(0.25*tmpMaxFP)
                    tmp75 = round(0.75*tmpMaxFP)
                    
                    # integration
                    rampup = integrate.trapz(tmpForce[0:tmpMax],tmpDispDat[0:tmpMax])
                    rampdown = -1 * integrate.trapz(tmpForce[tmpMax:-1],tmpDispDat[tmpMax:-1])
            
                    EnergyLoss = rampup - rampdown
                    PercentReturn.append( ( 1-(EnergyLoss/rampup))*100 )
                    Stiffness.append( tmpFDDat[tmp75] - tmpFDDat[tmp25]  )
                
                except Exception as e: print(e)
                
    
            arrayValues = [[int(entry.split('_')[1]),np.mean(PercentReturn), np.mean(Stiffness)]]
            col_vals = ['SpecimenNumber','PercentReturn', 'Stiffness']
            
            outcomes = pd.DataFrame(
                data = arrayValues,
                columns = col_vals)
            
    
            if save_on == 1:
                if os.path.exists(outfileName) == False:
                    
                    outcomes.to_csv(outfileName, mode='a', header=True, index = False)
                
                else:
                    outcomes.to_csv(outfileName, mode='a', header=False, index = False) 

