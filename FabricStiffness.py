# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:02:21 2023

@author: Adam.Luftglass
"""



# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 13:56:32 2023

Analyzing fabric stiffness testing results from Zwick Roell. File contains 
load data from one test of fabric stiffness

"""

import pandas as pd
import os

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import butter,filtfilt


#filename = askopenfilename()
fPath = 'C:\\Users\\adam.luftglass\\OneDrive - Boa Technology Inc\\General\\Materials Testing\\Fabric Stiffness\\117\\'
fileExt = r".csv"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]



save_on = 1
outfileName = fPath + 'testResults.csv'


fr = 100 #Hertz

order = 2
cutoff = 20
nyq = 0.5 * fr

def butter_lowpass_filter(data, cutoffVal, fs, order):
    normal_cutoff = cutoffVal / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

for entry in entries:

    #entry = entries[0]
    if entry.split(' ')[0].split('_')[-1] == 'Channels':

        
        dat = pd.read_csv(fPath + entry, sep=';', header = 0, skiprows=[1])
        dat.Force = dat.Force*1000
        
       
        normal_cutoff = cutoff/nyq
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        
    
        ForceDat = dat.Force
        
        # low pass filter 
        FilteredForceDat = butter_lowpass_filter(ForceDat,cutoff, fr, 2)

        #Plot the force profile        
        fig,ax = plt.subplots(1,1)  
        plt.plot(FilteredForceDat)
        plt.ylabel('Force (N)')
        plt.xlabel('Centiseconds (cs)')
      
                  
        arrayValues = [[int(entry.split('_')[1]),np.min(FilteredForceDat)]]
        col_vals = ['SpecimenNumber','Peak Load']
        outcomes = pd.DataFrame( data = arrayValues,columns = col_vals)
                
            
        if save_on == 1:
            if os.path.exists(outfileName) == False:
                            
                           outcomes.to_csv(outfileName, mode='a', header=True, index = False)
                        
            else:
                            outcomes.to_csv(outfileName, mode='a', header=False, index = False) 
        
