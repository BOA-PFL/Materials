

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 13:56:32 2023

Analyzing stretch testing results from Zwick Roell. File contains cycles of
stretch testing providing force and displacement of the data

"""
import os
import pandas as pd
import numpy as np
import scipy.signal as sig
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from tkinter import messagebox, Tk, simpledialog

# Initialize Tkinter root window (hide it)
root = Tk()
root.withdraw()  # Hides the root window as we only need the dialog

# Ask the user to input the desired output file name using a dialog box
outfileName = simpledialog.askstring("Input", "Please enter the desired output file name (without extension):")
if not outfileName:
    print("No filename entered. Exiting...")
    exit()
outfileName += '.csv'

# Path to the directory containing the CSV files

fPath = 'C:\\Users\\adam.luftglass\\OneDrive - BOA Technology Inc\\General\\Materials Testing\\alphaflyspeedlandtesting\\'

outfilePath = os.path.join(fPath, outfileName)  # Combine path and file name
fileExt = r".csv"  # File extension of the target files
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]  # List of CSV files in the directory

check_data = 1  # Flag to check the data
save_on = 1  # Flag to save the results
badFileList = []  # List to track files with issues
fr = 500  # Sampling frequency in Hertz

order = 2  # Order of the filter
cutoff = 20  # Cutoff frequency for the filter
nyq = 0.5 * fr  # Nyquist frequency

def butter_lowpass_filter(data, cutoffVal, fs, order):
    """
    Apply a low-pass Butterworth filter to the data.
    """
    normal_cutoff = cutoffVal / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

# Process each file in the entries list
for entry in entries:

   #entry = entries[3]

    if entry.split(' ')[0].split('_')[-1] == 'Channels':
        # Read the CSV file
        dat = pd.read_csv(fPath + entry, sep=';', header=0, skiprows=[1])
        dat.Force = dat.Force * 1000  # Convert force to correct units



        # Extract force and displacement data
        ForceDat = dat.Force *-1
        DispDat = dat.Displacement * -1

        ForceDisp = ForceDat / DispDat


        # Apply low-pass filter to the data
        FilteredForceDat = butter_lowpass_filter(ForceDat, cutoff, fr, 2)
        filteredDatDisp = butter_lowpass_filter(DispDat, cutoff, fr, 2)
        FilteredForceDispDat = FilteredForceDat / filteredDatDisp

        # Identify local minima in the force data


        locs_min, _ = sig.find_peaks(-1 * FilteredForceDat, distance=2500, prominence = 1500)
        pks_min = np.array(FilteredForceDat[locs_min])
        
        
        
        
        locs_max,_ = sig.find_peaks(FilteredForceDat, distance=2500, prominence = 1500)
        pks_max = np.array(FilteredForceDat[locs_max])
        
      
        
        
        if locs_max[0] < locs_min[0]:
            locs_max = locs_max[1:]
        if locs_max[-1] < locs_min[-1]:
            locs_min = locs_min[:-1]
            
          
                  
        answer = True
        if check_data == 1:
             
             plt.figure()
             plt.plot(FilteredForceDat)
             plt.plot(locs_min, FilteredForceDat[locs_min], 'x')
             plt.plot(locs_max, FilteredForceDat[locs_max], 'x')
             
             answer = messagebox.askyesno("Question", "Is data clean?")
             
        if answer == False:
            plt.close('all')
            print('Adding file to bad file list')
            badFileList.append(entry)
        
        if answer == True:
            plt.close('all')
            print('Calculating Stiffness and Energy Return')
         
            PercentReturn = []
            Stiffness = []
            Specimen = []
            
            for i, val in enumerate(locs_min[:-1]):
                 try:
                     
                     #i = 0
                     tmpForce = FilteredForceDat[locs_min[i]:locs_max[i]]
                     tmpDispDat = filteredDatDisp[locs_min[i]:locs_max[i]]
                     tmpFDDat = FilteredForceDispDat[locs_min[i]:locs_max[i]]
            
                     tmpMax = np.argmax(tmpForce)
                     tmpMaxFP = np.argmax(tmpDispDat)
                     tmp20 = round(0.2 * tmpMaxFP)
                     tmp40 = round(0.4 * tmpMaxFP)
                     tmp60 = round(0.6 * tmpMaxFP)
            
                     tmpForceRange = tmpForce[tmp20:tmp60]
                     tmpDispRange = tmpDispDat[tmp20:tmp60]
            
                     stiff = np.gradient(tmpForceRange, tmpDispRange)
            
                     rampup = np.trapz(tmpForce, tmpDispDat)
                     rampdown = -1 * np.trapz(FilteredForceDat[locs_max[i]:locs_min[i+1]], filteredDatDisp[locs_max[i]:locs_min[i+1]])
            
                     EnergyLoss = rampup - rampdown
                     PercentReturn.append((1 - (EnergyLoss / rampup)) * 100)
                     Stiffness.append(np.mean(stiff))
                     Specimen.append(int(entry.split('_')[1]))
                 except Exception as e:
                     print(e)
            
           
            outcomes = pd.DataFrame({'Specimen Number':list(Specimen), 'Percent Return': list(PercentReturn), 'Stiffness':list(Stiffness)})
            
            if save_on == 1:
                 if not os.path.exists(outfilePath):
                     outcomes.to_csv(outfilePath, mode='a', header=True, index=False)
                 else:
                     outcomes.to_csv(outfilePath, mode='a', header=False, index=False)
            
