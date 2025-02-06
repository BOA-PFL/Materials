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

import pandas as pd  # For data manipulation and analysis
import os  # For interacting with the operating system

import numpy as np  # For numerical computations
from matplotlib import pyplot as plt  # For plotting data
from scipy.signal import butter, filtfilt  # For signal processing

# Path to the directory containing the CSV files
fPath = 'C:\\Users\\adam.luftglass\\OneDrive - Boa Technology Inc\\General\\Materials Testing\\Fabric Stiffness\\117\\'
fileExt = r".csv"  # File extension of the target files
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]  # List of CSV files in the directory

save_on = 1  # Flag to save the results
outfileName = fPath + 'testResults.csv'  # Output file name for results

fr = 100  # Sampling frequency in Hertz

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
    if entry.split(' ')[0].split('_')[-1] == 'Channels':
        # Read the CSV file
        dat = pd.read_csv(fPath + entry, sep=';', header=0, skiprows=[1])
        dat.Force = dat.Force * 1000  # Convert force to correct units

        normal_cutoff = cutoff / nyq
        b, a = butter(order, normal_cutoff, btype='low', analog=False)

        # Extract force data
        ForceDat = dat.Force

        # Apply low-pass filter to the force data
        FilteredForceDat = butter_lowpass_filter(ForceDat, cutoff, fr, 2)

        # Plot the filtered force data
        fig, ax = plt.subplots(1, 1)
        plt.plot(FilteredForceDat)
        plt.ylabel('Force (N)')
        plt.xlabel('Centiseconds (cs)')

        # Prepare results for saving
        arrayValues = [[int(entry.split('_')[1]), np.min(FilteredForceDat)]]
        col_vals = ['SpecimenNumber', 'Peak Load']
        outcomes = pd.DataFrame(data=arrayValues, columns=col_vals)

        # Save results to the output file
        if save_on == 1:
            if os.path.exists(outfileName) == False:
                outcomes.to_csv(outfileName, mode='a', header=True, index=False)
            else:
                outcomes.to_csv(outfileName, mode='a', header=False, index=False)
