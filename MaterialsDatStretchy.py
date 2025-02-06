

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
from scipy.signal import butter, filtfilt
import scipy.signal as sig
from scipy import integrate
from tkinter import messagebox

# Path to the directory containing the CSV files
fPath = 'C:\\Users\\adam.luftglass\\OneDrive - Boa Technology Inc\\General\\Materials Testing\\Swatch Creation\\92624\\'
fileExt = r".csv"  # File extension of the target files
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]  # List of CSV files in the directory

check_data =0  # Flag to check the data
save_on = 1  # Flag to save the results
outfileName = fPath + 'testResults.csv'  # Output file name for results
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
    #entry = entries[2]
    if entry.split(' ')[0].split('_')[-1] == 'Channels':
        # Read the CSV file
        dat = pd.read_csv(fPath + entry, sep=';', header=0, skiprows=[1])
        dat.Force = dat.Force * 1000  # Convert force to correct units

        normal_cutoff = cutoff / nyq
        b, a = butter(order, normal_cutoff, btype='low', analog=False)

        # Extract force and displacement data
        ForceDat = dat.Force
        DispDat = dat.Displacement
        ForceDisp = ForceDat / DispDat

        # Apply low-pass filter to the data
        FilteredForceDat = butter_lowpass_filter(ForceDat, cutoff, fr, 2)
        filteredDatDisp = butter_lowpass_filter(DispDat, cutoff, fr, 2)
        FilteredForceDispDat = FilteredForceDat / filteredDatDisp

        # Identify local minima in the force data
        locs, _ = sig.find_peaks(-1 * FilteredForceDat, distance=200)
        pks = np.array(FilteredForceDat[locs])
        adaptiveThresh = pks[np.where((pks > 4) & (pks < 10))].mean()
        padding = 0.8
        FilteredForceDat[FilteredForceDat < adaptiveThresh + padding] = adaptiveThresh

        minima = []
        stops = []
        for i, val in enumerate(FilteredForceDat):
            if i == (len(FilteredForceDat) - 1):
                pass
            elif (val == adaptiveThresh) & (FilteredForceDat[i + 1] > adaptiveThresh):
                minima.append(i)
            elif (val == adaptiveThresh) & (FilteredForceDat[i - 1] > adaptiveThresh):
                stops.append(i)

        # Handle files with no detected minima
        if len(minima) == 0:
            print('No minima detected, adding to badFileList; check file ' + entry)
            badFileList.append(entry)
        else:
            n = 110
            del minima[:n]
            del stops[:n]
            try:
                if stops[0] < minima[0]:
                    stops.pop(0)
                if minima[-1] > stops[-1]:
                    minima.pop(-1)
            except Exception as e:
                print(e)
                answer = False
                badFileList.append(entry)

            if check_data == 1:
                fig, ax = plt.subplots(1, 1)
                for i, val in enumerate(minima):
                    try:
                        if i == len(minima) - 1:
                            ax.text(1.5, 35, 'Up', fontsize=12, color='blue')
                            ax.text(1.5, 30, 'Down', fontsize=12, color='orange')
                        else:
                            start = minima[i]
                            mid = minima[i] + round((stops[i] - minima[i]) / 2)
                            end = stops[i]

                            ax.plot(filteredDatDisp[start:mid], FilteredForceDat[start:mid], color='blue', label='up')
                            ax.plot(filteredDatDisp[mid:end], FilteredForceDat[mid:end], color='orange', label='down')
                    except Exception as e:
                        print(e)
                        answer = False
                        badFileList.append(entry)
                answer = messagebox.askyesno("Question", "Is data clean?")

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
            stiff = []

            # Calculate energy loss, percent return, and stiffness for each cycle
            for i, val in enumerate(minima):
                try:
                    tmpForce = FilteredForceDat[minima[i]:stops[i]]
                    tmpDispDat = filteredDatDisp[minima[i]:stops[i]]
                    tmpFDDat = FilteredForceDispDat[minima[i]:stops[i]]

                    tmpMax = np.argmax(tmpForce)
                    tmpMaxFP = np.argmax(tmpDispDat)
                    tmp20 = round(0.2 * tmpMaxFP)
                    tmp40 = round(0.4 * tmpMaxFP)
                    tmp60 = round(0.6 * tmpMaxFP)

                    tmpForceRange = tmpForce[tmp20:tmp60]
                    tmpDispRange = tmpDispDat[tmp20:tmp60]

                    stiff = np.gradient(tmpForceRange, tmpDispRange)

                    rampup = integrate.trapz(tmpForce[0:tmpMax], tmpDispDat[0:tmpMax])
                    rampdown = -1 * integrate.trapz(tmpForce[tmpMax:-1], tmpDispDat[tmpMax:-1])

                    EnergyLoss = rampup - rampdown
                    PercentReturn.append((1 - (EnergyLoss / rampup)) * 100)
                    Stiffness.append((tmpForce[tmp40] - tmpForce[tmp20]) / (tmpDispDat[tmp40] - tmpDispDat[tmp20]))

                    plt.plot(FilteredForceDat[tmp20:tmp60])
                except Exception as e:
                    print(e)

            # Prepare results for saving
            arrayValues = [[int(entry.split('_')[1]), np.mean(PercentReturn), np.mean(Stiffness), np.mean(stiff)]]
            col_vals = ['SpecimenNumber', 'PercentReturn', 'Stiffness OG', 'Stiffness New']

            outcomes = pd.DataFrame(
                data=arrayValues,
                columns=col_vals)

            # Save results to the output file
            if save_on == 1:
                if os.path.exists(outfileName) == False:
                    outcomes.to_csv(outfileName, mode='a', header=True, index=False)
                else:
                    outcomes.to_csv(outfileName, mode='a', header=False, index=False)
