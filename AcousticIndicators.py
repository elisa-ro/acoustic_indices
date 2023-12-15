"""
My Script

This script calculates relevant acoustic ecoindices for nature sounds. 
Namely, the ACI, ADI, Bioacoustic Index, Acoustic Evenness Index, Temporal 
Entropy, Temporal Median and NDSI. More information on each of the indexes
can be found in the report: 
https://github.com/elisa-ro/acoustic_indices/blob/main/Report.ipynb

Usage:
    python3 AcousticIndicators.py directoryName

Example usage: 
    python3 AcousticIndicators.py C:\\Users\\elisa\\Downloads\\sounds\\ 

Parameters:
    --param1 (str): This is the directory name where the .wav files exist. 
    This is a required parameter.

Result:
    --result (file): A csv file in the abovementioned directory is created. 
    Each row represents an acoustic file (.wav) in that directory, with the 
    calculated indices.
"""

from IPython.display import Audio, display

import os
import csv

import maad
from maad import sound
from maad import features

import numpy as np
import pydub
import soundfile
from pydub import AudioSegment
import pprint

from scipy.signal import spectrogram
import scipy.io.wavfile as wavfile
from scipy import signal
import matplotlib.pyplot as plt

from acoustics import Signal
import ipynb
import sys
_MIN_ = sys.float_info.min

def ADI(path, printt = False):
    #import preprocess as cmi
    s, fs = maad.sound.load(path)
    Sxx, tn, fn, ext = maad.sound.spectrogram (s, fs, mode='amplitude')  
    ADI  = acousticDiversity(Sxx,fn)
    if printt:
        print('ADI : %2.2f ' %ADI)
    return ADI

def acousticDiversity(Sxx, fn, fmin=0, fmax=20000, bin_step=1000, 
                            dB_threshold=-50, index="shannon"):
    N = np.floor((fmax-fmin)/bin_step)
    # convert into dB and normalization by the max
    Sxx_dB = maad.util.amplitude2dB(Sxx/(np.amax(Sxx))) 
    # Score for each frequency in the frequency bandwith
    s_sum = []
    for ii in np.arange(0,N):
        f0 = int(fmin+bin_step*(ii))
        f1 = int(f0+bin_step)
        s,_ = maad.features.alpha_indices._score(Sxx_dB[maad.util.index_bw(fn,(f0,f1)),:], threshold=dB_threshold, axis=0)
        s_sum.append(np.mean(s))
    
    s = np.asarray(s_sum)
    # Entropy
    if index =="shannon":
        ADI = shannonEntropy(s)
    elif index == "simpson":
        s = s/sum(s)
        s = s**2
        ADI = 1-sum(s)
    elif index == "invsimpson":
        s = s/sum(s)
        s = s**2
        ADI = 1/sum(s)   
    
    return ADI

def shannonEntropy(datain, axis=0):
    """
    Shannon Entropy
    
    Parameters
    ----------
    datain : ndarray of floats
        Vector or matrix containing the data
    
    axis : integer, optional, default is 0
        entropy is calculated along this axis.
    Returns
    -------    
    Hs : ndarray of floats
        Vector or matrix of Shannon Entropy
    """
    # length of datain along axis
    n = datain.shape[axis]
    Hs = entropy(datain, axis=axis) * np.log(n)
    return Hs


def entropy(x, axis=0):
    """
    Compute the entropy of a vector (waveform) or matrix (spectrogram).
    
    Parameters
    ----------
    x : ndarray of floats
        x is a vector (1d) or a matrix (2d)
    axis : int, optional, default is 0
        select the axis where the entropy is computed
        if x is a vector, axis=0
        if x is a 2d ndarray, axis=0 => rows, axis=1 => columns
                
    Returns
    -------
    H : float or ndarray of floats
        entropy of x
        
    """
    # force x to be ndarray
    x = np.asarray(x)
    
    if isinstance(x, (np.ndarray)) == True:
        if x.ndim > axis:
            if x.shape[axis] == 0: 
                print ("WARNING: x is empty") 
                H = None 
            elif x.shape[axis] == 1:
                print("entered here because shape is 1")
                H = 0 # null entropy
            elif x.any() == 0: # test if there are only zeros
                if x.ndim == 1 : # case vector
                    H = 1 # entropy = 1
                else : # case matrix
                    if axis == 0 : H = np.ones(x.shape[1]) # entropy = 1
                    if axis == 1 : H = np.ones(x.shape[0]) # entropy = 1
            else:
                # if datain contains negative values -> rescale the signal between 
                # between posSitive values (for example (0,1))
                if np.min(x)<0:
                    x = linear_scale(x,minval=0,maxval=1)
                # length of datain along axis
                n = x.shape[axis]
                # Tranform the signal into a Probability mass function (pmf)
                # Sum(pmf) = 1
                if axis == 0 :
                    pmf = x/np.sum(x,axis)
                elif axis == 1 :                     
                    pmf = (x.transpose()/np.sum(x,axis)).transpose()
                pmf[pmf==0] = _MIN_
                #normalized by the length : H=>[0,1]
                H = -np.sum(pmf*np.log(pmf),axis)/np.log(n)
        else:
            print ("WARNING :axis is greater than the dimension of the array")    
            H = None 
    else:
        print ("WARNING: x must be ndarray")   
        H = None 

    return H

def compute(path, start = None, end = None, printt = False):
    if printt:
        s = Signal.from_wav(path)
        s.plot_spectrogram()
    s, fs = maad.sound.load(path)

    #acoustic complexity index
    Sxx, tn, fn, ext = maad.sound.spectrogram (s, fs, mode='amplitude')  
    _, _ , ACI  = maad.features.acoustic_complexity_index(Sxx)
    if printt:
        print('ACI : %2.0f ' %ACI)
    
    #acoustic diversity index
    adi = ADI(path, printt)

    #biodiveristy index
    BI = maad.features.bioacoustics_index(Sxx,fn)
    
    #acoustic evenness index
    AEI  = maad.features.acoustic_eveness_index(Sxx,fn)
    
    #temporal entropy
    TE = features.temporal_entropy(s)
    
    #temporal median
    TM = features.temporal_median(s)
    
    if printt:
        print('BI Soundecology : %2.2f ' %BI)
        print('AEI : %2.2f ' %AEI)
        print('Temporal entropy: %2.2f ' %TE)
        print('Temporal median: %2.2f '%TM)    

    #NDSI
    Sxx_power, tn, fn, ext = maad.sound.spectrogram (s, fs)  
    NDSI, ratioBA, antroPh, bioPh  = maad.features.soundscape_index(Sxx_power,fn,R_compatible=None)
    if printt:
        print('NDSI MAAD/Seewave: %2.2f ' %NDSI)
    return ACI, adi, BI, AEI, TE, TM, NDSI

def main():
    header = ['name', 'ACI', 'ADI', 'BI','AEI','temporal entropy','temporal median','NDSI']
    data = []

    directory = sys.argv[1]

    for filename in os.scandir(directory):
        print(filename.path)
        # checking if it is a file
        if '.WAV' in filename.path or '.wav' in filename.path:
            aci, adi, bi, aei, te, tm, ndsi = compute(filename.path)
            data.append([filename.path, aci, adi, bi, aei, te, tm, ndsi])

    path = directory + '/bioacoustic_indices.csv'
    with open(path, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)
        # write multiple rows
        writer.writerows(data)

if __name__ == "__main__":
    main()

 