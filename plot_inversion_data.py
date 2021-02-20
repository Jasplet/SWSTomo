#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:27:43 2021

Contains functions to make plots of ScS, SKS and SKKS waveforms prior to inversion
and after corrections for the predicted splitting have been made. 
i.e if obs(x,t) = ani*u(x,t) we are plotting obs and u (assuming our model of ani is correct)

@author: ja17375
"""
import obspy 
from obspy import UTCDateTime
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import splitwavepy as sw
import pandas as pd
import numpy as np
from correct_waveforms import event_relative_time

DATA_DIR = '/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/ppv_model/wavecorr/corrs'
FIG_DIR = '/Users/ja17375/SWSTomo/Figures'

def read_path_as_pair(fstem):
    ''' 
    Reads the 3-component waveform (SAC) data 
    '''
    st = obspy.read(f'{DATA_DIR}/{fstem}.BH?')        
    north = st[1].data 
    east = st[0].data
    metadata = st[0].stats
    pair = sw.Pair(north, east, delta=metadata.delta)      
    return pair, metadata



def add_waveforms(pair, ert, ax):
    '''function to draw the raw input waveforms'''
    
    time = pair.t() + ert
    r = pair.data()[0]
    t = pair.data()[1]
    ax.plot(time, r, "-", label='BHN')
    ax.plot(time, t, "-", label='BHE')

def set_yrange(ydata, axs):
    ymin = rounddown(ydata.min())
    ymax = roundup(ydata.max())
    yrange = [ymin - abs(ymin)*0.1, 
              ymax + abs(ymax)*0.1]
    for ax in axs:
        ax.set_ylim(yrange)
    
    return yrange

def roundup(x):
    return np.ceil(x /100) * 100

def rounddown(x):
    return np.floor(x /100) * 100

if __name__ == '__main__':

    paths = pd.read_csv('/Users/ja17375/SWSTomo/Inversions/HQ_phases_on_fast_anom.sdb', delim_whitespace=True) 
    n = len(paths) 
    fig, axes = plt.subplots(nrows=n, ncols=2,
                               gridspec_kw={'width_ratios': [3, 1]},
                               figsize=(10, 20), constrained_layout=True)
    
    for i, path in paths.iterrows():
        #filename = f'{path.STAT}_{path.DATE}_*_{path.PHASE}'
        filename = f'{path.STAT}_{path.DATE}_*corr.00000002'
        pair, metadata = read_path_as_pair(filename)
        station = metadata.station
        ert = event_relative_time(metadata)
        pair.set_window(path.WBEG - ert, path.WEND - ert)
        ax_in_tr = axes[i, 0]
        ax_in_pm = axes[i, 1]
        pair._ptr(ax_in_tr)
        pair._ppm(ax_in_pm)
        ax_in_tr.set_title(f'Station {path.STAT}, {path.PHASE} recorded {path.DATE} {path.TIME} ')
        ax_in_tr.set_xlim([25, 100])
    
    print(ert)
    fig.savefig(f'{FIG_DIR}/ppv_corrected_waveforms.png')
