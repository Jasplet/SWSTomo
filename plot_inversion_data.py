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

DATA_DIR = '/Users/ja17375/SWSTomo/data'
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

    paths = pd.read_csv('/Users/ja17375/SWSTomo/Inversions/Joint7Phases.sdb', delim_whitespace=True) 
    n = len(paths)
    pathc = '/Users/ja17375/SWSTomo/Inversions/Dom1160/Joint7Phase/CorrectedPhases'
    
    fig = plt.figure(figsize=(20, 25))
    gs = GridSpec(n, 4, figure=fig, width_ratios=[3,1,3,1])
    for i, path in paths.iterrows():
        filename = f'{path.STAT}_{path.DATE}_{path.TIME}*_{path.PHASE}'
        file_corr = f'{path.STAT}_{path.PHASE}_corrected.pair'
        pair, metadata = read_path_as_pair(filename)
        station = metadata.station
        pairc = sw.load(f'{pathc}/{file_corr}')
        ert = event_relative_time(metadata)
        pair.set_window(path.WBEG - ert, path.WEND - ert)
                   
        ax_in_tr = fig.add_subplot(gs[i, 0])
        ax_in_pm = fig.add_subplot(gs[i, 1])
        ax_out_tr = fig.add_subplot(gs[i, 2])
        ax_out_pm = fig.add_subplot(gs[i, 3])
        pair._ptr(ax_in_tr)
        pair._ppm(ax_in_pm)
        pairc._ptr(ax_out_tr)
        pairc._ppm(ax_out_pm)
        ax_in_tr.text(0.075, 0.9, f'{station} ({path.PHASE})', transform=ax_in_tr.transAxes, fontsize=14)
        if i == 0:
            ax_in_tr.set_title('Uncorrected Traces', fontsize=16)
            ax_out_tr.set_title('Corrected Traces', fontsize=16)
        
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/path_data_test.png')