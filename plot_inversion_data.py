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

    
def add_windows(metadata, axs, yrange):
    '''adds analysis window ranges to ax'''
    windows = [metadata.sac['user0'], metadata.sac['user1'], 
               metadata.sac['user2'], metadata.sac['user3']
               ]
    for ax in axs:
        ax.vlines(windows, ymin = yrange[0], ymax=yrange[1], colors='black', linestyle='solid') 
        ax.set_xlim([windows[0]-5, windows[3]+5])



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
    
    paths = ['HUMO_2008321_170232_SKS',
         'COR_2008321_170232_SKS',
         '116A_2006360_122621_SKKS',
         'K20A_2009003_223342_SKKS',
         'L24A_2009003_194355_SKKS',
         'DAN_2003174_121231_ScS',
         'RDM_2003174_121231_ScS']
    n = len(paths)
    fig = plt.figure(figsize=(14, 18))
    gs = GridSpec(n, 2, figure=fig)
    for i, path in enumerate(paths):
        corr = corrs[corrs.STAT == path.split('_')[0]]
        pair, metadata = read_path(path)
        ert = event_relative_time(metadata)
        ax_in = fig.add_subplot(gs[i, 0])
        ax_out = fig.add_subplot(gs[i,1])
        add_waveforms(pair, ert, ax_in)
        if i == 0:
            ax_in.legend(['BHN', 'BHE'])
        # pairc = remove_splitting(pair, corr)
        # add_waveforms(pairc, ert, ax_out)
        yrange = set_yrange(pair.data(), (ax_in, ax_out))
        add_windows(metadata, (ax_in,ax_out), yrange)
    
    plt.savefig(f'{FIG_DIR}/path_data_test.png')