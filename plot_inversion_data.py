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
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import splitwavepy as sw
import pandas as pd
import numpy as np
from correct_waveforms import event_relative_time

MODEL = 'EllipTI'
#DATA_DIR = f'/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/{MODEL}/wavecorr/corrs'
DATA_DIR = '/Users/ja17375/SWSTomo/data'
FIG_DIR = '/Users/ja17375/SWSTomo/Figures/Waveforms'

def read_path_as_pair(fstem, full_path=False):
    ''' 
    Reads the 3-component waveform (SAC) data 
    '''
    if not full_path:
        st = obspy.read(f'{DATA_DIR}/{fstem}.BH?')        
    else:
        st = obspy.read(f'{fstem}.BH?')
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
    
def predict_tt(evdp, dist, phase='ScS'):
    model = TauPyModel()
    # uses IASP91 by default
    arrivals = model.get_travel_times(source_depth_in_km=evdp,
                                      distance_in_degree=dist,
                                      phase_list=[phase])
    tt = arrivals[0].time
    return tt

def plot_individual_input_waveforms(paths):
    
    for i, path in paths.iterrows():
        fig = plt.figure(figsize=(7,3))
        ax = fig.add_subplot(111)
        filename = f'{path.STAT}_{path.DATE}_*_{path.PHASE}'
        #filename = f'{path.STAT}_{path.DATE}_*corr.00000002'
        pair, metadata = read_path_as_pair(filename)
        station = metadata.station
        ert, evt = event_relative_time(metadata)
        ax.plot(pair.t()+ert, pair.x, label= 'BHN')
        ax.plot(pair.t()+ert, pair.y, label= 'BHE')
        ax.axvline(path.WBEG, linewidth=1, color='k')
        ax.axvline(path.WEND, linewidth=1, color='k')
        datestring = evt.strftime('%Y-%m-%d %H:%M:%S')
        ax.set_title(f'Event {datestring}. Station {path.STAT}. Phase {path.PHASE}')
        tt = predict_tt(path.EVDP, path.DIST, path.PHASE)
        ax.set_xlim([path.WBEG - 20, path.WEND + 20])
        ax.legend(framealpha=0.75, loc=1)     
        ax.set_xlabel('Time (relative to event time) (s)')
        plt.tight_layout(w_pad=1.25)
        fig.savefig(f'{FIG_DIR}/{station}_input_waveform.png')
    
def plot_all_particle_motions(paths, bw=True):

    corr_stem = '/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/'        
    models  = ['Input','EllipTI', 'pv_100_001', 'ppv_001_100',
            'ppv_010_100']    
    fig, axs = plt.subplots(nrows=11, ncols=5, figsize=(12, 20))
    for i, path in paths.iterrows():
        ax_row = axs[i, :]
        for i, model in enumerate(models):
            if model == 'Input':
                ddir = '/Users/ja17375/SWSTomo/data'
                filename = f'{ddir}/{path.STAT}_{path.DATE}_*_{path.PHASE}'    
            else:
                ddir = f'{corr_stem}/{model}/wavecorr/corrs'
                filename = f'{ddir}/{path.STAT}_{path.DATE}_*corr.00000002'
               
            pair, metadata = read_path_as_pair(filename, full_path=True)
            ert = event_relative_time(metadata)[0]
            pair.set_window(path.WBEG - ert, path.WEND - ert)
            ax = ax_row[i]
            if bw:
                x, y = pair.chop().x, pair.chop().y 
                # x is North and y is East Component
                lim = np.abs(pair.data()).max() * 1.1
                ax.plot(y, x, 'k')
                ax.set_aspect('equal')
                ax.set_xlim([-lim, lim])
                ax.set_ylim([-lim, lim])
                
                ax.axes.xaxis.set_ticklabels([])
                ax.axes.yaxis.set_ticklabels([])
                ext = '_BW'
                
            else:
                pair._ppm(ax)
                ext = ''
            ax.set_xlabel(None)
            ax.set_ylabel(None)
            
    
    plt.tight_layout(pad=1.05)
    fig.savefig(f'{FIG_DIR}/ParticleMotions/All_particle_motions{ext}.png',dpi=500)

if __name__ == '__main__':
    
    #~/SWSTomo/Inversions/HQ_phases_NoIRON_fa_anom.sdb 
    paths = pd.read_csv('/Users/ja17375/SWSTomo/Inversions/HQ_phases_on_fast_anom.sdb', delim_whitespace=True) 
    plot_all_particle_motions(paths)
 #   plot_individual_input_waveforms(paths)
#
        
    plt.show()
