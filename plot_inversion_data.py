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
DATA_DIR = '/Users/ja17375/Projects/Epac_fast_anom/data'
FIG_DIR = '/Users/ja17375/Projects/Epac_fast_anom/Figures'

def read_path_as_pair(fstem, full_path=False, rotate=False):
    ''' 
    Reads the 3-component waveform (SAC) data 
    '''
    if not full_path:
        st = obspy.read(f'{DATA_DIR}/{fstem}.BH?')        
    else:
        st = obspy.read(f'{fstem}.BH?')
    
    metadata = st[0].stats
    
    if rotate:
        # rotate to RTZ
        baz = st[0].stats.sac['baz']
        st.rotate('NE->RT', back_azimuth=baz)
        trans = st[1].data
        radial = st[0].data
        pair = sw.Pair(trans, radial, delta=metadata.delta)
    else:
        north = st[1].data 
        east = st[0].data
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

def plot_input_waveforms(paths):
    
    fig, axs = plt.subplots(nrows=6, ncols=2, figsize = (12,14))
    for i, path in paths.iterrows():
        if i <= 5:
            ax = axs[i,0]
        else:
            ax = axs[i-6,1]

        filename = f'{path.STAT}_{path.DATE}_*_{path.PHASE}'
        #filename = f'{path.STAT}_{path.DATE}_*corr.00000002'
        st = obspy.read(f'{DATA_DIR}/{filename}.BH?')  
        metadata = st[0].stats
       #Roate to RTZ
        baz = metadata.sac['baz']
        st.rotate('NE->RT', back_azimuth=baz)
        trans = st.select(channel='BHT')[0].data
        radial = st.select(channel='BHR')[0].data       
        time = np.arange(metadata.npts)*metadata.delta
        ert, evt = event_relative_time(metadata)
        ax.plot(time+ert, trans, label= 'Transverse')
        ax.plot(time+ert, radial, label= 'Radial')
        ax.axvline(path.WBEG, linewidth=1, color='k')
        ax.axvline(path.WEND, linewidth=1, color='k')
        datestring = evt.strftime('%Y-%m-%d %H:%M:%S')
        ax.set_title(f'Event {datestring}. Station {path.STAT}. Phase {path.PHASE}')
        ax.set_xlim([path.WBEG - 20, path.WEND + 20])
        ax.legend(framealpha=0.75, loc=1)     
        ax.set_xlabel('Time (relative to event time) (s)')
        plt.tight_layout(w_pad=1.25)
        fig.savefig(f'{FIG_DIR}/HQ_input_waveform.eps',format='eps',dpi=500)
    

def plot_all_particle_motions(paths, bw=True):

    corr_stem = '/Users/ja17375/Projects/Epac_fast_anom/HQ_data/ScS_fix_test'        
    models  = ['Input','ellipTI', 'pv_100_001', 'ppv_001_100',
            'ppv_010_100']    
    fig, axs = plt.subplots(nrows=11, ncols=5, figsize=(7, 11))
    i = 0
    for ind, path in paths.iterrows():
        ax_row = axs[i, :]
        i+=1
        for j, model in enumerate(models):
            if model == 'Input':
                ddir = '/Users/ja17375/Projects/Epac_fast_anom/data'
                filename = f'{ddir}/{path.STAT}_{path.DATE}_*_{path.PHASE}'    
            else:
                ddir = f'{corr_stem}/{model}/wavecorr/corrs'
                filename = f'{ddir}/{path.STAT}_{path.DATE}_*corr.00000002'
               
            pair, metadata = read_path_as_pair(filename, full_path=True)
            ert = event_relative_time(metadata)[0]
            pair.set_window(path.WBEG - ert, path.WEND - ert)
            ax = ax_row[j]
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
    fig.savefig(f'{FIG_DIR}/ScS_particle_motions{ext}.png',dpi=500)

def plot_all_particle_motions_ls(paths, bw=True):
    ''' 
    Plots all particle motions but in a landscape grid
    '''
    corr_stem = '/Users/ja17375/Projects/Epac_fast_anom/'        
    models  = ['Input','ellipTI', 'pv_100_001', 'ppv_001_100',
            'ppv_010_100']    
    fig, axs = plt.subplots(nrows=5, ncols=11, figsize=(3, 5))
    i = 0
    for ind, path in paths.iterrows():
        ax_col = axs[:, i]
        i +=1
        for j, model in enumerate(models):
            if model == 'Input':
                ddir = '/Users/ja17375/Projects/Epac_fast_anom/data'
                filename = f'{ddir}/{path.STAT}_{path.DATE}_*_{path.PHASE}'    
            else:
                ddir = f'{corr_stem}/{model}/wavecorr/corrs'
                filename = f'{ddir}/{path.STAT}_{path.DATE}_*corr.00000002'
               
            pair, metadata = read_path_as_pair(filename, full_path=True)
            ert = event_relative_time(metadata)[0]
            pair.set_window(path.WBEG - ert, path.WEND - ert)
            ax = ax_col[j]
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
    fig.savefig(f'{FIG_DIR}/All_particle_motions_ls_{ext}.png',dpi=500)

def plot_ScS_waveform_test(paths, bw=True):
    '''Plots the input waveform particle motions and the corrected particle motions of ppv[100](001) with and without our ScS fix'''  
    fig, axs = plt.subplots(nrows=11, ncols=3, figsize=(3,11))
    # First col (Input Data)
    ddirs = ['/Users/ja17375/Projects/Epac_fast_anom/data',
             '/Users/ja17375/Projects/Epac_fast_anom/ppv_001_100/wavecorr/corrs',
             '/Users/ja17375/Projects/Epac_fast_anom/ScS_fix_test/wavecorr/corrs']
    for j, ddir in enumerate(ddirs):
        for i, path in paths.iterrows():
            if j == 0:
                filename = f'{ddir}/{path.STAT}_{path.DATE}_*_{path.PHASE}'    
            else:
                filename = f'{ddir}/{path.STAT}_{path.DATE}_*corr.00000002'
               
            pair, metadata = read_path_as_pair(filename, full_path=True)
            ert = event_relative_time(metadata)[0]
            pair.set_window(path.WBEG - ert, path.WEND - ert)
            ax = axs[i,j]
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
    axs[0,0].set_title('Input Data')
    axs[0,1].set_title('ppv1 old')
    axs[0,2].set_title('ppv1 ScS fix')
    
    plt.tight_layout(pad=1.05)
    fig.savefig(f'{FIG_DIR}/Compare_ScS_fix_results_pm.png',dpi=500)

if __name__ == '__main__':
    
    #~/SWSTomo/Inversions/HQ_phases_NoIRON_fa_anom.sdb 
    paths = pd.read_csv('/Users/ja17375/Projects/Epac_fast_anom/HQ_data/HQ_phases_on_fast_anom.sdb', delim_whitespace=True) 
   # sample_paths = paths[paths.STAT.isin(['COR','K20A','DAN'])]
    plot_all_particle_motions(paths)
  #  plot_input_waveforms(paths)
    #plot_ScS_waveform_test(paths)
        
    plt.show()
