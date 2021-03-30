#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 11:33:08 2021

@author: ja17375
"""
import obspy 
from obspy import UTCDateTime
import splitwavepy as sw
import numpy as np
import pandas as pd
import pickle
DATA_DIR = '/Users/ja17375/SWSTomo/data'

class Corrector:
    '''Corrects expected splitting for a single shear-wave'''
    def __init__(self, path, corrections, phase, window):
        
        self.pair, self.metadata = read_path_as_pair(path)
        self.phase = phase #S-wave phase type, informs corrections to be made
        self.corrections = corrections
        self.correction_applied = False
        self.window = window - event_relative_time(self.metadata)
        self.pair.set_window(self.window[0], self.window[1])
        print(f'Apply corrections to event {path}')

    def apply_corrections(self):
        '''Splitting Corrections are applied in reverse propagation order (reciver, Lowermost mantle, Source) '''
        # Remove Station splitting
        self.pair.unsplit(self.corrections['r_fast'], self.corrections['r_lag'])
        # Remove Inverted D'' 
        self.pair.unsplit(self.corrections['d_fast'], self.corrections['d_lag'])
        # Remove Source side splitting ( if phase is ScS and SKS/SKKS do not need source side corrections )
        if self.phase == 'ScS':      
            self.pair.unsplit(self.corrections['s_fast'], self.corrections['s_lag'])
        self.correction_applied = True
        
    def measure_splitting(self):
        '''Measure splitting of data (using splitwavepy's implementation of eigenvalue minimisation)
    
        '''
        eigm = sw.EigenM(self.pair, lags=(4.,))
        
        lam2_norm = (eigm.lam2 / eigm.lam1).min()
        if self.correction_applied is True:
            print('Meausured residual (post-corr) splitting is')
            self.residual_eigm = eigm
            self.residual_lam2norm = lam2_norm
        else:
            print('Measured apparent (pre-corr) splitting is')
            self.apparent_eigm = eigm
            self.lam2norm = lam2_norm
        
        print(f'fast = {eigm.fast} +/- {eigm.dfast}')
        print(f'dt = {eigm.lag} +/- {eigm.dlag}')
        print(f' lam2/lam1 = {lam2_norm:4.3f}')
        eigm = 0 
    
    def save_correction(self, filename):
        if self.correction_applied is True:
            outfile = f'/Users/ja17375/SWSTomo/Inversions/Dom1160/Joint7Phase/CorrectedPhases/{filename}'
            with open(outfile, 'wb') as f:
                pickle.dump(self.pair,f)    
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

def event_relative_time(st_stats):
    '''
    Use obspy's UTSDateTime to work out the time since source event
    
    '''
    start = st_stats.starttime 
    sacstat = st_stats.sac
    
    startdate = UTCDateTime(start)
    eventtime = UTCDateTime(year=sacstat['nzyear'], julday=sacstat['nzjday'], hour=sacstat['nzhour'],
                            minute=sacstat['nzmin'], second=sacstat['nzsec'], microsecond=sacstat['nzmsec'])
    rel_start = startdate - eventtime
    return rel_start, eventtime

def save_pair(Pair, filename):
    '''
    Dumps pair to pickle object. Borrowed from SplitWavePy as pip install version appears to have it missing
    '''


if __name__ == '__main__':
   
    paths = pd.read_csv('/Users/ja17375/SWSTomo/Inversions/Joint7Phases.sdb', delim_whitespace=True) 
    corrs = pd.read_csv('/Users/ja17375/SWSTomo/Inversions/Dom1160/Predicted_Splitting.txt',
                        delim_whitespace=True)

        
    misfit = 0 
    for i, path in paths.iterrows():
        filename = f'{path.STAT}_{path.DATE}_{path.TIME}*_{path.PHASE}'
        outfile = f'{path.STAT}_{path.PHASE}_corrected.pair'
        icorr = corrs[corrs.STAT == path.STAT]
        phase = path.PHASE
        window = np.array([path.WBEG, path.WEND])
        if phase == 'ScS':
            corr = {'r_fast':icorr.FAST_R.values[0], 'r_lag': icorr.TLAG_R.values[0],
                    'd_fast':icorr.FAST_D.values[0], 'd_lag': icorr.TLAG_D.values[0],
                    's_fast':icorr.FAST_S.values[0], 's_lag': icorr.TLAG_S.values[0]}
        else:
            corr = {'r_fast':icorr.FAST_R.values[0], 'r_lag': icorr.TLAG_R.values[0],
                    'd_fast':icorr.FAST_D.values[0], 'd_lag': icorr.TLAG_D.values[0]}
        C = Corrector(filename, corr, phase, window)
        C.apply_corrections()
        C.measure_splitting()
        C.save_correction(outfile)
        misfit += C.residual_lam2norm

    print(f'Misfit (sum(lam2) is {misfit}')
