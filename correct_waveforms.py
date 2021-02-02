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
DATA_DIR = '/Users/ja17375/SWSTomo/data'

class Corrector:
    '''Corrects expected splitting for a single shear-wave'''
    def __init__(self, path, corrections, phase):
        
        self.pair, self.metadata = read_path_as_pair(path)
        self.phase = phase #S-wave phase type, informs corrections to be made
        self.corrections = corrections
        self.correction_applied = False
        self.initial_splitting = self.measure_splitting()
        window = self.retrieve_window()
        self.pair.set_window(window[0], window[1])

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
        
    def retrieve_window(self):
        rel_time_corr = event_relative_time(self.metadata)
        wind_start = self.metadata.sac['user0'] - rel_time_corr
        wind_end = self.metadata.sac['user3'] - rel_time_corr     
        return (wind_start, wind_end)
        
    def measure_splitting(self):
        '''Measure splitting of data (using splitwavepy's implementation of eigenvalue minimisation)
    
        '''
        eigm = sw.EigenM(self.pair, lags=(4.,))
        
        lam2_norm = (eigm.lam2 / eigm.lam1).min()
        print('Meausured splitting is')
        print(f'fast = {eigm.fast} +/- {eigm.dfast}')
        print(f'dt = {eigm.lag} +/- {eigm.dlag}')
        print(f' lam2/lam1 = {lam2_norm:4.3f}')
        if self.correction_applied is True:
            self.residual_eigm = eigm
            self.residual_lam2norm = lam2_norm
        else:
            self.apparent_eigm = eigm
            self.lam2norm = lam2_norm
    
    def save_correction(self, filename):
        if self.correction_applied is True:
            self.pair.save(f'/Users/ja17375/SWSTomo/Inversions/Dom1160/Joint7Phase/CorrectedPhases{filename}')
    
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
    return rel_start

if __name__ == '__main__':
    paths = ['HUMO_2008321_170232_SKS']
      # 'COR_2008321_170232_SKS',
      # '116A_2006360_122621_SKKS',
      # 'K20A_2009003_223342_SKKS',
      # 'L24A_2009003_194355_SKKS',
      # 'DAN_2003174_121231_ScS',
      # 'RDM_2003174_121231_ScS']
    
    corrs = pd.read_csv('/Users/ja17375/SWSTomo/Inversions/Dom1160/Predicted_Splitting.txt',
                        delim_whitespace=True)
    for i, path in enumerate(paths):
        icorr = corrs[corrs.STAT == path.split('_')[0]]
        phase = path.split('_')[-1]
        if phase == 'ScS':
            corr = {'r_fast':icorr.FAST_R.values[0], 'r_lag': icorr.TLAG_R.values[0],
                    'd_fast':icorr.FAST_D.values[0], 'd_lag': icorr.TLAG_D.values[0],
                    's_fast':icorr.FAST_S.values[0], 's_lag': icorr.TLAG_S.values[0]}
        else:
            corr = {'r_fast':icorr.FAST_R.values[0], 'r_lag': icorr.TLAG_R.values[0],
                    'd_fast':icorr.FAST_D.values[0], 'd_lag': icorr.TLAG_D.values[0]}
        C = Corrector(path, corr, phase)
        C.apply_corrections()
        C.measure_splitting()
        C.apparent_eigm.plot()
        C.residual_eigm.plot()
        