#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:59:09 2021

@author: ja17375
"""
import pandas as pd 
from calc_aoi import slw2aoi, get_rayparam
import numpy as np
date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}

def make_synth_sdb(sdb_in):
    '''
    Reads a SDB of real data, chops off splitting measurements and adds params needed to make 
    splitting predictions needed for synthetic tests
    '''
    sdb = pd.read_csv(sdb_in, converters=date_time_convert, delim_whitespace=True)
    inc = np.zeros((len(sdb), 1))
    dist = np.zeros((len(sdb), 1))
    for i, row in sdb.iterrows():
        inc[i], dist[i] = calc_incl_dist(row.EVDP, row.DIST, row.PHASE)
    
    sdb['INCL'] = inc
    sdb['LENGTH'] = dist
    sdb_out = sdb[['STAT', 'DATE', 'TIME', 'PHASE', 'EVLA', 'EVLO', 'EVDP', 'STLA', 'STLO',
       'DIST', 'BAZ', 'AZI', 'INCL', 'LENGTH']]
    return sdb_out
    
def calc_incl_dist(evdp, gcarc, phase):
    '''calculates inclination of phase and path length in D'' using calc_aoi routines'''
    rp = get_rayparam(evdp, gcarc, phase)
    aoi = slw2aoi(2890., rp, 's') # 2890. is assumed depth to CMB in km
    dom_h = 250. # [km] assumed height of D'' (and therefore our domain)
    
    dist = dom_h / np.cos(np.deg2rad(aoi))
    inc = 90 - aoi  # Inclination = 90 - aoi
    if phase == 'ScS':
        return inc, dist*2
    else:
        return inc, dist

if __name__ == '__main__':
    sdb =  make_synth_sdb('/Users/ja17375/SWSTomo/Inversions/HQ_phases_on_fast_anom.sdb')
    sdb.to_csv('/Users/ja17375/SWSTomo/Inversions/Synth_base.sdb', sep=' ',index=False)