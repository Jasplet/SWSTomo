#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 11:43:56 2021

@author: ja17375

Script to make KDE pairwise plots
"""

from EnsembleVisualiser import Ensemble


# Free Beta test
# Yeah so this didnt work too well - we dont have the coverage for FreeBeta here [7/12/21]
# ppv1_fb = Ensemble('/Users/ja17375/Projects/Epac_fast_anom/FreeBeta/ppv',fname='ppv_010_100_ensemble.out', strmax=0.5, dims='all')
# ppv1_fb.pairwise_plot(fname='ppv1_fb_wieghted_KDE', params=['alpha', 'beta', 'gamma','strength'])

path = '/Users/ja17375/Projects/Matisse_Synthetics/ppv1/ideal/Noise01/results'

for cand in ['ellipTI', 'pv_001_100', 'ppv_100_001', 'ppv_100_010']:
    fname = f'{cand}_ensemble.out'
    if cand == 'ellipTI':    
        #Handle that ellipTI has a different max strength
        Ens = Ensemble(f'{path}',fname=fname, dims='ags')
    else:
        Ens = Ensemble(f'{path}',fname=fname, dims='ags', strmax=0.5)
        
    Ens.pairwise_plot(params=['alpha','gamma','strength'],fname=f'{cand}_synthetics_n01_kde')
    

