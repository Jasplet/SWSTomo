#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 15:45:49 2021

@author: ja17375

Reads Matisse Ensembles files in working dir, finds best fitting models and returns
"""

import pandas as pd
import argparse
import numpy as np
from EnsembleVisualiser import Ensemble, restore_params
#path is the directiory to read from!

def read_best_model(path, candidate, strmax):
    '''Reads BestFitting.out file from Matisee and translates it back from param vectors
    Args:
        path - absolute path to xxx_bestfitting.out file
        canddiate - candidate model used (ellipTI, ppv_100_001, ppv_100_010 or pv_001_100)
        strmax - maximum strength value 
    
    '''
    model = np.loadtxt(f'{path}/{candidate}_bestmodel.out')
    alpha = np.around(restore_params(model[0], 0, 90), decimals=3)
    gamma = np.around(restore_params(model[1], -180, 180), decimals=3)
    strength = np.around(restore_params(model[2], 0, strmax), decimals=3)    
    misfit = np.around( 1/model[-1], decimals=3)
    
    return alpha, gamma, strength, misfit

if __name__ == '__main__':

    bazs = ['ideal', 'real']    
    noises = ['02', '01', '075', '005']
    store = {'candidate': [], 'noise_level':[],'backazimuthal_coverage': [],
             'alpha': [], 'gamma': [], 'strength': [], 'misfit': []}
    for baz in bazs:
        for noise in noises:            
            noisedir = f'Noise{noise}'
            path = f'/Users/ja17375/Projects/Matisse_Synthetics/ppv1/{baz}/{noisedir}/results'
            if noise in ['01', '02']:
                #higher noise levels -n 0.1 and -n 0.2 in sacsplitwave
                n = float(noise)/10
            elif noise == '005':
                n = 0.05
            elif noise == '075':
                n = 0.075
            else:
                raise ValueError()
            for cand in ['ellipTI', 'ppv_100_001', 'ppv_100_010', 'pv_001_100']:    
                if cand == 'ellipTI':
                    strmax = 0.02
                else:
                    strmax = 0.5
                a, g, s, mf = read_best_model(path, cand, strmax)
                store['candidate'].append(cand)
                store['noise_level'].append(n)
                store['backazimuthal_coverage'].append(baz)
                store['alpha'].append(a)
                store['gamma'].append(g)
                store['strength'].append(s)
                store['misfit'].append(mf)
                
    df = pd.DataFrame(store)
    df.to_csv('/Users/ja17375/Projects/Matisse_Synthetics/BestFittingModels.txt', sep=' ', index=False)