#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 17:11:53 2021

@author: ja17375
Code to geterate bootstrapped MTS pathset XML files 
"""
from pathset import PathSetter
import os
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from EnsembleVisualiser import Ensemble
from plot_bootstraps import add_kde

PATH = '/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/Bootstrapping'

def make_bootstrap_samples(Setter, n):
    '''
    Generate the Pathset file for the ith Bootstrap sample
    '''
    for i in range(1, n+1):
        
        samp_name = f'Sample_{i}'
        os.mkdir(f'{PATH}/Pathsets/{samp_name}')
        Setter.gen_PathSet_XML(stations='random_sample', lower_dom_no='fast_anom',
                            phases=['ScS','SKS','SKKS'],fname=f'{samp_name}/Pathset')

def plot_candidates(param='misfit'):
    if param == 'strength':
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8), sharey=False, sharex=False)
        bestfits = [0.008, 0.484, 0.143, 0.245]
        ymax = [250, 15, 15, 15]
    elif param == 'alpha':   
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8), sharey=True, sharex=True)
        bestfits = [74.491, 71.502, 60.673, 34.404]
        ymax = [0.2, 0.2, 0.2, 0.2]
    elif param == 'gamma':
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8), sharey=True, sharex=True)
        bestfits = [-140.557, -150.188, 38.264, -12.267 ]
        ymax =[0.08, 0.08, 0.08, 0.08]
    elif param == 'misfit':
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8), sharey=True, sharex=True)
        bestfits = [0.254, 0.289, 0.433, 0.328]  
        ymax = [8, 8, 8, 8]
    plt.style.use('default')
    # load candidates
    ellip = Ensemble(f'{PATH}/Results/', strmax=0.02, fname='ellip_bootstraps.out')
    br = Ensemble(f'{PATH}/Results/', strmax=0.5, fname='pv_100_001_bootstraps.out')
    ppv1 = Ensemble(f'{PATH}/Results/', strmax=0.5, fname='ppv_001_100_bootstraps.out')
    ppv2 = Ensemble(f'{PATH}/Results/', strmax=0.5, fname='ppv_010_100_bootstraps.out')
    candidates = [ellip, br, ppv1, ppv2]
    titles = ['Elliptical TI', 'Bridgmanite [001](100)',
              'Post-perovskite [100](001)', 'Post-perovskite [100](010)']
    

    for i, cand in enumerate(candidates):
        cand_param = cand.models[param]
        cand_bestfit = bestfits[i]
        ax = axs.ravel()[i]
        
        add_kde(ax, cand, param)
        ax.vlines(cand_bestfit, ymin=0, ymax=ymax[i], colors='red')
        ax.set_title(titles[i])
        print(titles[i])
        print(f' Sample mean: {cand_param.mean()} +/- {cand_param.std()}')
        ax.text(0.005, ymax[i]*0.85, r'$\mu$ = {:4.0f} $\pm$ {:4.0f}'.format(cand_param.mean(), cand_param.std()))
        ax.text(0.005, ymax[i]*0.8, r'$\chi$ = {:4.3f}'.format(cand_bestfit))
        ax.set_ylim([0, ymax[i]])
    axs.ravel()[0].set_ylabel('Density')
    axs.ravel()[2].set_ylabel('Density')
    axs.ravel()[2].set_xlabel(param.capitalize())
    axs.ravel()[3].set_xlabel(param.capitalize())
                                                            
    plt.tight_layout()
    plt.savefig(f'/Users/ja17375/SWSTomo/Figures/Candidate_models_bootstrapped_{param}.png', dpi=400)
    plt.show()

if __name__ == '__main__':
    
    # outdir = f'{PATH}/Pathsets'
    # Setter = PathSetter(phasefile='/Users/ja17375/SWSTomo/Inversions/HQ_phases_on_fast_anom.sdb',
    # odir=outdir, model='EP_fast_anom_Model.xml')
    
    #make_bootstrap_samples(Setter, 500)
    plot_candidates(param='alpha')
