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

PATH = '/Users/ja17375/Projects/Epac_fast_anom/Bootstrapping'

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
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(7,7), sharey=False, sharex=False)
        bestfits = [0.006, 0.392, 0.277, 0.252]
        ymax = [250, 15, 15, 15]
    elif param == 'alpha':   
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(7,7), sharey=True, sharex=True)
        bestfits = [81.6, 84.4, 27.2, 28.7]
        ymax = [0.2, 0.2, 0.2, 0.2]
    elif param == 'gamma':
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(7,7), sharey=True, sharex=True)
        bestfits = [-140.6, -144.1, -31.0, -14.6]
        ymax =[0.08, 0.08, 0.08, 0.08]
    elif param == 'misfit':
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(7,7), sharey=True, sharex=True)
        bestfits = [0.302, 0.361, 0.306, 0.311]  
        ymax = [8, 8, 8, 8]
    plt.style.use('default')
    # load candidates
    ellip = Ensemble(f'{PATH}/Results/ScS_fix', strmax=0.02, fname='ellip_bootstraps.out')
    br = Ensemble(f'{PATH}/Results/ScS_fix', strmax=0.5, fname='pv_100_001_bootstraps.out')
    ppv1 = Ensemble(f'{PATH}/Results/ScS_fix', strmax=0.5, fname='ppv_001_100_bootstraps.out')
    ppv2 = Ensemble(f'{PATH}/Results/ScS_fix', strmax=0.5, fname='ppv_010_100_bootstraps.out')
    candidates = [ellip, br, ppv1, ppv2]
    titles = ['Elliptical TI', 'Bridgmanite [001](100)',
              'Post-perovskite [100](001)', 'Post-perovskite [100](010)']
    

    for i, cand in enumerate(candidates):
        cand_param = cand.models[param]
        cand_bestfit = bestfits[i]
        ax = axs.ravel()[i]
        
        add_kde(ax, cand, param, bestfits[i])
        ax.axvline(cand_bestfit, color='red')
        ax.set_title(titles[i])
        print(titles[i])
        print(f' Sample mean: {cand_param.mean()} +/- {cand_param.std()}')
            
        ax.text(cand_bestfit*0.25, ymax[i]*0.85, r'$\mu$ = {:4.3f} $\pm$ {:4.3f}'.format(cand_param.mean(), cand_param.std()))
        ax.text(cand_bestfit*0.25, ymax[i]*0.78, r'$\chi$ = {:4.3f}'.format(cand_bestfit))
        ax.set_ylim([0, ymax[i]])
    axs.ravel()[0].set_ylabel('Density')
    axs.ravel()[2].set_ylabel('Density')
    axs.ravel()[2].set_xlabel(param.capitalize())
    axs.ravel()[3].set_xlabel(param.capitalize())
                                                            
    plt.tight_layout()
    plt.savefig(f'/Users/ja17375/Projects/Epac_fast_anom/Candidate_models_bootstrapped_{param}_ScS_fix.png', dpi=400)
    plt.show()

def plot_all_cands_all_params(save=False):
    '''
    Makes an unholy 4x4 plot of bootstrapping results 
    '''    
    ellip = Ensemble(f'{PATH}/Results/ScS_fix', strmax=0.02, fname='ellip_bootstraps.out')
    ellip_best = [74, -141, 0.008, 0.254]
    br = Ensemble(f'{PATH}/Results/ScS_fix', strmax=0.5, fname='pv_100_001_bootstraps.out')
    br_best = [72, -150, 0.48, 0.289]
    ppv1 = Ensemble(f'{PATH}/Results/ScS_fix', strmax=0.5, fname='ppv_001_100_bootstraps.out')
    ppv1_best = [61, 38, 0.14, 0.433]
    ppv2 = Ensemble(f'{PATH}/Results/ScS_fix', strmax=0.5, fname='ppv_010_100_bootstraps.out')
    ppv2_best = [34, -12, 0.25, 0.328]
    candidates = [ellip, br, ppv2, ppv1]
    c_fits = [ellip_best, br_best, ppv2_best, ppv1_best]
    fig, axs = plt.subplots(nrows=4, ncols=4, figsize = (8.5,8.5))
    titles = ['Elliptical TI', 'Bridgmanite [001](100)',
              'Post-perovskite [100](001)', 'Post-perovskite [100](010)']
#   Row 1 (Alpha)
    for i, cand in enumerate(candidates):
        for j, param in enumerate(['alpha', 'gamma', 'strength', 'misfit']):
            cand_param = cand.models[param]
            bestfit = c_fits[i][j]
            ax_ij = axs[i,j]
            add_kde(ax_ij, cand, param, bestfit)
            print(f' Sample mean: {cand_param.mean()} +/- {cand_param.std()}')
            # ax_ij.text(0.005, ymax[i]*0.85, r'$\mu$ = {:4.0f} $\pm$ {:4.0f}'.format(cand_param.mean(), cand_param.std()))
            # ax_ij.text(0.005, ymax[i]*0.8, r'$\chi$ = {:4.3f}'.format(cand_bestfit))
            if j == 0 :
                ax_ij.set_ylabel('Density')
            
    axs[3,0].set_xlabel(r'$\alpha (\degree)$')
    axs[3,1].set_xlabel(r'$\gamma (\degree)$')
    axs[3,2].set_xlabel('Strength')
    axs[3,3].set_xlabel('Misfit')
    plt.tight_layout()
    if save:
        plt.savefig('/Users/ja17375/Projects/Epac_fast_anom/All_bootstrapped_params_ScSfix.eps', dpi=400)
    
    plt.show()
if __name__ == '__main__':
    
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=10)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    # outdir = f'{PATH}/Pathsets'
    # Setter = PathSetter(phasefile='/Users/ja17375/Projects/Epac_fast_anom/HQ_phases_on_fast_anom.sdb',
    # odir=outdir, model='/Users/ja17375/Projects/Epac_fast_anom/Models/EP_fast_anom_Model.xml')
    
    # make_bootstrap_samples(Setter, 500)
    # plot_all_cands_all_params()
    for param in ['alpha', 'gamma', 'strength', 'misfit']:
        plot_candidates(param)
