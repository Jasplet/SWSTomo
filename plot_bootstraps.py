#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 12:26:01 2021

@author: ja17375

Reworked (and hopefulyl more generalised) plotting of boot strapping results
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from EnsembleVisualiser import Ensemble


def plot_bootstrapped_samples(Samples,outfile=None):
    '''
    Function to read and control plotting of bootstrap sampling results.
    Makes a 4-panel plot showing the KDE for alpha, gamma, strength and model misfit
    
    Args:
        Samples - An Ensemble pbject containing the bootstrap samples
    Returns
    -------
    None.

    '''
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8), sharey=False, sharex=False)
    
    # Top left box = alpha param
    add_kde(axs[0,0], Samples, 'alpha')
    axs[0,0].set_xlabel(r'$\alpha$ ($\degree$)')
    axs[0,0].set_ylabel('Prob. Density')
    # Top right box = gamma param
    add_kde(axs[0,1], Samples, 'gamma')
    axs[0,1].set_xlabel(r'$\gamma$ ($\degree$)')
    # Bottom left box = strength
    add_kde(axs[1,0], Samples, 'strength')
    axs[1,0].set_xlabel('strength')
    #Bottom right box = misfit
    add_kde(axs[1,1], Samples, 'misfit')
    axs[1,1].set_xlabel('misfit')
    if outfile:
        plt.savefig(outfile, dpi=250)
    
    plt.show()
    
def add_kde(ax, cand, param):
    '''Make a KDE plot for input values. 
    Which model parameter (alpha, gamma, strength, misfit) is determined by param'''
    if param == 'alpha':
        x = np.linspace(0, 90, 1000)
        lim = [0, 90]
    elif param == 'gamma':
        x = np.linspace(-180, 180, 2000)
        lim=[-180, 180]
    elif param == 'strength':
        x = np.linspace(0, cand.config['strength_max'], 1000)
        lim=[0, cand.config['strength_max']]
    elif param == 'misfit':
        x = np.linspace(0, 1.5, 1000)
        lim=[0, 0.8]
    values = cand.models[param]
    kde = stats.gaussian_kde(values)
    y = kde(x)
    # ymax = np.ceil(y.max())
    ax.plot(x, y, 'k')

    sig_min = values.mean()-values.std()
    sig_max = values.mean()+values.std()
    sig_x = np.linspace(sig_min, sig_max, 200)
    ax.fill_between(sig_x, 0, kde(sig_x), alpha=0.5)
    #ax.vlines(values.mean(), ymin=0, ymax=250,colors='k')
    ax.set_xlim(lim)
    
if __name__ == '__main__':
        path = '/Users/ja17375/Projects/Matisse_Synthetics/ppv1'
        Samples = Ensemble(f'{path}/ideal/Noise01/Bootstrapping',
                           fname='ppv_100_001_Ideal_noise01_bootstraps.out',
                           strmax=0.5)
        plot_bootstrapped_samples(Samples,
                                  outfile=f'{path}/ideal/Noise01/Bootstrapping/ppv_100_001_bootstraps.png')
    