#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 11:57:46 2021

@author: ja17375
"""
import matplotlib.pyplot as plt
FIG_DIR = '/Users/ja17375/SWSTomo/Figures'

def scatter_as(Ensemble, n=None, fname=None):
    '''
    scatter plot for 2-D ensemble with alpha and strength as free paramaters

    Parameters
    ----------
    models : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
        
    if n:
        #if n is given plot n-best models
        models_sort = Ensemble.models.sort_values(by=['misfit'])
        nmodels = models_sort.iloc[0:n]
        # sort again so lowest msfits are at the end and therefore drawn on top!
        models = nmodels.sort_values(by=['misfit'], ascending=False)
    else:
        # if no n plot all models
        models = Ensemble.models.sort_values(by=['misfit'])

    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)
    C = ax.scatter(x=models.strength, y=models.alpha, c=models.misfit, cmap='viridis_r')
    ax.set_title(f'2-D inversion results, {len(models)} models')
    ax.set_xlabel('Strength')
    ax.set_ylabel(r'Alpha $(\degree$)')
    fig.colorbar(C)
    
    if fname:
        plt.savefig(f'{FIG_DIR}/fname', format='png')
    else:
        plt.show()
        
def scatter_gs(Ensemble, n=None, fname=None):
    '''
    scatter plot for 2-D ensemble with gamma and strength as free paramaters
    '''
        
    if n:
        #if n is given plot n-best models
        models_sort = Ensemble.models.sort_values(by=['misfit'])
        nmodels = models_sort.iloc[0:n]
        # sort again so lowest msfits are at the end and therefore drawn on top!
        models = nmodels.sort_values(by=['misfit'], ascending=False)
    else:
        # if no n plot all models
        models = Ensemble.models.sort_values(by=['misfit'])
        n = len(models)

    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)
    C = ax.scatter(x=models.strength, y=models.gamma, c=models.misfit, cmap='viridis_r')
    ax.set_title(f'2-D inversion results, {n} models')
    ax.set_xlabel('Strength')
    ax.set_ylabel(r'Gamma $(\degree$)')
    fig.colorbar(C)
    
    if fname:
        plt.savefig(f'{FIG_DIR}/fname', format='png')
    else:
        plt.show()