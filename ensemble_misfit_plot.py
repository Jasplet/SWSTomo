#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 16:23:29 2021

@author: ja17375
"""
import matplotlib.pyplot as plt
import numpy as np 
import EnsembleVisualiser
FIG_DIR = '/Users/ja17375/SWSTomo/Figures'

def plot_nbest_models(Ensemble, n, full_range=False):
    
    models_sort = Ensemble.models.sort_values(by=['misfit'])
    pmodels = models_sort.iloc[0:n]
    # sort again so lowest msfits are at the end and therefore drawn on top!
    models = pmodels.sort_values(by=['misfit'], ascending=False)
    fig1 = plot_3d_models(models)
    (ax_3d) = fig1.axes[0]
    fig2 = plot_2d_views(models)
    (ax1, ax2, ax3, cbar) = fig2.get_axes()
    
    if full_range:
        ax_3d.set_xlim([Ensemble.model_config['alpha_min'],
                        Ensemble.model_config['alpha_max']])
        ax_3d.set_ylim([Ensemble.model_config['gamma_min'],
                        Ensemble.model_config['gamma_max']])
        ax_3d.set_zlim([Ensemble.model_config['strength_min'],
                        Ensemble.model_config['strength_max']])
        
        ax1.set_xlim([Ensemble.model_config['strength_min'],
                      Ensemble.model_config['strength_max']])
        ax1.set_ylim([Ensemble.model_config['alpha_min'],
                  Ensemble.model_config['alpha_max']])
        ax2.set_xlim([Ensemble.model_config['strength_min'],
                  Ensemble.model_config['strength_max']])
        ax2.set_ylim([Ensemble.model_config['gamma_min'],
                  Ensemble.model_config['gamma_max']]) 
        ax3.set_xlim([Ensemble.model_config['alpha_min'],
          Ensemble.model_config['alpha_max']]) 
        ax3.set_ylim([Ensemble.model_config['gamma_min'],
          Ensemble.model_config['gamma_max']]) 
    # Add best fitting model to both figs
    best_model = Ensemble.find_best_fitting(ret=True)
    ax_3d.plot(xs=best_model.alpha, ys=best_model.gamma,
                 zs=best_model.strength, marker='x', color='red')
    ax1.plot(best_model.strength, best_model.alpha, marker='x', color='red')
    ax2.plot(best_model.strength, best_model.gamma, marker='x', color='red')
    ax3.plot(best_model.alpha, best_model.gamma, marker='x', color='red')
    plt.show()
    
def plot_3d_models(models):
    '''
    Makes a 3-D (x,y,z) plot showing the n best fitting models
    Args:
        Ensemble (Class) - The ensemble object holding all the models 
        n (float) - number of models to plot
    Returns:
        A figure object containgng the 3-D plot
    '''
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    C = ax.scatter(models.alpha, models.gamma, models.strength, c=models.misfit,s=6, cmap='viridis_r')
    ax.set_xlabel(r'$\alpha \degree $')
    ax.set_ylabel(r'$\gamma \degree $')    
    ax.set_zlabel('strength')
    cbar = fig.colorbar(C, ax=ax)
    cbar.ax.set_ylabel('Model misfit')
    return fig

def plot_2d_views(models):
    ''' 
    Makes a 3-panel figure showing 2-D views of the input models (i.e along each axis)
    '''
    fig = plt.figure(figsize=(20, 6))
    # plot alpha v strength
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.scatter(x=models.strength, y=models.alpha, c=models.misfit, s=6, cmap='viridis_r')
    ax1.set_xlabel('Strength')
    ax1.set_ylabel(r'Alpha ($\degree$)')
    # plot gammav strength
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.scatter(x=models.strength, y=models.gamma, c=models.misfit, s=6, cmap='viridis_r')
    ax2.set_xlabel('Strength')
    ax2.set_ylabel(r'Gamma ($\degree$)')
    # plot gamma v alpha
    ax3 = fig.add_subplot(1, 3, 3)
    C = ax3.scatter(x=models.alpha, y=models.gamma, c=models.misfit, s=6, cmap='viridis_r')
    ax3.set_xlabel('Alpha ($\degree$)')
    ax3.set_ylabel('Gamma ($\degree$)')
    cbar = fig.colorbar(C)
    cbar.ax.set_ylabel(r'$\Sigma (\lambda_2)$')
    return fig
    

if __name__ == '__main__':
    
    Ensemble = EnsembleVisualiser.Ensemble('/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/')
    m = Ensemble.find_best_fitting(ret=True)
    plot_2d_sections(Ensemble.models, m)
