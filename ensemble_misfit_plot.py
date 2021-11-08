#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 16:23:29 2021

@author: ja17375
"""
import matplotlib.pyplot as plt
import numpy as np 
import EnsembleVisualiser
FIG_DIR = '/Users/ja17375/Projects/Epac_fast_anom/Figures'

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
   
def plot_2d_sections(Ensemble, ndf, fname=None):
    '''
    Plots 2d sections through the input models (through each axis)
    '''
    m = Ensemble.find_best_fitting(ret=True)
    Ensemble.find_fcrit(ndf)
    
    fig = plt.figure(figsize=(20, 6))
    #make left plot (alpha v strength) - Gamma is fixed
    sec1 = Ensemble.make_section(m[1], axis='gamma')
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.scatter(x=sec1.strength, y=sec1.alpha, c=sec1.misfit, cmap='magma_r',vmin=m[3], vmax=1)
    ax1.set_xlabel('Strength')
    ax1.set_ylabel(r'Alpha ($\degree$)')
    ax1.set_title(f'2D section through gamma = {m[1]:5.3f}')
                
    # make middle plot (gamma v strength )
    sec2 = Ensemble.make_section(m[0], axis='alpha')
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.scatter(x=sec2.strength, y=sec2.gamma, c=sec2.misfit, cmap='magma_r',vmin=m[3], vmax=1)
    ax2.set_xlabel('Strength')
    ax2.set_ylabel(r'Gamma ($\degree$)')
    ax2.set_title(f'2D section through alpha = {m[0]:5.3f}')
    
    #make right plot (alpha v gamma)
    sec3 = Ensemble.make_section(m[2], axis='strength')
    ax3 = fig.add_subplot(1, 3, 3)
    C = ax3.scatter(x=sec3.alpha, y=sec3.gamma, c=sec3.misfit, cmap='magma_r',vmin=m[3], vmax=1)
    ax3.set_xlabel(r'Alpha ($\degree$)')
    ax3.set_ylabel(r'Gamma ($\degree$)')
    ax3.set_title(f'2D section through strength = {m[2]:5.3f}')
    
    #Add model we are slicing through to each section
    ax1.plot(m[2], m[0], marker='x', color='red')
    ax2.plot(m[2], m[1], marker='x', color='red')
    ax3.plot(m[0], m[1], marker='x', color='red')
    fig.suptitle('2-D sections through model space, intersecting at best fitting model (+/- 0.1%)',
                 fontsize=14)
    fig.colorbar(C, ax=ax3)
    ax1.set_xlim([Ensemble.model_config['strength_min'], Ensemble.model_config['strength_max']])
    ax1.set_ylim([Ensemble.model_config['alpha_min'], Ensemble.model_config['alpha_max']])
    ax2.set_xlim([Ensemble.model_config['strength_min'], Ensemble.model_config['strength_max']])
    ax2.set_ylim([Ensemble.model_config['gamma_min'], Ensemble.model_config['gamma_max']])
    ax3.set_xlim([Ensemble.model_config['alpha_min'], Ensemble.model_config['alpha_max']])
    ax3.set_ylim([Ensemble.model_config['gamma_min'], Ensemble.model_config['gamma_max']])
    if fname:
        plt.savefig(f'{FIG_DIR}/{fname}')
    else:
        plt.show() 

if __name__ == '__main__':
    
    E = EnsembleVisualiser.Ensemble('/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/Creasy_tensors/ppv2/SimpleShear_100',
                                           strmax=1)
    m = E.find_best_fitting(ret=True)
    plot_2d_sections(E, m)
