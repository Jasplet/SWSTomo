#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 16:17:28 2021

@author: ja17375
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats 
import resampler
from EnsembleVisualiser import Ensemble

FIG_DIR = '/Users/ja17375/Projects/Epac_fast_anom/Figures/ModelSlices/'

def ftest(lam2min,ndf,k=2,alpha=0.05):
    """
    returns lambda2 value at 100(1-alpha)% confidence interval
    by default alpha = 0.05 = 95% confidence interval
    following Silver and Chan (1991) [modifications by walsh et al., 2013]
    As we are dealing with traces that have alreayd been passed through SHEBA,
    we do not need to check (or calculate) degrees of freedom as this has alreay
    been done.

    Needed for pair_stack to calc lam2alpha for SKS and SKKS
    """
    F = stats.f.ppf(1-alpha,k,ndf)
    lam2alpha = lam2min * ( 1 + (k/(ndf-k)) * F)
    return lam2alpha 

def plot_sections(Ensemble, model, config, err):
      
    clevels = [0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    fig = plt.figure(figsize=(22,6)) 
    fcrit = Ensemble.fcrit
    print(fcrit)
    # 2D slice through alpha axis
    ax1 = fig.add_subplot(131)
    g_int = model[1]
    s, a, slc = resampler.resample_2d_as_slice(Ensemble.models, g_int, config)
    ss, aa = np.meshgrid(s, a)
    CS1 = ax1.contour(ss, aa, slc, levels=clevels, colors='k')
    ax1.contour(ss, aa, slc, levels=[fcrit], linewidths=3., colors='k')
    ax1.clabel(CS1, inline=True)
    ax1.set_xlim([config['strength_min'], config['strength_max']])
    ax1.set_ylim([0,90])
    ax1.set_xticks(np.linspace(config['strength_min'], config['strength_max'], 6))
    ax1.set_yticks(np.linspace(0,90,7))
    ax1.set_xlabel('Strength', fontsize=14)
    ax1.set_ylabel(r'$\alpha$ (°)', fontsize=14)
    ax1.set_title(r'Slice through $\gamma$ = {:4.3f}°'.format(model[1]), fontsize=16)
    ax1.tick_params(labelsize=16)
    
    ax2 = fig.add_subplot(132)
    a_int = model[0]
    s, g, slc = resampler.resample_2d_gs_slice(Ensemble.models, a_int, config)
    ss, gg = np.meshgrid(s, g)
    CS2 = ax2.contour(ss,gg, slc, levels=clevels, colors='k')
    ax2.contour(ss, gg, slc, levels=[fcrit], linewidths=3., colors='k')
    ax2.clabel(CS2, inline=True)
    ax2.set_xlim([config['strength_min'], config['strength_max']])
    ax2.set_ylim([-180, 180])
    ax2.set_xticks(np.linspace(config['strength_min'], config['strength_max'], 6))
    ax2.set_yticks(np.linspace(config['gamma_min'], config['gamma_max'], 7))
    ax2.set_xlabel('Strength', fontsize=14)
    ax2.set_ylabel(r'$\gamma$ (°)',fontsize=14)
    ax2.set_title(r'Slice through $\alpha$ = {:4.0f}°'.format(model[0]))
    ax2.tick_params(labelsize=16)
    
    ax3 = fig.add_subplot(133)
    s_int = model[2]
    a, g, slc = resampler.resample_2d_ag_slice(Ensemble.models, s_int, config)
    aa, gg = np.meshgrid(a, g)
    CS3 = ax3.contour(aa, gg, slc, levels=clevels, colors='k')
    ax3.contour(aa, gg, slc, levels=[fcrit], linewidths=3., colors='k')
    ax3.clabel(CS3, inline=True)
    ax3.set_ylim([-180,180])
    ax3.set_yticks(np.linspace(config['gamma_min'], config['gamma_max'], 7))
    ax3.set_xlim([0,90])
    ax3.set_xticks(np.linspace(0,90,7))
    ax3.set_xlabel(r'$\alpha$ (°)',fontsize=14)
    ax3.set_ylabel(r'$\gamma$ (°)',fontsize=14)
    ax3.set_title(f'Slice through strength = {model[2]:4.3f}')
    ax3.tick_params(labelsize=16)
    
    #Add best model
    ax1.errorbar(s_int, a_int, fmt='bx', yerr=err[0], xerr=err[2], markersize=15)
    ax2.errorbar(s_int, g_int, fmt='bx', yerr=err[1], xerr=err[2], markersize=15)
    ax3.errorbar(a_int, g_int, fmt='bx', yerr=err[1], xerr=err[0], markersize=15)
    # fig.colorbar(C, ax=[ax1,ax2,ax3])
    fig.tight_layout()
          
    # SMALL_SIZE = 12
    # MEDIUM_SIZE = 16
    
    # plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    # plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    # plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    # plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    # plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    # plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    
    return fig

if __name__ == '__main__':

     # Add sample steps to Ensemble config
    samp_config = {'alpha_min': 0, 'alpha_max':90, 'alpha_step':1,
                   'gamma_min':-180, 'gamma_max':180, 'gamma_step': 2,
                   'strength_min':0, 'strength_max': 0.02, 'strength_step': 0.0005}  
    path = '/Users/ja17375/Projects/Epac_fast_anom/HQ_data/ScS_fix_test'
    #Elliptical model
    E = Ensemble(f'{path}/ellipTI', strmax=samp_config['strength_max'])
    E.find_fcrit(ndf=239)
    mod = E.find_best_fitting(ret=True)
    TI_err = [6, 2, 0.002]
    fig1 = plot_sections(E, mod, samp_config, TI_err)
    fig1.savefig(f'{FIG_DIR}/Ellip_misfit_2D_slice.eps')
    # Bridgmanite Model
    samp_config['strength_max'] = 0.5
    samp_config['strength_step'] = 0.005
    E = Ensemble(f'{path}/pv_100_001', strmax=samp_config['strength_max'])
    E.find_fcrit(ndf=239)
    mod = E.find_best_fitting(ret=True)
    bm_err = [5, 5, 0.06]
    fig2 = plot_sections(E, mod, samp_config, bm_err)
    fig2.savefig(f'{FIG_DIR}/pv_100_001_misfit_2D_slice.eps')
    # ppv 100_001 Model
    samp_config['strength_max'] = 0.5
    samp_config['strength_step'] = 0.005
    E = Ensemble(f'{path}/ppv_001_100', strmax=samp_config['strength_max'])
    E.find_fcrit(ndf=239)
    mod = E.find_best_fitting(ret=True)
    ppv1_err = [2, 4, 0.06]
    fig3 = plot_sections(E, mod, samp_config, ppv1_err)
    fig3.savefig(f'{FIG_DIR}/ppv_001_100_misfit_2D_slice.eps')
    # ppv 100_010 Model
    samp_config['strength_max'] = 0.5
    E = Ensemble(f'{path}/ppv_010_100', strmax=samp_config['strength_max'])
    E.find_fcrit(ndf=239)
    ppv2_err = [2, 2, 0.02]
    mod = E.find_best_fitting(ret=True)
    fig4 = plot_sections(E, mod, samp_config, ppv2_err)
    fig4.savefig(f'{FIG_DIR}/ppv_010_100_misfit_2D_slice.eps')