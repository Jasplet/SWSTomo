#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 16:17:28 2021

@author: ja17375
"""
import matplotlib.pyplot as plt
import numpy as np
from l2stats import ftest
import resampler
from EnsembleVisualiser import Ensemble

FIG_DIR = '/Users/ja17375/SWSTomo/Figures/ModelSlices/Contoured'
 

def plot_sections(Ensemble, model, config, err):
    
    clevels = [0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1.0,2.0]
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
    ax1.clabel(CS1, inline=True, fontsize=12)
    ax1.set_xlim([config['strength_min'], config['strength_max']])
    ax1.set_ylim([0,90])
    ax1.set_xticks(np.linspace(config['strength_min'], config['strength_max'], 6))
    ax1.set_yticks(np.linspace(0,90,7))
    ax1.set_xlabel('Strength',fontsize=16)
    ax1.set_ylabel(r'$\alpha$ (°)', fontsize=16)
    ax1.set_title(r'Slice through $\gamma$ = {:4.3f}°'.format(model[1]), fontsize=16)
    ax1.tick_params(labelsize=16)
    
    ax2 = fig.add_subplot(132)
    a_int = model[0]
    s, g, slc = resampler.resample_2d_gs_slice(Ensemble.models, a_int, config)
    ss, gg = np.meshgrid(s, g)
    CS2 = ax2.contour(ss,gg, slc, levels=clevels, colors='k')
    ax2.contour(ss, gg, slc, levels=[fcrit], linewidths=3., colors='k')
    ax2.clabel(CS2, inline=True, fontsize=12)
    ax2.set_xlim([config['strength_min'], config['strength_max']])
    ax2.set_ylim([-180, 180])
    ax2.set_xticks(np.linspace(config['strength_min'], config['strength_max'], 6))
    ax2.set_yticks(np.linspace(config['gamma_min'], config['gamma_max'], 7))
    ax2.set_xlabel('Strength', fontsize=16)
    ax2.set_ylabel(r'$\gamma$ (°)', fontsize=16)
    ax2.set_title(r'Slice through $\alpha$ = {:4.3f}°'.format(model[0]), fontsize=16) 
    ax2.tick_params(labelsize=16)
    
    ax3 = fig.add_subplot(133)
    s_int = model[2]
    a, g, slc = resampler.resample_2d_ag_slice(Ensemble.models, s_int, config)
    aa, gg = np.meshgrid(a, g)
    CS3 = ax3.contour(aa, gg, slc, levels=clevels, colors='k')
    ax3.contour(aa, gg, slc, levels=[fcrit], linewidths=3., colors='k')
    ax3.clabel(CS3, inline=True, fontsize=12)
    ax3.set_ylim([-180,180])
    ax3.set_yticks(np.linspace(config['gamma_min'], config['gamma_max'], 7))
    ax3.set_xlim([0,90])
    ax3.set_xticks(np.linspace(0,90,7))
    ax3.set_xlabel(r'$\alpha$ (°)', fontsize=16)
    ax3.set_ylabel(r'$\gamma$ (°)', fontsize=16)
    ax3.set_title(f'Slice through strength = {model[2]:4.3f}', fontsize=16)
    ax3.tick_params(labelsize=16)
    
    #Add best model
    ax1.errorbar(s_int, a_int, 'kx', yerr=err[0], xerr=err[2], markersize=10)
    ax2.plot(s_int, g_int, 'kx', yerr=err[1], xerr=err[2], markersize=10)
    ax3.plot(a_int, g_int, 'kx', yerr=err[1], xerr=err[0], markersize=10)
    # fig.colorbar(C, ax=[ax1,ax2,ax3])
    fig.tight_layout()
    
    return fig

if __name__ == '__main__':
     # Add sample steps to Ensemble config
    samp_config = {'alpha_min': 0, 'alpha_max':90, 'alpha_step':1,
                   'gamma_min':-180, 'gamma_max':180, 'gamma_step': 2,
                   'strength_min':0, 'strength_max': 0.02, 'strength_step': 0.0005}   
    #Elliptical model
    rundir = '/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/EllipTI/'
    E = Ensemble(rundir, strmax=samp_config['strength_max'])
    E.find_fcrit(ndf=239)
    mod = E.find_best_fitting(ret=True)
    TI_err = [4.777, 8.630, 0.002]
    fig1 = plot_sections(E, mod, samp_config, TI_err)
    fig1.savefig(f'{FIG_DIR}/Ellip_misfit_2D_slice.png')
    # Bridgmanite Model
    samp_config['strength_max'] = 0.5
    samp_config['strength_step'] = 0.005
    rundir = '/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/pv_100_001/'
    E = Ensemble(rundir, strmax=samp_config['strength_max'])
    E.find_fcrit(ndf=239)
    mod = E.find_best_fitting(ret=True)
    bm_err = [4.517, 31.934, 0.053]
    fig2 = plot_sections(E, mod, samp_config, bm_err)
    fig2.savefig(f'{FIG_DIR}/pv_001_100_misfit_2D_slice.png')
    # ppv 100_001 Model
    samp_config['strength_max'] = 0.5
    samp_config['strength_step'] = 0.005
    rundir = '/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/ppv_001_100/'
    E = Ensemble(rundir, strmax=samp_config['strength_max'])
    E.find_fcrit(ndf=239)
    mod = E.find_best_fitting(ret=True)
    ppv1_err = [9.463, 35.148, 0.071]
    fig3 = plot_sections(E, mod, samp_config)
    fig3.savefig(f'{FIG_DIR}/ppv_100_001_misfit_2D_slice.png')
    # ppv 100_010 Model
    samp_config['strength_max'] = 0.5
    rundir = '/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/ppv_010_100/'
    E = Ensemble(rundir, strmax=samp_config['strength_max'])
    E.find_fcrit(ndf=239)
    ppv2_err = [4.352, 60.047, 0.051]
    mod = E.find_best_fitting(ret=True)
    fig4 = plot_sections(E, mod, samp_config)
    fig4.savefig(f'{FIG_DIR}/ppv_100_010_misfit_2D_slice.png')