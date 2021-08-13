#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 17:21:47 2021

@author: ja17375
"""
import numpy as np
from multiprocessing import Pool
import EnsembleVisualiser
import matplotlib.pyplot as plt

def resample_2d_gs_slice(models, a_int, config):
    '''Resample models to make a 2-D slice that intersect at a givien input model (m_int) '''
    y_samps = np.arange(config['gamma_min'], config['gamma_max']+ config['gamma_step'],
                        config['gamma_step'])
    x_samps = np.arange(config['strength_min'], config['strength_max'] + config['strength_step'],
                        config['strength_step'])
    # Make slice through alpha axis
    slc = np.zeros((x_samps.size, y_samps.size))
    for i, s in enumerate(x_samps):
        # Flip y samples to fix indexing starting at top left corner
        for j, g in enumerate(y_samps):
            node = [a_int, g, s]
            misfit = ellip_dist_avg(models, node, config)
            slc[i, j] = misfit 
    return x_samps, y_samps, slc.T

def resample_2d_as_slice(models, g_int, config):
    '''Resample models to make 2-D splice along the gamma axis fixed at a given g_int'''
    y_samps = np.arange(config['alpha_min'], config['alpha_max'] + config['alpha_step'],
                        config['alpha_step'])
    x_samps = np.arange(config['strength_min'], config['strength_max'] + config['strength_step'],
                        config['strength_step'])
    # Make slice through alpha axis
    slc = np.zeros((x_samps.size, y_samps.size))
    for i, s in enumerate(x_samps):
        # Flip y samples to fix indexing starting at top left corner
        for j, a in enumerate(y_samps):
            node = [a, g_int, s]
            misfit = ellip_dist_avg(models, node, config)
            slc[i, j] = misfit 
    return x_samps, y_samps, slc.T

def resample_2d_ag_slice(models, s_int, config):
    '''Resample models to make a 2D slice along the strength axis fixed at a given s_int 
    (so slice is gamma v alpha ''' 
    y_samps = np.arange(config['gamma_min'], config['gamma_max'] + config['gamma_step'],
                        config['gamma_step'])
    x_samps = np.arange(config['alpha_min'], config['alpha_max'] + config['alpha_step'],
                        config['alpha_step'])
    # Make slice through alpha axis
    print(f'Slicing through s = {s_int}')
    slc = np.zeros((x_samps.size, y_samps.size))
    for i, a in enumerate(x_samps):
        # Flip y samples to fix indexing starting at top left corner
        for j, g in enumerate(y_samps):
#                print(f'{i}: {a}, {j}: {g}')
            node = [a, g, s_int]
            misfit = ellip_dist_avg(models, node, config)
            slc[i, j] = misfit 
    return x_samps, y_samps, slc.T
    
def ellip_dist_avg(models, node, config):
    '''Finds all models within the sample ellipse of the requested node and returns the weighted average misfit '''
    a_ellip = config['alpha_step']*1.25
    g_ellip = config['gamma_step']*1.25
    s_ellip = config['strength_step']*1.25
    xa = (models.alpha - node[0])/a_ellip   
    ya = (ang_dist(models.gamma, node[1]))/g_ellip
    za = (models.strength - node[2])/s_ellip
    r = xa**2 + ya**2 + za**2 
    # Add radial distance to models DataFrame
    models['r'] = r
    mods_in_elli = models[models.r <= 1]
    if mods_in_elli.empty is False:
        mf_avg = np.average(mods_in_elli.misfit, weights=mods_in_elli.r)
    else:
        mf_avg = heal(models, node, config)
                
    return mf_avg
    
def ang_dist(mod_ang, node_ang):
    mod_ang360 = mod_ang + 180
    d = ((mod_ang360 - node_ang) % 360) -180 
    return abs(d)

def heal(models, node, hconfig):
    '''Heal points in slice which are NaNs, this is done by recursivly expanding the ellipse until
       we find some model points to average. Ensuring smoothness across model 
    ''' 
    hconfig['alpha_step'] = hconfig['alpha_step']*1.1
    hconfig['gamma_step'] = hconfig['gamma_step']*1.1
    hconfig['strength_step'] = hconfig['strength_step']*1.1
    print(f'Healing node {node}')  
    print('New Steps {} {} {}'.format(hconfig['alpha_step'], hconfig['gamma_step'],
                                      hconfig['strength_step']))
    mf = ellip_dist_avg(models, node, hconfig)
    return mf 
if __name__ == '__main__':
    
    rundir = '/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/ppv_model'
    Ensemble = EnsembleVisualiser.Ensemble(rundir, strmax=0.25)
    config = Ensemble.model_config
    config['gamma_min'] = -180
    config['gamma_max'] = 180 
    config['alpha_min'] = 0
    config['alpha_max'] = 90
    config['alpha_step'] = 10
    config['gamma_step'] = 21
    config['strength_step'] = 0.005
    a, g, slc = resample_2d_ag_slice(Ensemble.models, 0.1036, config)
    aa, gg = np.meshgrid(a, g)
    C = plt.contourf(aa, gg, slc, levels=20, cmap='viridis_r')
    plt.colorbar(C)
    plt.show()