#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:12:12 2020

@author: ja17375
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from mayavi import mlab

class Ensemble:
    '''
    Ensemble of models produced by Matisse for a domain
    (this will get more complicated with multiple domains, which may required some
     retooling - i.e figuring out what column is what as the output Ensembles 
     have no column headings)
    '''
    def __init__(self, rundir):
        '''
        Read raw ensemble and transform it back from the normalise parameter space
        '''
        self.model_config = {'alpha_min': 0, 'alpha_max': 90, 'gamma_min': -180,
                             'gamma_max': 180, 'strength_min': 0, 'strength_max': 0.02}
        
        raw_ensemble = self.read_ensemble(rundir,'MTS_dev_ensemble.out')
        
        alpha = self.restore_params(raw_ensemble[:,0],
                            self.model_config['alpha_min'], self.model_config['alpha_max'])
        gamma = self.restore_params(raw_ensemble[:,1],
                            self.model_config['gamma_min'], self.model_config['gamma_max'])
        strength = self.restore_params(raw_ensemble[:,2],
                            self.model_config['strength_min'], self.model_config['strength_max'])
        
        self.models = pd.DataFrame({'alpha': alpha, 'gamma': gamma,
                                   'strength':strength, 'misfit': raw_ensemble[:,-1]})
        # Store modes as a DataFrame for ease of reference. Misfit index as [:,-1] in order to
        # "future proof" as misfit is always last column of MTS_ensemble.out 

    def read_ensemble(self, path, fname='MTS_ensemble.out'):
        '''function to read the MTS_emsemble from the input run directory'''
        raw_ensemble = np.loadtxt(f'{path}/{fname}')
        
        return raw_ensemble

    def restore_params(self, param, p_min, p_max):
        ''' 
        Restores true value of parameters in each model from their normalised state 
        '''
        p_range = p_max - p_min 
        true_params = p_min + (p_range * param)
        
        return true_params
    
    def evaluate_2d_kde(Ensemble, params, nsamps=50j):
    
        x_param = params[0]
        y_param = params[1]
        xmin = Ensemble.model_config[f'{x_param}_min']
        xmax = Ensemble.model_config[f'{x_param}_max']
        ymin = Ensemble.model_config[f'{y_param}_min']
        ymax = Ensemble.model_config[f'{y_param}_max']
        
        x_samp, y_samp = np.mgrid[xmin:xmax:nsamps, ymin:ymax:nsamps]
        samps = np.vstack([x_samp.ravel(), y_samp.ravel()])
        models = np.vstack([Ensemble.models[x_param].values, Ensemble.models[y_param].values])
        # stack model points for alpha, gamma together into a 2d array for evalution in KDE
        kernel = stats.gaussian_kde(models)
        # do the gaussian KDE for our model values, now we need to evaluate PDF at sample points
        pdf = np.reshape(kernel(samps).T, x_samp.shape)
    return pdf 

    def evaluate_3d_kde(Ensemble, nsamps=50j):
        ''' Function to make the full 3D PDF from model Ensemble. As we are using all 3 params
        (alpha, gamma, strength) we do not need to worry about handling different combos'''
        a_samp, g_samp, s_samp = np.mgrid[
            Ensemble.model_config['alpha_min'] : Ensemble.model_config['alpha_max'] : nsamps,
            Ensemble.model_config['gamma_min'] : Ensemble.model_config['gamma_max'] : nsamps,
            Ensemble.model_config['strength_min'] : Ensemble.model_config['strength_max'] : nsamps    
            ]
        samples = np.vstack([a_samp.ravel(), g_samp.ravel(), s_samp.ravel()])
        models = np.vstack([Ensemble.models.alpha.values,
                            Ensemble.models.gamma.values,
                            Ensemble.models.strength.values
                            ])
        kernel = stats.gaussian_kde(models)
        self.pdf = np.reshape(kernel(samples), a_samp.shape)
        self.samples = samples
        
def plot_ags_models(Ensemble):
    '''
    Simple 3-D plot of all models in Ensemble, assuming the 3 pametees alpha, gamma, strength
    have been inverted for 
    '''    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    C = ax.scatter(Ensemble.models.alpha, Ensemble.models.gamma, Ensemble.models.strength, marker='.',
               c=Ensemble.models.misfit, cmap='viridis')
    ax.set_xlim([Ensemble.model_config['alpha_min'], Ensemble.model_config['alpha_max']])
    ax.set_ylim([Ensemble.model_config['gamma_min'], Ensemble.model_config['gamma_max']])
    ax.set_zlim([Ensemble.model_config['strength_min'], Ensemble.model_config['strength_max']])
    plt.colorbar(C)


    
def plot_alpha_gamma_pdf(Ensemble, pdf):
    '''Plots 2D pdf for the parameters  alpha, gamma from model ensemble'''
    # sample PDF at sample points, reshape to fit grid, take transeverse or array so its the right way round
    fig, ax = plt.subplots()
    C =  ax.contourf(pdf, cmap=plt.cm.gist_earth_r,
              extent=[Ensemble.model_config['alpha_min'], Ensemble.model_config['alpha_max'],
                      Ensemble.model_config['gamma_min'], Ensemble.model_config['gamma_max']],
                      )
    plt.colorbar(C)
#    ax.plot(Ensemble.models.alpha, Ensemble.models.gamma, '.')
    
def plot_3d_pdf(Ensemble):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    alpha = Ensemble.samples[0,:]
    gamam = Ensemble.samples[1,:]
    strength = Ensemble.sample[2,:]
    ax.scatter(a_samp, g_samp, s_samp, c=pdf)

def plot_3d_isosurface(Ensemble, pdf, nsamps):
    
    a_samp, g_samp, s_samp = np.mgrid[
    Ensemble.model_config['alpha_min'] : Ensemble.model_config['alpha_max'] : nsamps,
    Ensemble.model_config['gamma_min'] : Ensemble.model_config['gamma_max'] : nsamps,
    Ensemble.model_config['strength_min'] : Ensemble.model_config['strength_max'] : nsamps    
    ]
    mlab.contour3d(a_samp, g_samp, s_samp, pdf, opacity=0.5)
    mlab.axes()
    mlab.show()
    
def plot_pdf_slice():
    '''plot a 2d slice through a 3d PDF'''
    pass

def plot_jelly_bowl():
    '''plot a spherical projection of PDF using strenght as radius and alpha, gamma as theta, phi'''
    pass
    
    
    
    
    
