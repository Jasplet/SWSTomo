#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:12:12 2020

@author: ja17375
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import RegularGridInterpolator
from skimage import measure

class Ensemble:
    '''
    Ensemble of models produced by Matisse for a domain
    (this will get more complicated with multiple domains, which may required some
     retooling - i.e figuring out what column is what as the output Ensembles 
     have no column headings)
    '''
    def __init__(self, rundir, read=True, fname='MTS_ensemble.out'):
        '''
        Read raw ensemble and transform it back from the normalise parameter space
        '''
        self.model_config = {'alpha_min': 0, 'alpha_max': 90, 'gamma_min': -180,
                             'gamma_max': 180, 'strength_min': 0, 'strength_max': 0.02}
        self.rundir = rundir
        if read:
            raw_ensemble = self.read_ensemble(rundir,fname)
            
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

    def read_ensemble(self, path, fname):
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

    def evaluate_3d_kde(self, nsamps=50):
        ''' Function to make the full 3D PDF from model Ensemble. As we are using all 3 params
        (alpha, gamma, strength) we do not need to worry about handling different combos'''
        self.alpha_samples = np.linspace(self.model_config['alpha_min'],
                                         self.model_config['alpha_max'],nsamps)
        self.gamma_samples = np.linspace(self.model_config['gamma_min'],
                                         self.model_config['gamma_max'],nsamps)
        self.strength_samples = np.linspace(self.model_config['strength_min'],
                                            self.model_config['strength_max'],nsamps)
        a_samp, g_samp, s_samp = np.meshgrid(
                                    self.alpha_samples,
                                    self.gamma_samples,
                                    self.strength_samples)  
        samples = np.vstack([a_samp.ravel(), g_samp.ravel(), s_samp.ravel()])
        models = np.vstack([self.models.alpha.values,
                            self.models.gamma.values,
                            self.models.strength.values
                            ])
        kernel = stats.gaussian_kde(models)
        self.pdf_3d = np.reshape(kernel(samples), (nsamps, nsamps, nsamps))

    def slice_from_3d_pdf(self, axis, value, n, norm=True):
        '''Interpolates 3D PDF along a fixed axis. 
        I.e generates a 2D slice of the 3D volume '''        
        if axis == 'alpha':
            if (value >= self.model_config['gamma_min']) & (value <= self.model_config['gamma_max']):
                g = np.linspace(self.model_config['gamma_min'],
                                self.model_config['gamma_max'], n)
                s = np.linspace(self.model_config['strength_min'],
                                self.model_config['strength_max'], n)
                gg, ss = np.meshgrid(g, s)
                aa = np.ones(gg.shape) *value
            else:
                raise ValueError(f'Value {value} for alpha is outside limits')  
                
        elif axis == 'gamma':
            if (value >= self.model_config['gamma_min']) & (value <= self.model_config['gamma_max']):
                a = np.linspace(self.model_config['alpha_min'],
                                self.model_config['alpha_max'], n)
                s = np.linspace(self.model_config['strength_min'],
                             self.model_config['strength_max'], n)
                aa, ss = np.meshgrid(a, s)
                gg = np.ones(aa.shape) * value
            else:
                raise ValueError(f'Value {value} for gamma is outside limits')   
                                 
        elif axis == 'strength':
            if (value >= self.model_config['strength_min']) & (value <= self.model_config['strength_max']):
                a = np.linspace(self.model_config['alpha_min'],
                                self.model_config['alpha_max'], n)
                g = np.linspace(self.model_config['gamma_min'],
                                self.model_config['gamma_max'], n)
                aa, gg = np.meshgrid(a, g)
                ss = np.ones(aa.shape) * value
            else:
                raise ValueError(f'Value {value} for strength is outside limits')
        else:    
            raise ValueError(f'{axis} not recognised')
        
        if norm:
            pdf = self.pdf_3d / self.pdf_3d.max()
        else:
            pdf = self.pdf_3d
        interp_function = RegularGridInterpolator((self.alpha_samples,
                                                   self.gamma_samples,
                                                   self.strength_samples),
                                                  pdf)
        points =  np.vstack((aa.flatten(), gg.flatten(), ss.flatten())).T # Stack samples points together. Taking transerve makes it a row 
        slc = interp_function(points).reshape(n, n)   
        return points, slc

    def read_3d_pdf(self, fileID):
        ''' Reads existing 3D PDF '''
        self.pdf_3d = np.load(fileID)
        # check pdf is evenly sampled
        assert self.pdf_3d.shape[0] == self.pdf_3d.shape[1]
        assert self.pdf_3d.shape[0] == self.pdf_3d.shape[2]
        self.nsamps = self.pdf_3d.shape[0]
        # Now we know sampling is even, generate sample values
        self.alpha_samples = np.linspace(self.model_config['alpha_min'],
                                         self.model_config['alpha_max'],self.nsamps)
        self.gamma_samples = np.linspace(self.model_config['gamma_min'],
                                         self.model_config['gamma_max'],self.nsamps)
        self.strength_samples = np.linspace(self.model_config['strength_min'],
                                            self.model_config['strength_max'],self.nsamps)

        
    def save_3d_pdf(self, outfile):
        ''' Saves 3-D PDF as Numpy Array'''
        np.save(f'{self.rundir}/{outfile}', self.pdf_3d)
    
    
    
    
