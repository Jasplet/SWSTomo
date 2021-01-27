#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:12:12 2020

@author: ja17375
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator

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
            # raw ensemble also has our function f(x) that approximates (or is proportional to) the true PDF
            # evaluated for each model. f(x) = 1 / sum(lam2) for all paths in the inversion
            # to get the model misfit (i.e the sum of the lam2 residuals) we need to take the reciprocal
            misfit = 1 / raw_ensemble[:, -1]
            self.models = pd.DataFrame({'alpha': alpha, 'gamma': gamma,
                                       'strength':strength, 'misfit': misfit})
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

    def find_most_likely(self, ret=False):
        '''
        finds the most likely (i.e. maximum) value in model ensemble. 
        this method requires that the object has self.pdf_3d (either read in or generated) 
        '''

        idxmax = self.pdf_3d.argmax()
        # Translate single index to the indicies of the 3D array
        ags = np.unravel_index(idxmax, self.pdf_3d.shape)
        # Make an array holding the most likely model
        self.most_likely_model = [self.alpha_samples[ags[0]],
                                  self.gamma_samples[ags[1]],
                                  self.strength_samples[ags[2]]
                                  ]
        if ret:    
            return self.most_likely_model
        elif not ret:
            print('The most likely model for the given Model Ensemble is:')
            print(f'Alpha = {self.most_likely_model[0]:.3f}')
            print(f'Gamma = {self.most_likely_model[1]:.3f}')
            print(f'Strength = {self.most_likely_model[2]:.3f}')
 
    def find_best_fitting(self, ret=False):
        '''
        Finds the best fitting model from the model ensemble 
        '''
        imin = self.models.misfit.idxmin()
        self.best_fitting_model = self.models.iloc[imin]
        if ret:
            return self.best_fitting_model
        elif not ret:
            print('The best fitting model is:')
            print(f'Alpha = {self.best_fitting_model[0]:.3f}')
            print(f'Gamma = {self.best_fitting_model[1]:.3f}')
            print(f'Strength = {self.best_fitting_model[2]:.3f}')
            print(f'Model misfit = {self.best_fitting_model[3]:.3f}')
    
    def slice_3d_volume(self, axis, value, n, method):
        '''
        Interpolates 3D models along a give axis to return a regularaly sampled 2-D slice
        This can be done to interpolate either model misfit or the 3D PDF (if the PDF has been generated)
        
        '''        
        if method == 'pdf':
            pdf = self.pdf_3d / self.pdf_3d.max()
            interp_function = RegularGridInterpolator((self.alpha_samples,
                                                   self.gamma_samples,
                                                   self.strength_samples),
                                                      pdf)
        elif method == 'misfit':
            misfit = self.models.misfit
            models = np.vstack([self.models.alpha,
                               self.models.gamma,
                               self.models.strength]).T
            interp_function = LinearNDInterpolator(models, misfit)
        
        aa, gg, ss = self.gen_2d_param_grid(axis, value, n)
        points =  np.vstack((aa.flatten(), gg.flatten(), ss.flatten())).T # Stack samples points together. Taking transerve makes it a row 
        slc = interp_function(points).reshape(n, n)   
        return points, slc

    def gen_2d_param_grid(self, axis, value, n):
        '''
        Function that makes regular grids of models parameters for where one parameter (or axis) of the model 
        is fixed.
        This is to facilitate the interpolation of 2-D slices through a 3-D volume
        '''
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

        return aa, gg, ss
    
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
    
    
