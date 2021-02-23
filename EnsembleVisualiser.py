#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:12:12 2020

@author: ja17375
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from l2stats import ftest
FIG_DIR = '/Users/ja17375/SWSTomo/Figures'


class Ensemble:
    '''
    Ensemble of models produced by Matisse for a domain
    (this will get more complicated with multiple domains, which may required some
     retooling - i.e figuring out what column is what as the output Ensembles 
     have no column headings)
    '''
    def __init__(self, rundir, read=True, dims='ags', strmax=0.02, fname='MTS_ensemble.out'):
        '''
        Read raw ensemble and transform it back from the normalise parameter space
        '''
        self.model_config = {'alpha_min': 0, 'alpha_max': 90, 'gamma_min': -180,
                             'gamma_max': 180, 'strength_min': 0, 'strength_max': strmax}

        self.rundir = rundir
        if read:
            raw_ensemble = self.read_ensemble(rundir,fname)
            if dims == 'ags':
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
                print('Model Config set as:')
                print('{} < alpha < {}'.format(self.model_config['alpha_min'], self.model_config['alpha_max']))
                print('Beta = fixed')
                print('{} < gamma < {}'.format(self.model_config['gamma_min'], self.model_config['gamma_max']))
                print('{} < strength < {}'.format(self.model_config['strength_min'],
                                                  self.model_config['strength_max']))
                      
            elif dims == 'ag':
                #2-D alpha and gamma only
                alpha = self.restore_params(raw_ensemble[:,0],
                    self.model_config['alpha_min'], self.model_config['alpha_max'])
                gamma = self.restore_params(raw_ensemble[:,1],
                                    self.model_config['gamma_min'], self.model_config['gamma_max'])
                misfit = 1 / raw_ensemble[:, -1]
                self.models = pd.DataFrame({'alpha': alpha, 'gamma': gamma,
                                            'misfit': misfit})
            elif dims == 'as':
                alpha = self.restore_params(raw_ensemble[:,0],
                    self.model_config['alpha_min'], self.model_config['alpha_max'])
                strength = self.restore_params(raw_ensemble[:,1],
                                    self.model_config['strength_min'], self.model_config['strength_max'])
                misfit = 1 / raw_ensemble[:, -1]
                self.models = pd.DataFrame({'alpha': alpha, 'strength': strength,
                                           'misfit': misfit})
            elif dims == 'gs':
                gamma = self.restore_params(raw_ensemble[:,0],
                                    self.model_config['gamma_min'], self.model_config['gamma_max'])
                strength = self.restore_params(raw_ensemble[:,1],
                                    self.model_config['strength_min'], self.model_config['strength_max'])
                misfit = 1 / raw_ensemble[:, -1]
                self.models = pd.DataFrame({'gamma': gamma,
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

    def find_fcrit(self, ndf):
        '''Find the f-test 95% confidence value  
        
        Args:
            ndf (int) - degrees of freedom of input data (ndf is calculated by sheba for each path in advance)
        '''
#       number of model dims is 1 less than cols of ensemble 
        k = self.models.shape[1] - 1 
        best = self.find_best_fitting(ret=True)
        # last col of model is always misfit
        self.fcrit = ftest(best[-1], ndf,k)
        
        
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
            print(f'Alpha = {self.best_fitting_model[0]:.5f}')
            print(f'Gamma = {self.best_fitting_model[1]:.5f}')
            print(f'Strength = {self.best_fitting_model[2]:.5f}')
            print(f'Model misfit = {self.best_fitting_model[3]:.5f}')
       
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
    
    def plot_2d_sections(self, ndf, fname=None):
        '''
        Plots 2d sections through the input models (through each axis)
        '''
        m = self.find_best_fitting(ret=True)
        self.find_fcrit(ndf)
        
        fig = plt.figure(figsize=(20, 6))
        #make left plot (alpha v strength) - Gamma is fixed
        sec1 = self.make_section(m[1], axis='gamma')
        (xi, yi, asgrid) = self.grid_section(sec1, axis='gamma')
        ax1 = fig.add_subplot(1, 3, 1)
        # ax1.contour(xi, yi, asgrid, levels=[self.fcrit])
        # ax1.contourf(xi, yi, asgrid, levels=10, cmap='viridis_r')
        ax1.scatter(x=sec1.strength, y=sec1.alpha, c=sec1.misfit, cmap='viridis_r',vmin=m[3])
        ax1.set_xlabel('Strength')
        ax1.set_ylabel(r'Alpha ($\degree$)')
        ax1.set_title(f'2D section through gamma = {m[1]:5.3f}')
                    
        # make middle plot (gamma v strength )
        sec2 = self.make_section(m[0], axis='alpha')
        # (xi, yi, gsgrid) = self.grid_section(sec2, axis='alpha')
        ax2 = fig.add_subplot(1, 3, 2)
        # ax2.contourf(xi, yi, gsgrid, levels=10, cmap='viridis_r')
        # ax2.contour(xi, yi, gsgrid, [self.fcrit])
        ax2.scatter(x=sec2.strength, y=sec2.gamma, c=sec2.misfit, cmap='viridis_r',vmin=m[3])
        ax2.set_xlabel('Strength')
        ax2.set_ylabel(r'Gamma ($\degree$)')
        ax2.set_title(f'2D section through alpha = {m[0]:5.3f}')
        
        #make right plot (alpha v gamma)
        sec3 = self.make_section(m[2], axis='strength')
        # (xi, yi, aggrid) = self.grid_section(sec3, axis='strength')
        ax3 = fig.add_subplot(1, 3, 3)
        # ax3.contourf(xi, yi, aggrid, levels=10, cmap='viridis_r')
        # ax3.contour(xi, yi, aggrid, [self.fcrit])
        C = ax3.scatter(x=sec3.alpha, y=sec3.gamma, c=sec3.misfit, cmap='viridis_r',vmin=m[3])
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
        if fname:
            plt.savefig(f'{FIG_DIR}/{fname}')
        else:
            plt.show()
    
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
    
    def grid_misfit_volume(self):
        '''Grids the models so we can make pretty plots (i hope) '''
        a = self.models.alpha.values
        g = self.models.gamma.values
        s = self.models.strength.values
        mfit = self.models.misift.values
        xi, yi, zi = np.mgrid[self.model_config['alpha_min']:self.model_config['alpha_max']:100j,
                              self.model_config['gamma_min']:self.model_config['gamma_max']:100j,
                              self.model_config['strength_min']:self.model_config['strength_max']:100j
                              ]
        grid = griddata((a, g, s), mfit, (xi, yi, zi), method='linear')
    
    def grid_section(self, sec, axis):
        '''  
        Takes a 2-D section and makes a regularly spaced grid (for contouring)
        '''
        if axis not in ['alpha', 'gamma', 'strength']:
            raise ValueError(f'Axis {axis} undefined')
        elif axis == 'alpha':
            # alpha is fixed so slice is in gamma/str
            xi, yi = np.mgrid[self.model_config['strength_min']:self.model_config['strength_max']:2000j,
                              self.model_config['gamma_min']:self.model_config['gamma_max']:2000j]
            grid = griddata((sec.strength, sec.gamma), sec.misfit, (xi, yi), method='linear')
        elif axis == 'gamma':
            xi, yi = np.mgrid[self.model_config['strength_min']:self.model_config['strength_max']:2000j,
                              self.model_config['alpha_min']:self.model_config['alpha_max']:2000j]
            grid = griddata((sec.strength, sec.alpha), sec.misfit, (xi, yi), method='linear')       
        elif axis == 'strength':
            xi, yi = np.mgrid[self.model_config['alpha_min']:self.model_config['alpha_max']:2000j,
                              self.model_config['gamma_min']:self.model_config['gamma_max']:2000j]
            grid = griddata((sec.alpha, sec.gamma), sec.misfit, (xi, yi), method='linear')
            
        return xi, yi, grid
    
    def make_section(self, value, axis):
        '''
    
        '''
        if axis not in ['alpha', 'gamma', 'strength']:
            raise ValueError(f'Axis {axis} undefined')
        elif axis == 'alpha':
            cols = ['gamma', 'strength', 'misfit']
            assert np.less_equal(value, self.model_config['alpha_max'])
            assert np.greater_equal(value, self.model_config['alpha_min'])
        elif axis == 'gamma':
            cols = ['alpha', 'strength', 'misfit']
            assert np.less_equal(value, self.model_config['gamma_max'])
            assert np.greater_equal(value, self.model_config['gamma_min'])
        elif axis == 'strength':
            cols = ['alpha', 'gamma', 'misfit']
            assert np.less_equal(value, self.model_config['strength_max'])
            assert np.greater_equal(value, 0)
            
        slc = self.models[(np.isclose(self.models[axis], value, rtol=1.e-2)) & (self.models['misfit'] <= 1)]
        section = slc[cols]
        
        return section  

if __name__ == '__main__':
    ppv = Ensemble('/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/ppv_model/', strmax=0.5)
    ppv.make_section