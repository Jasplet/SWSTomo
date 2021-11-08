#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:12:12 2020

@author: ja17375
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import griddata, RegularGridInterpolator
import matplotlib.pyplot as plt
from l2stats import ftest
FIG_DIR = '/Users/ja17375/SWSTomo/Figures'


def restore_params(param, p_min, p_max):
    ''' 
    Restores true value of parameters in each model from their normalised state 
    '''
    p_range = p_max - p_min 
    true_params = p_min + (p_range * param)
    
    return true_params

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
        self.config = {'alpha_min': 0, 'alpha_max': 90, 'gamma_min': -180,
                             'gamma_max': 180, 'strength_min': 0, 'strength_max': strmax}

        self.rundir = rundir
        if read:
            raw_ensemble = self.read_ensemble(rundir,fname)
            if dims == 'ags':
                alpha = restore_params(raw_ensemble[:,0],
                                    self.config['alpha_min'], self.config['alpha_max'])
                gamma = restore_params(raw_ensemble[:,1],
                                    self.config['gamma_min'], self.config['gamma_max'])
                strength = restore_params(raw_ensemble[:,2],
                                    self.config['strength_min'], self.config['strength_max'])
                # raw ensemble also has our function f(x) that approximates (or is proportional to) the true PDF
                # evaluated for each model. f(x) = 1 / sum(lam2) for all paths in the inversion
                # to get the model misfit (i.e the sum of the lam2 residuals) we need to take the reciprocal
                misfit = 1 / raw_ensemble[:, -1]
                self.models = pd.DataFrame({'alpha': alpha, 'gamma': gamma,
                                           'strength':strength, 'misfit': misfit})
                print('Model Config set as:')
                print('{} < alpha < {}'.format(self.config['alpha_min'], self.config['alpha_max']))
                print('Beta = fixed')
                print('{} < gamma < {}'.format(self.config['gamma_min'], self.config['gamma_max']))
                print('{} < strength < {}'.format(self.config['strength_min'],
                                                  self.config['strength_max']))
                      
            elif dims == 'ag':
                #2-D alpha and gamma only
                alpha = self.restore_params(raw_ensemble[:,0],
                    self.config['alpha_min'], self.config['alpha_max'])
                gamma = self.restore_params(raw_ensemble[:,1],
                                    self.config['gamma_min'], self.config['gamma_max'])
                misfit = 1 / raw_ensemble[:, -1]
                self.models = pd.DataFrame({'alpha': alpha, 'gamma': gamma,
                                            'misfit': misfit})
            elif dims == 'as':
                alpha = self.restore_params(raw_ensemble[:,0],
                    self.config['alpha_min'], self.config['alpha_max'])
                strength = self.restore_params(raw_ensemble[:,1],
                                    self.config['strength_min'], self.config['strength_max'])
                misfit = 1 / raw_ensemble[:, -1]
                self.models = pd.DataFrame({'alpha': alpha, 'strength': strength,
                                           'misfit': misfit})
            elif dims == 'gs':
                gamma = self.restore_params(raw_ensemble[:,0],
                                    self.config['gamma_min'], self.config['gamma_max'])
                strength = self.restore_params(raw_ensemble[:,1],
                                    self.config['strength_min'], self.config['strength_max'])
                misfit = 1 / raw_ensemble[:, -1]
                self.models = pd.DataFrame({'gamma': gamma,
                                           'strength':strength, 'misfit': misfit})
        # Store modes as a DataFrame for ease of reference. Misfit index as [:,-1] in order to
        # "future proof" as misfit is always last column of MTS_ensemble.out 

    def read_ensemble(self, path, fname):
        '''function to read the MTS_emsemble from the input run directory'''
        raw_ensemble = np.loadtxt(f'{path}/{fname}')
        
        return raw_ensemble



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
        self.alpha_samples = np.linspace(self.config['alpha_min'],
                                         self.config['alpha_max'],nsamps)
        self.gamma_samples = np.linspace(self.config['gamma_min'],
                                         self.config['gamma_max'],nsamps)
        self.strength_samples = np.linspace(self.config['strength_min'],
                                            self.config['strength_max'],nsamps)
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
        self.alpha_samples = np.linspace(self.config['alpha_min'],
                                         self.config['alpha_max'],self.nsamps)
        self.gamma_samples = np.linspace(self.config['gamma_min'],
                                         self.config['gamma_max'],self.nsamps)
        self.strength_samples = np.linspace(self.config['strength_min'],
                                            self.config['strength_max'],self.nsamps)
     
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
        
    def save_3d_pdf(self, outfile):
        ''' Saves 3-D PDF as Numpy Array'''
        np.save(f'{self.rundir}/{outfile}', self.pdf_3d)

    def param_kde(self, ax, p, param):
        '''
        Make a Kernel Density estimate of the parameter and plot on fig axis
        '''
        kde = stats.gaussian_kde(p)
        if param in ['alpha','gamma','strength']:
            mx = self.config[f'{param}_max']
            mn = self.config[f'{param}_min']
        elif param == 'misfit':
            mx = 1
            mn = 0 
        else:
            raise ValueError
        
        xx = np.linspace(mn, mx, 1000)
        y = kde(xx)
        ax.plot(xx, y)
        ax.set_xlim([mn, mx])
        ax.set_ylim([0, y.max()*1.15])
        return y.max()
    
    def plot_param_hists(self,candidate):
        '''
        Makes a 4-panel histogram of model parameters and misfit
        '''  
        if candidate == 'ellip':
            title = 'Elliptical TI'
            fout = 'ellip_TI'
        elif candidate == 'br':
            title = 'Bridgmanite (100)[001]'
            fout = 'br_100_001'
            bestmodel = [42.85864, 22.26776,
                         0.25486, 0.42735]
        elif candidate == 'ppv1':
            title = 'Post-Perovskite (001)[100]'
            fout = 'ppv_001_100'
            bestmodel = [22.80220, -84.61130,
                         0.10352, 0.44248]
        elif candidate == 'ppv2':
            title = 'Post-Perovskite (010)[100]'
            fout = 'ppv_010_100'
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(7,7))
        
        # ax1.hist(self.models.alpha,bins=20)
        a_kde = self.param_kde(ax1, self.models.alpha, 'alpha')
        ax1.vlines(bestmodel[0], ymin=0, ymax=a_kde*1.15)
        ax1.set_title('Alpha')
        ax1.set_ylabel('Density')
        ax1.set_xticks([0, 30, 60, 90])
        # ax2.hist(self.models.gamma,bins=20)
        g_kde = self.param_kde(ax2, self.models.gamma, 'gamma')
        ax2.vlines(bestmodel[1], ymin=0, ymax=g_kde*1.15)
        ax2.set_title('Gamma') 
        ax2.set_ylabel('Density')
        ax2.set_xticks([-180,-120,-60,0,60,120,180])
        # ax3.hist(self.models.strength,bins=20)
        s_kde = self.param_kde(ax3, self.models.strength, 'strength')
        ax3.vlines(bestmodel[2], ymin=0, ymax=s_kde*1.15)
        ax3.set_title('Strength')
        ax3.set_ylabel('Density')
        # ax4.hist(self.models.misfit,bins=20)
        m_kde = self.param_kde(ax4, self.models.misfit, 'misfit')
        ax4.vlines(bestmodel[3], ymin=0, ymax=m_kde*1.15)
        ax4.set_title('Model Misfit')
        ax4.set_ylabel('Density')
        fig.suptitle(f'{title} bootstrapped KDEs')
        plt.tight_layout()
        plt.savefig(f'{FIG_DIR}/{fout}_bootstrap_KDE.png',dpi=400)
        plt.show()
if __name__ == '__main__':
    ppv = Ensemble('/Users/ja17375/SWSTomo/Inversions/Epac_fast_anom/ppv_model/', strmax=0.5)
    ppv.plot_2d_sections(239)