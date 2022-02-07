#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:12:12 2020

@author: ja17375
"""
import numpy as np
import pandas as pd
from scipy import stats
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
from l2stats import ftest
FIG_DIR = '/Users/ja17375/Projects/Epac_fast_anom/Figures'


def restore_params(param, p_min, p_max):
    ''' 
    Restores true value of parameters in each model from their normalised state 
    
    Args:
        param : array-like
            normalised array of sampled values for a given parameter 
        p_min : float
            minimum (expanded) value of parameter (e.g. pmin=-180 for gamma)
        p_max : float
            maximum value for restored parameters
    Returns:
        true_params : array-like 
            array containing the sampled values restored to their true units.
    '''
    p_range = p_max - p_min 
    true_params = p_min + (p_range * param)
    
    return true_params

class Ensemble:
    '''
    A class to hold an ensemble of models produced by Matisse for a single anisotropic domain
    (this will get more complicated with multiple domains, which may required some
    retooling - i.e figuring out what column is what as the output Ensembles 
    have no column headings)
     
    Attributes
    -------------
        rundir : str
            path to local directory that contains results from the relevent MTS run
        dims : str
            number (and name) of dimensions used in inversion
        config : dict
            parameter configuration of inversion, stores minimum and maximum values for each parameter dimension
        models : DataFrame
            the full ensemble of models sampled by Matisse
        fcrit : float
            the F-test 95% confidence value for model misfit 
        best_fitting_model : DataFrame
            row selected from models which has the lowest misfit 
    Methods
    ----------
    read_ensemble(path, fname) :
        Reads in the raw model ensemble 
    find_fcrit(ndf) :
        Finds the f-test 95% confidence value for misfit values
    find_best_fitting(ret=False) :
        Selects the model from the ensemble with the lowest misfit value
    plot_2d_marginal_histograms(self, xparam, yparam, fname=None) : 
        Plots a 2-D histogram with 1-D marginals on the edge for a input pair of parameter
    plot_all_2d_histograms(fstem=None)
    
 
    '''
    def __init__(self, rundir, dims='ags', strmax=0.02, fname='MTS_ensemble.out'):
        '''
        Constructs EnsembleVisualiser. 
        Reads in raw ensemble data and converts it from normalised
        parameter vecotrs back to the real data units. 
        Configures in/outpiut directories and sets default values for different parameter combinations
        
        Parameters
        -----------
        rundir : str
            path to local directory that contains results from the relevent MTS run
        dims : str, optional, default='ags'
            number (and name) of dimensions used in inversion. Each parameter is string should
            be encoded by its first letter (e.g., alpha, gamma = 'ag')
        strmax : float, optional, default=0.02
            maximum value of strength parameter allowed in inversions
        fname : str, optional, default='MTS_ensemble.out'
            name of the file that contains the model ensemble output by Matisse
            
        '''
        self.config = {'alpha_min': 0, 'alpha_max': 90,
                       'beta_min': 0, 'beta_max': 90,
                       'gamma_min': -180,
                             'gamma_max': 180, 'strength_min': 0, 'strength_max': strmax}

        self.rundir = rundir
        self.dims = dims
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
            self.models = pd.DataFrame({'alpha': alpha,'gamma': gamma,
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
        elif dims == 'all':
            alpha = restore_params(raw_ensemble[:,0],
                                self.config['alpha_min'], self.config['alpha_max'])
            beta = restore_params(raw_ensemble[:,1], 
                                  self.config['beta_min'], self.config['beta_max'])
            gamma = restore_params(raw_ensemble[:,2],
                                self.config['gamma_min'], self.config['gamma_max'])
            strength = restore_params(raw_ensemble[:,3],
                                self.config['strength_min'], self.config['strength_max'])
            # raw ensemble also has our function f(x) that approximates (or is proportional to) the true PDF
            # evaluated for each model. f(x) = 1 / sum(lam2) for all paths in the inversion
            # to get the model misfit (i.e the sum of the lam2 residuals) we need to take the reciprocal
            misfit = 1 / raw_ensemble[:, -1]
            self.models = pd.DataFrame({'alpha': alpha, 
                                        'beta': beta, 
                                        'gamma': gamma,
                                       'strength':strength, 'misfit': misfit})
            print('Model Config set as:')
            print('{} < alpha < {}'.format(self.config['alpha_min'], self.config['alpha_max']))
            print('{} < beta < {}'.format(self.config['beta_min'], self.config['beta_max']))
            print('{} < gamma < {}'.format(self.config['gamma_min'], self.config['gamma_max']))
            print('{} < strength < {}'.format(self.config['strength_min'],
                                              self.config['strength_max']))
                      
    def read_ensemble(self, path, fname):
        ''' Reads in the raw model ensemble
        
        Parameters
        -----------
        path : str
            full path to model ensemble file
        fname : str
            name of model ensemble to read
            
        Returns 
        -----------
        raw_ensemble : array-like
            raw model ensemble in normalise parameter vector form
        '''
        raw_ensemble = np.loadtxt(f'{path}/{fname}')
        
        return raw_ensemble

    def find_fcrit(self, ndf):
        '''Find the f-test 95% confidence value  
        
        Parameters
        ------------
            ndf : int
                degrees of freedom of input data. The degrees of freedom for each waveform is calculated
                by Sheba for each path prior to inversions. The ndf for the dataset is the sum of these.
        '''
#       number of model dims is 1 less than cols of ensemble 
        k = self.models.shape[1] - 1 
        best = self.find_best_fitting(ret=True)
        # last col of model is always misfit
        self.fcrit = ftest(best[-1], ndf,k)
        
    def find_best_fitting(self, ret=False):
        '''
        Selects the model from the ensemble with the lowest misfit value.
        
        This model is assigned to self.best_fitting_model.
        
        Parameters
        ----------
        ret : boolean, optional
            Option for returning best fitting model. Default is False, which 
        '''
        imin = self.models.misfit.idxmin()
        self.best_fitting_model = self.models.iloc[imin]
        if ret:
            return self.best_fitting_model
        else:
            print('The best fitting model is:')
            print(f'Alpha = {self.best_fitting_model[0]:.5f}')
            print(f'Gamma = {self.best_fitting_model[1]:.5f}')
            print(f'Strength = {self.best_fitting_model[2]:.5f}')
            print(f'Model misfit = {self.best_fitting_model[3]:.5f}')
         
    def plot_2d_marginal_histograms(self, xparam, yparam, fname=None):
        '''
        Plots a 2-D histogram with 1-D marginals on the edge for a input pair of parameter
        
        Parameters
        ----------
        xparam : str
            name of model parameter to draw on x axis
        yparam : str
            name of model parameter to draw on y axis
        fname : str, optional, default=None
            filename to save figure as. if provided figure will be saved to the FIG_DIR
        '''
        params = ['alpha', 'gamma', 'beta', 'strength']
        if (xparam not in params) or (yparam not in params):
            raise ValueError(f'Unexpected parameters {xparam}, {yparam}')
        # fin best fitting model
        bestmodel = self.find_best_fitting(ret=True)
        xbest = bestmodel[xparam]
        ybest = bestmodel[yparam]
        # Make figure object
        ticks = {'alpha':[0,15,30,46,60,75,90],
                 'gamma':[-180,-120,-60,0,60,120,180],
                 'strength':np.linspace(0, self.config['strength_max'], 6)
                 }
            
        fig = plt.figure(figsize = (10,10))
        axs = gs(ncols = 5,nrows = 4, hspace =0.3, wspace=0.3,
                 width_ratios=[1,1,1,1,0.2], height_ratios=[1,1,1,1])
        # axis object for 1-D marginal histogram for yparam 
        ax_y = fig.add_subplot(axs[:-1,0])
        # axis object for 1-D marginal histogram for xparam
        ax_x = fig.add_subplot(axs[-1,1:-1])
        ax_main = fig.add_subplot(axs[:-1,1:-1],sharey=ax_y,sharex=ax_x)
        # axis obecjt for the colorbar
        ax_c = fig.add_subplot(axs[:-1,-1])
        
        # Plot 2-D histogram
        h2d = ax_main.hist2d(x=self.models[xparam], y=self.models[yparam], bins=[80, 80],
                             cmap='plasma')
        # add best fitting model[]
        ax_main.plot(xbest, ybest, 'kx',markersize=20)
        ax_main.set_xlim([self.config[f'{xparam}_min'], self.config[f'{xparam}_max']])
        ax_main.set_ylim([self.config[f'{yparam}_min'], self.config[f'{yparam}_max']])
        # Plot 1-D histogram for yparam
        ax_y.hist(self.models[yparam], bins=80,
                  orientation='horizontal')
        ax_y.set_ylabel(yparam.capitalize(), fontsize=14)
        ax_y.set_yticks(ticks[yparam])
        ax_y.invert_xaxis()
        
        # Plot 1-D histogram for xparam
        ax_x.hist(self.models[xparam], bins=80)
        ax_x.set_xticks(ticks[xparam])
        ax_x.invert_yaxis()
        ax_x.set_xlabel(xparam.capitalize(), fontsize=14)
        
        # Add colorbar
        plt.colorbar(h2d[3], cax=ax_c, aspect=20)
        plt.setp(ax_main.get_xticklabels(),visible=False)
        plt.setp(ax_main.get_yticklabels(),visible=False)
        
        if fname:
            plt.savefig(f'{FIG_DIR}/{fname}', dpi=500)
        else:
            plt.show()
        
    def plot_all_2d_histograms(self, fstem=None):
        '''
        Plot all 3 2d histograms using plot_2d_histograms
        if fstem is provided then figures are saved to FIGDIR
        
        Parameters
        ----------
        fstem : str, optional, defualt=None
            stem of filename to use when outputting figures
            
        '''
        
        if fstem:
            self.plot_2d_marginal_histograms('gamma', 'alpha',fname=f'{fstem}_gamma_alpha.png')
            self.plot_2d_marginal_histograms('gamma', 'strength', fname=f'{fstem}_gamma_str.png')
            self.plot_2d_marginal_histograms('alpha', 'strength', fname=f'{fstem}_alpha_str.png')
        else:
            self.plot_2d_marginal_histograms('gamma', 'alpha')
            self.plot_2d_marginal_histograms('gamma', 'strength')
            self.plot_2d_marginal_histograms('alpha', 'strength')
        
    def _add_2d_kde(self, ax, xparam, yparam):
        '''
        Makes a 2-D kernel desnity estimation plot for MTS Ensemble data
        
        Parameters:
        ----------
        ax : matplotlib.axes object
            axis handle to draw KDE onto
        xparam : str
            name of model parameter to draw on x axis
        yparam : str
            name of model parameter to draw on y axis
        
        '''
        # Get the correct x,y limits depending on input params
        xlim = [self.config[f'{xparam}_min'], self.config[f'{xparam}_max']]
        ylim = [self.config[f'{yparam}_min'], self.config[f'{yparam}_max']]
        # Set up samples 100 in both X and Y directions
        X, Y = np.mgrid[xlim[0]:xlim[1]:100j, ylim[0]:ylim[1]:100j]
        samples = np.vstack([X.ravel(), Y.ravel()])
        # Get params out of DF as numpy arrays
        x_models = self.models[xparam].values
        y_models = self.models[yparam].values
        # Stack params for KDE
        models = np.vstack([x_models, y_models])
        # 'make' KDE function
        kernal = stats.gaussian_kde(models)
        pdf = kernal(samples).T
        # Now we have our KDE, reshape the array for plotting
        Z = np.reshape(pdf, X.shape)
        ax.contourf(X, Y, Z, cmap=plt.cm.gist_earth_r, levels=8, vmin=1e-5)
        ax.contour(X, Y, Z, levels=10)
        # Add location of best fitting model
        x_best = self.best_fitting_model[xparam]
        y_best = self.best_fitting_model[yparam]
        ax.plot(x_best, y_best, 'x', color='red' )

    def _add_1d_kde(self, ax, param):
        '''
        Make a KDE plot for a set of model parameters
        
        Parameters:
        ----------
        ax : matplotlib.axes object
            axis handle to draw KDE onto
        param : str
            name of model parameter to draw on x axis
        '''
        lim = [self.config[f'{param}_min'], self.config[f'{param}_max']]
        x = np.linspace(lim[0], lim[1], 1000)
        ticks = np.linspace(lim[0], lim[1], 4) # Auto control ticks on axes labels
        values = self.models[param].values
        kernal = stats.gaussian_kde(values)
        pdf = kernal(x)
        ax.plot(x, pdf, 'k')
        # Add best fitting model
        ax.axvline(self.best_fitting_model[param], color='red', markersize=3)
        ax.set_xticks(ticks)
        ax.set_ylim([0, np.max(pdf)*1.1])
        
if __name__ == '__main__':
    ppv1 = Ensemble('/Users/ja17375/Projects/Epac_fast_anom/HQ_data/ScS_fix_test/ppv_001_100', strmax=0.5)
    ppv1.plot_all_2d_histograms('ppv1_2D_hist')
    ppv2 = Ensemble('/Users/ja17375/Projects/Epac_fast_anom/HQ_data/ScS_fix_test/ppv_010_100', strmax=0.5)
    ppv2.plot_all_2d_histograms('ppv2_2D_hist')
    # pv = Ensemble('/Users/ja17375/Projects/Epac_fast_anom/HQ_data/ScS_fix_test/pv_100_001', strmax=0.5)
    # ellipTI = Ensemble('/Users/ja17375/Projects/Epac_fast_anom/HQ_data/ScS_fix_test/ellipTI')
