#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 11:55:07 2020

@author: ja17375
"""
import EnsembleVisualiser
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import mplstereonet
from numpy import deg2rad, rad2deg, sin, cos
import numpy as np

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
    alpha = Ensemble.samples_3d[0,:]
    gamma = Ensemble.samples_3d[1,:]
    strength = Ensemble.samples_3d[2,:]
    ax.scatter(alpha, gamma, strength, c=Ensemble.pdf_3d)

def plot_pdf_slice(Ensemble, param, value, nsamps = 100):
    '''plot a 2d slice through a for a requested param. Interpolatre and resample 3-d PDF'''
    points, slc = Ensemble.slice_from_3d_pdf(param, value, nsamps)

    fig = plt.figure()
    if param == 'alpha':
        # Plot gamma v strength, polar plot
        gg = deg2rad(points[:,1].reshape(nsamps, nsamps))
        ss = points[:,2].reshape(nsamps, nsamps)
        ax = fig.add_subplot(111, projection='polar')
        ax.contourf(gg, ss, slc)
    elif param == 'gamma':
        # Plot alpha v strength, plot plot 
        aa = deg2rad(points[:,0].reshape(nsamps, nsamps))
        ss = points[:,2].reshape(nsamps, nsamps)
        ax = fig.add_subplot(111, projection='polar')
        ax.contourf(aa, ss, slc)
    elif param == 'strength':
        # Plot alpha and gamma, stereonet/ sterographic projection plot
        # Remember that alpha == dip and gamma == strike
        # To project dip onto polar plot using Lamber equal area:
        # R = 2*cos(alpha/2)   
        # alpha = 2*arccos(R/2)
        gg = deg2rad(points[:,1].reshape(nsamps, nsamps))
        a = points[:,0].reshape(nsamps, nsamps)
        aa = (2*cos(deg2rad(a)/2))
        ax = fig.add_subplot(111, projection='polar')
        ax.contourf(gg, aa, slc)
        
    else:
        raise ValueError(f'{param} not recognised')

    return ax

def plot_jelly_bowl(Ensemble):
    '''plot a spherical projection of PDF using strenght as radius and alpha, gamma as theta, phi'''
    
    alpha = Ensemble.samples_3d[0,:]
    gamma = Ensemble.samples_3d[1,:]
    strength = Ensemble.samples_3d[2,:]
    pdf = Ensemble.pdf_3d.reshape(alpha.shape) / Ensemble.pdf_3d.max()
    x, y, z = ags_to_xyz(alpha, gamma, strength)
    
    fig = go.Figure()
    fig.add_trace(go.Isosurface(
        x=x,
        y=y,
        z=z, 
        value=pdf,
        isomin=0.5,
        isomax=0.95
        ))
    fig.show()

def ags_to_xyz(alpha, gamma, s):
    '''Function that maps the natively "spherical" alpha, gamma, strength samples to
       cartesian co-ordinates for plotting
    '''
    radius = s 
    theta = deg2rad(alpha) 
    phi = deg2rad(gamma)
    
    x = radius * sin(theta) * cos(phi)
    y = radius * sin(theta) * sin(phi)
    z = radius * cos(theta)
    
    return x, y, z

def slice_grid(Ensemble, nsamps, axis):
    '''
    slices 3D PDF from Ensemble alogn the desired axis, returning a 2D slice
    '''

def plot_2d_slice_through_str(Ensemble, s, nsamps):
    '''
    plot a 2-d slice throught the model ensemble with a fixed strength.
    i.e plot likelihoods for alpha, gamma for a given stength param
    '''
    #

def plot_2d_slice_through_alpha(Ensemble, alpha, nsamps):
    '''
    plot a 2-d slice through model PDF 
    '''
    alpha_samples = alpha.reshape(50,50,50)[:,0,0]
    ax = plt.subplot(111, projection='polar')

def plot_3d_isosurface(Ensemble):
    
    alpha = Ensemble.samples_3d[0,:]
    gamma = Ensemble.samples_3d[1,:]
    strength = Ensemble.samples_3d[2,:]
    pdf = Ensemble.pdf_3d.reshape(alpha.shape) / Ensemble.pdf_3d.max()
    # "Normalise" the PDF
    fig= go.Figure() 
    fig.update_layout(
                    xaxis_title_text='alpha  (dip of TI plane)',
                    yaxis_title_text='gamma (azi of TI plane)',
                    

                    title_text = '3-D PDF' 
                    )
    fig.add_trace(go.Isosurface(
    x=alpha,
    y=gamma,
    z=strength,
    value=pdf,
    cmin=0,
    cmax=1,
    isomin=0.7, 
    isomax=0.95,
    opacity=0.5,
    surface_count=2,
    caps=dict(x_show=False, y_show=False, z_show=False),
    ))

    
    fig.add_trace(go.Isosurface(
    x=alpha,
    y=gamma,
    z=strength,
    value=pdf,
    cmin=0,
    cmax=1,
    isomin=0, 
    isomax=1,
    surface_count=2,
    caps=dict(x_show=False, y_show=False, z_show=False),
    slices_z=dict(show=True, locations=[0.0075]),
    slices_x=dict(show=True, locations=[20, 75])
    ))
    

    fig.show()

if __name__ == '__main__':

    mode = input('Enter Mode to run in [1] plot or [2] gen and save PDF > ')
    if mode == '1':
        Ensemble = EnsembleVisualiser.Ensemble(
            '/Users/ja17375/SWSTomo/BlueCrystal/Dom1160/StackedCorrs/Alpha0_90/AddingDAN',
            read=False)
        Ensemble.read_3d_pdf('/Users/ja17375/SWSTomo/BlueCrystal/Dom1160/StackedCorrs/Alpha0_90/AddingDAN/MTS_3D_pdf.npy')
       # Ensemble.evaluate_3d_kde(nsamps=100j)
        plot_3d_isosurface(Ensemble)
        #plot_jelly_bowl(Ensemble)
        
    elif mode =='2':
        Ensemble = EnsembleVisualiser.Ensemble(
        '/Users/ja17375/SWSTomo/BlueCrystal/Dom1160/StackedCorrs/Alpha0_90/Free_ags',
        'MTS_ensemble.out')
        Ensemble.evaluate_3d_kde(nsamps=50j)
        Ensemble.save_3d_pdf('MTS_3D_pdf.npy')
        print('Saved 3D PDF')

        