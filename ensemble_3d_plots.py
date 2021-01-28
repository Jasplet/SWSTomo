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
from skimage import measure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
FIG_DIR = '/Users/ja17375/SWSTomo/Figures'

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

    alpha, gamma, strength = np.meshgrid(Ensemble.alpha_samples,
                                     Ensemble.gamma_samples,
                                     Ensemble.strength_samples)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(alpha, gamma, strength, c=Ensemble.pdf_3d)
    plt.show()

def plot_model_slices(Ensemble, model, method, nsamps = 100, save=False):
    '''plot a 2d slicethrough PDF for a inout model point'''
    # Check model
    if len(model) == 3:
        print('Plotting 2-D sections through model')
        print(f'alpha = {model[0]}')
        print(f'gamma = {model[1]}')
        print(f'strength = {model[2]}')
    elif model == 'most-likely':
        Ensemble.find_most_likely()
        model = Ensemble.most_likely_model
        print('Plotting ...')
    else:
        raise ValueError('Model not understood')
    
    fig = plt.figure(figsize= (16, 8))        
    # Plot gamma v strength, polar plot (Alpha is fixed to mode value)
    points, slc = Ensemble.slice_3d_volume('alpha', model[0], nsamps, method)
    gg = deg2rad(points[:,1].reshape(nsamps, nsamps))
    ss = points[:,2].reshape(nsamps, nsamps)
    ax1 = fig.add_subplot(131, projection='polar')
    ax1.plot(deg2rad(model[1]), model[2], 'kx')
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)
    C1 = ax1.contour(gg, ss, slc, 10)
    ax1.clabel(C1, C1.levels, inline=True, fontsize=10)
    ax1.set_title(f'Slice through alpha axis. alpha = {model[0]:0.3f}')

    # Plot alpha v strength, polar plot (Gamma is fixed)
    points, slc = Ensemble.slice_3d_volume('gamma', model[1], nsamps, method)
    aa = deg2rad(points[:,0].reshape(nsamps, nsamps))
    ss = points[:,2].reshape(nsamps, nsamps)
    ax2 = fig.add_subplot(132, projection='polar')
    C2 = ax2.contour(aa, ss, slc, 10)
    ax2.clabel(C2, C2.levels, inline=True, fontsize=10)
    ax2.plot(deg2rad(model[0]), model[2], 'kx')
    ax2.set_theta_zero_location('N')
    ax2.set_theta_direction(-1)
    ax2.set_thetamin(0)
    ax2.set_thetamax(90)
    ax2.set_title(f'Slice through Gamma axis, gamma = {model[1]:0.3f}')
 
    # Plot alpha and gamma, stereonet/ sterographic projection plot
    # Remember that alpha == dip and gamma == strike
    # To project dip onto polar plot using Lambert equal area:
    # R = 2*cos(alpha/2)   
    # alpha = 2*arccos(R/2)
    points, slc = Ensemble.slice_3d_volume('strength', model[2], nsamps, method)
    gg = deg2rad(points[:,1].reshape(nsamps, nsamps))
    a = points[:,0].reshape(nsamps, nsamps)
    aa = (2*cos(deg2rad(a)/2))
    ax3 = fig.add_subplot(133, projection='polar')
    C3 = ax3.contour(gg, aa, slc, 10)
    ax3.clabel(C3, C3.levels, inline=True, fontsize=10)
    ax3.plot(deg2rad(model[1]), (2*cos(deg2rad(model[0])/2)), 'kx')

    ax3.set_theta_zero_location('N')
    ax3.set_theta_direction(-1)
    ax3.set_title(f'Slice through strength axis, s = {model[2]:0.3f}')
    
    if save is True:
        fname = input('Enter filename (without extension) to save as >')
        plt.savefig(f'{FIG_DIR}/{fname}.png', dpi=500)
    plt.show()


def ags_to_xyz(alpha, gamma, s):
    '''Function that maps the natively "spherical" alpha, gamma, strength samples to
       cartesian co-ordinates for plotting
    '''
    radius = s 
    theta = deg2rad(90 - alpha) 
    phi = deg2rad(gamma) 
    x = radius * sin(theta) * cos(phi)
    y = radius * sin(theta) * sin(phi)
    z = radius * cos(theta) 
    return x, y, z
    
def isosurface_view(Ensemble, level):
    '''
    Uses the marching cubes method to find an isosurface in the 3D PDF.
    This isosurface is plotted using matplotlib's 3D toolkit 
    
    '''
    # Normalise PDF to max likelihood
    pdf = Ensemble.pdf_3d / Ensemble.pdf_3d.max()
    a_spacing = Ensemble.alpha_samples[1] - Ensemble.alpha_samples[0]
    g_spacing = Ensemble.gamma_samples[1] - Ensemble.gamma_samples[0]
    s_spacing = Ensemble.strength_samples[1] - Ensemble.strength_samples[0]
    space = (a_spacing, g_spacing, s_spacing)
    verts, faces, normals, values = measure.marching_cubes(pdf, level, spacing=space)
    verts[:,1] -= 180 # Map gamma verticies back to -180 - 180
    xx, yy, zz = ags_to_xyz(verts[:,0], verts[:,1], verts[:,2])
    sphe_verts = np.vstack((xx, yy, zz)).T
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    mesh = Poly3DCollection(sphe_verts[faces])
    mesh.set_edgecolor('k')
    ax.add_collection3d(mesh)
    add_str_scale(ax)
    # Add the most likely model 
    model = Ensemble.find_most_likely(ret=True)
    xm, ym, zm = ags_to_xyz(model[0], model[1], model[2])
    ax.plot([0, xm], [0, ym], [0, zm],'-')
    ax.plot(xm, ym, zm, 'x')
    # ax.scatter(xx, yy, zz, c=values)
    ax.set_xlim([-0.02, 0.02])
    ax.set_xlabel('X')
    ax.set_ylim([-0.02, 0.02])
    ax.set_ylabel('Y')
    ax.set_zlim([0, 0.02])
    # Set initial view angle to be top down and so that 0 degrees points North
    ax.view_init(90, 180)
    plt.savefig(f'test_isosurface_view_level{level}.png', dpi=400)
    plt.show()

def cartesian_isosurface_view(Ensemble, level, save=False):
    '''
    '''
    # Normalise PDF to max likelihood
    pdf = Ensemble.pdf_3d / Ensemble.pdf_3d.max()
    a_spacing = Ensemble.alpha_samples[1] - Ensemble.alpha_samples[0]
    g_spacing = Ensemble.gamma_samples[1] - Ensemble.gamma_samples[0]
    s_spacing = Ensemble.strength_samples[1] - Ensemble.strength_samples[0]
    space = (a_spacing, g_spacing, s_spacing)
    verts, faces, normals, values = measure.marching_cubes(pdf, level, spacing=space)
    verts[:,1] -= 180 # Map gamma verticies back to -180 - 180 
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    mesh = Poly3DCollection(verts[faces])
    mesh.set_edgecolor('k')
    ax.add_collection3d(mesh)
    # Add the most likely model 
    model = Ensemble.find_most_likely(ret=True)
    ax.plot([0, model[0]], [0, model[1]], [0, model[2]],'-')
    ax.plot(model[0], model[1], model[2], 'x')    
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Gamma')
    ax.set_zlabel('Strength')
    ax.set_xlim([0, 90])
    ax.set_ylim([-180, 180])
    ax.set_zlim([0, 0.02])
    
    if save is True:
        uid = input('Enter a inversion uid to add to fiename')
        plt.savefig(f'{FIG_DIR}/{uid}_cart_isosurf_lvl_{level}.png', dpi=500)
    plt.show()

def add_str_scale(ax):
    '''
    draws circles on 3-D pyplot axis to show strength parameter (radius of PDF)
    '''
    t = np.linspace(0, 2*np.pi, 50)
    strs = [0.005, 0.01, 0.015, 0.02] # strenght params to use
    for s in strs:
        x = s*np.cos(t)
        y = s*np.sin(t)
        z = np.zeros(50)
        ax.plot(x, y, z,'r-')
        ax.text(s*cos(np.pi/4), s*sin(np.pi/4),0, s=str(s))        
    ax.plot([0, 0.02], [0, 0], [0, 0], 'k-')

def plot_3d_isosurface(Ensemble):
    '''
    Uses the plotly package to draw a more interactive isosuface view of the 
    3-D PDF 
    '''
    alpha, gamma, strength = np.meshgrid(Ensemble.alpha_samples,
                                         Ensemble.gamma_samples,
                                         Ensemble.strength_samples)
    pdf = Ensemble.pdf_3d.flatten() / Ensemble.pdf_3d.max()
    # "Normalise" the PDF
    fig= go.Figure() 
    fig.update_layout(
                    xaxis_title_text='alpha  (dip of TI plane)',
                    yaxis_title_text='gamma (azi of TI plane)',
                    title_text = '3-D PDF' 
                    )
    fig.add_trace(go.Isosurface(
    x=alpha.flatten(),
    y=gamma.flatten(),
    z=strength.flatten(),
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

def plot_param_histograms(Ensemble):
    '''
    makes a histogram for each model parameter 
    '''
    models = Ensemble.models
    fig = plt.figure(figsize=(21,7))
    # historgram for alpha
    abins = np.linspace(0, 91, 50)
    ax1 = fig.add_subplot(131)
    ax1.hist(models.alpha, bins=abins)
    ax1.set_xlim([0, 90])
    ax1.set_xlabel(r'alpha ($\degree$)')
    # histogram for gamma
    gbins = np.linspace(-180, 181, 50)
    ax2 = fig.add_subplot(132)
    ax2.hist(models.gamma, bins=gbins)
    ax2.set_xlim([-180, 180])
    ax2.set_xlabel(r'gamma ($\degree$)')
    # histogram for strength
    sbins = np.linspace(0, 0.02, 50)
    ax3 = fig.add_subplot(133)
    ax3.hist(models.strength, bins=sbins)
    ax3.set_xlim([0, 0.02])
    ax3.set_xlabel('strength parameter')

    plt.show()
   
def plot_param_v_misfit(Ensemble):
    '''
    
    '''    
    models = Ensemble.models
    fig = plt.figure(figsize=(21,7))
    ax1 = fig.add_subplot(131)
    ax1.scatter(x=models.alpha, y=models.misfit)
    # gamma
    ax2 = fig.add_subplot(132)
    ax2.scatter(x=models.gamma, y=models.misfit)
    # strength
    ax3 = fig.add_subplot(133)
    ax3.scatter(x=models.strength, y=models.misfit)
    
    plt.show()

def plot_nbest_models(Ensemble, n, full_range=False):
    
    models_sort = Ensemble.models.sort_values(by=['misfit'])
    pmodels = models_sort.iloc[0:n]
    # sort again so lowest msfits are at the end and therefore drawn on top!
    models = pmodels.sort_values(by=['misfit'], ascending=False)
    fig1 = plot_3d_nbest_models(models, n)
    (ax_3d) = fig1.axes[0]
    fig2 = plot_2d_views_nbest_model(models, n)
    (ax1, ax2, ax3, cbar) = fig2.get_axes()
    
    if full_range:
        ax_3d.set_xlim([Ensemble.model_config['alpha_min'],
                        Ensemble.model_config['alpha_max']])
        ax_3d.set_ylim([Ensemble.model_config['gamma_min'],
                        Ensemble.model_config['gamma_max']])
        ax_3d.set_zlim([Ensemble.model_config['strength_min'],
                        Ensemble.model_config['strength_max']])
        
        ax1.set_xlim([Ensemble.model_config['strength_min'],
                      Ensemble.model_config['strength_max']])
        ax1.set_ylim([Ensemble.model_config['alpha_min'],
                  Ensemble.model_config['alpha_max']])
        ax2.set_xlim([Ensemble.model_config['strength_min'],
                  Ensemble.model_config['strength_max']])
        ax2.set_ylim([Ensemble.model_config['gamma_min'],
                  Ensemble.model_config['gamma_max']]) 
        ax3.set_xlim([Ensemble.model_config['alpha_min'],
          Ensemble.model_config['alpha_max']]) 
        ax3.set_ylim([Ensemble.model_config['gamma_min'],
          Ensemble.model_config['gamma_max']]) 
    # Add best fitting model to both figs
    best_model = Ensemble.find_best_fitting(ret=True)
    ax_3d.plot(xs=best_model.alpha, ys=best_model.gamma,
                 zs=best_model.strength, marker='x', color='red')
    ax1.plot(best_model.strength, best_model.alpha, marker='x', color='red')
    ax2.plot(best_model.strength, best_model.gamma, marker='x', color='red')
    ax3.plot(best_model.alpha, best_model.gamma, marker='x', color='red')
    plt.show()
    
def plot_3d_nbest_models(models, n):
    '''
    Makes a 3-D (x,y,z) plot showing the n best fitting models
    Args:
        Ensemble (Class) - The ensemble object holding all the models 
        n (float) - number of models to plot
    Returns:
        A figure object containgng the 3-D plot
    '''
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    C = ax.scatter(models.alpha, models.gamma, models.strength, c=models.misfit,s=6, cmap='viridis_r')
    ax.set_xlabel(r'$\alpha \degree $')
    ax.set_ylabel(r'$\gamma \degree $')    
    ax.set_zlabel('strength')
    cbar = fig.colorbar(C, ax=ax)
    cbar.ax.set_ylabel('Model misfit')
    return fig

def plot_2d_views_nbest_model(models, n):
    ''' 
    Makes a 3-panel figure showing 2-D views of the n best fitting models
    '''
    fig = plt.figure(figsize=(20, 6))
    # plot alpha v strength
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.scatter(x=models.strength, y=models.alpha, c=models.misfit, s=6, cmap='viridis_r')
    ax1.set_xlabel('Strength')
    ax1.set_ylabel(r'Alpha ($\degree$)')
    # plot gammav strength
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.scatter(x=models.strength, y=models.gamma, c=models.misfit, s=6, cmap='viridis_r')
    ax2.set_xlabel('Strength')
    ax2.set_ylabel(r'Gamma ($\degree$)')
    # plot gamma v alpha
    ax3 = fig.add_subplot(1, 3, 3)
    C = ax3.scatter(x=models.alpha, y=models.gamma, c=models.misfit, s=6, cmap='viridis_r')
    ax3.set_xlabel('Alpha ($\degree$)')
    ax3.set_ylabel('Gamma ($\degree$)')
    cbar = fig.colorbar(C)
    cbar.ax.set_ylabel(r'$\Sigma{\lambda_2]$')
    return fig

def plot_3d_best_misfit(Ensemble, mf):
    '''
    Makes a 3-D (a, g, s) plot of models from Ensemble with misfit less than or equal to mf
    '''
    models = Ensemble.models[Ensemble.models.misfit <= mf]
    n = len(models)
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    C = ax.scatter(xs=models.alpha,
                   ys=models.gamma,
                   zs=models.strength,)
    ax.set_xlabel(r'$\alpha \degree $')
    ax.set_ylabel(r'$\gamma \degree $')    
    ax.set_zlabel('strength')
    ax.set_xlim([0,90])
    ax.set_ylim([-180,180])
    ax.set_title(f'The {n} models with misfit <= {mf:f}')
    cbar = fig.colorbar(C, ax=ax)
    cbar.ax.set_ylabel('Model misfit')
    plt.show()    

def plot_3d_spherical_model_misfit(Ensemble, n):
    '''
    Makes a 3-D (spherical) plot showing the n best fitting models
    Args:
        Ensemble (Class) - The ensemble object holding all the models 
        n (float) - number of models to plot
    Returns:
        A figure object containgng the 3-D plot
    '''
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    models_sort = Ensemble.models.sort_values(by=['misfit'],ascending=True)
    pmodels = models_sort.iloc[0:n]
    x, y, z = ags_to_xyz(pmodels.alpha, pmodels.gamma, pmodels.strength)
    C = ax.scatter(x, y, z, c=pmodels.misfit,s=6)
    ax.set_xlabel(r'$\alpha \degree $')
    ax.set_ylabel(r'$\gamma \degree $')    
    ax.set_zlabel('strength')
    cbar = fig.colorbar(C, ax=ax)
    cbar.ax.set_ylabel('Model misfit')
    plt.show()

if __name__ == '__main__':

    mode = input('Enter Mode to run in [1] plot or [2] gen and save PDF > ')
    if mode == '1':
        Ensemble = EnsembleVisualiser.Ensemble(
            '/Users/ja17375/SWSTomo/BlueCrystal/Dom1160/StackedCorrs/Alpha0_90/AddingDAN',
            read=False)
        Ensemble.read_3d_pdf('/Users/ja17375/SWSTomo/BlueCrystal/Dom1160/StackedCorrs/Alpha0_90/AddingDAN/D1160_doubleScS_pathlength_PDF.npy')
       # Ensemble.evaluate_3d_kde(nsamps=100j)
        model = [5.1138, 139.5269, 0.0089]
        isosurface_view(Ensemble, 0.5,model)
        #plot_jelly_bowl(Ensemble)
        plt.show()
        
    elif mode =='2':
        Ensemble = EnsembleVisualiser.Ensemble(
        '/Users/ja17375/SWSTomo/BlueCrystal/Dom1160/StackedCorrs/Alpha0_90/Free_ags',
        'MTS_ensemble.out')
        Ensemble.evaluate_3d_kde(nsamps=50j)
        Ensemble.save_3d_pdf('MTS_3D_pdf.npy')
        print('Saved 3D PDF')

        