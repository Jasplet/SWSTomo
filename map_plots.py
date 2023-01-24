#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:03:35 2021

@author: ja17375
"""

import matplotlib.pyplot as plt
import pygmt
import numpy as np
import pandas as pd
# import cartopy.feature as cfeature 
# import cartopy.crs as ccrs
#from plot_models import draw_trigonal_doms
from pathset import find_phases_in_domain
from sphe_trig import vincenty_dist

FIG_DIR = '/Users/ja17375/Projects/Epac_fast_anom/Figures'
THESIS_FIGS = '/Users/ja17375/Thesis-enclosing/Thesis/chapters/chapter02/Figs'
def map_data_at_station(phasefile,station,extent=[-150,-100,30,70],save=False):
    '''
    This function maps all the phases recorded at a given station. 
    This allows us to assess how good the coverage is (or isnt!)
    
    Args:
        phasefile [str] - path to the dataset being used
        station [str] - station to make the plot for
        save [bool] - switch for if you want to save the plot or just show it
        
    Returns:
        A map showing the input station and all phases (in the datafile) recorded at it,
        on the assumption they are all used to invert for the station correction
    '''
    domains = np.loadtxt('T3_global.bins',skiprows=1)
    dom1160 = domains[np.isin(domains[:,0],1160)] # Hard coded for now as we are only looking at D1160
    alldata = pd.read_csv(phasefile,delim_whitespace=True)
    data = alldata[alldata.STAT == station]
    scs = data[data.PHASE == 'ScS']
    sks = data[data.PHASE == 'SKS']
    skks = data[data.PHASE == 'SKKS']
    print(alldata)
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto')) #coastline data
    draw_trigonal_doms(ax, dom1160)

    ax.plot(data.STLO, data.STLA, linestyle='', marker = 'v', markersize=14, 
       transform=proj, color='red', markeredgecolor='black', label='Stations')
    ax.plot(scs.LOWMM_LON, scs.LOWMM_LAT, markeredgecolor='black',
            linestyle='',color='cornflowerblue', marker='o', transform=proj,
            label='ScS', markersize=8)
    ax.plot(sks.LOWMM_LON, sks.LOWMM_LAT, markeredgecolor='black',
            linestyle='', color='orange', marker='s', transform=proj,
            label='SKS', markersize=8)
    ax.plot(skks.LOWMM_LON, skks.LOWMM_LAT, markeredgecolor='black',
            linestyle='', color='mediumseagreen', marker='p', transform=proj,
            label='SKKS', markersize=8)
    ax.legend()
    ax.gridlines(draw_labels=True,linewidth=0.5)
    
    if save:
        plt.savefig(f'{FIG_DIR}/{station}_SC_phasemap.png',bbox_inches='tight')
        
def map_single_domain_phases(phasefile,dom_ID,stations,extent=[110,260,-10,65]):
    '''
    This function draws a map for a single domain (or potentially multiple domains in the future) along with the phase locations in D`` and the station locations
    '''
    domains = np.loadtxt('T3_global.bins',skiprows=1)
    doms2plot = domains[np.isin(domains[:,0],dom_ID)]
    data = find_phases_in_domain(phasefile,dom_ID)
    scs = data[data.PHASE == 'ScS']
    sks = data[data.PHASE == 'SKS']
    skks = data[data.PHASE == 'SKKS']
    stats = data[data.STAT.isin(stations)]
    # Draw figure
    proj0 = ccrs.PlateCarree(central_longitude=0)
    proj180 = ccrs.PlateCarree(central_longitude=180)
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111,projection=proj180)
    ax.set_extent(extent, crs=proj0)      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto')) #coastline data
    draw_trigonal_doms(ax, doms2plot)

    for i, stat in stats.iterrows():
        ax.text(stat.STLO+2, stat.STLA, stat.STAT, fontsize=12)
        ax.plot([stat.EVLO, stat.STLO], [stat.EVLA, stat.STLA], 'k-',
           transform=ccrs.Geodetic()) 
    ax.plot(data.EVLO, data.EVLA, linestyle='', marker='*', markersize=12,
            transform=proj0, color='blue')

    ax.plot(scs.LOWMM_LON, scs.LOWMM_LAT, markeredgecolor='black',
            linestyle='',color='cornflowerblue', marker='o', transform=proj0,
            label='ScS', markersize=8)
    ax.plot(sks.LOWMM_LON, sks.LOWMM_LAT, markeredgecolor='black',
            linestyle='', color='orange', marker='s', transform=proj0,
            label='SKS', markersize=8)
    ax.plot(skks.LOWMM_LON, skks.LOWMM_LAT, markeredgecolor='black',
            linestyle='', color='mediumseagreen', marker='p', transform=proj0,
            label='SKKS', markersize=8)
    ax.plot(data.STLO, data.STLA, linestyle='', marker = 'v', markersize=9,
                transform=proj0, color='grey')
    
    ax.plot(stats.STLO, stats.STLA, linestyle='', marker = 'v', markersize=14, 
       transform=proj0, color='red', markeredgecolor='black')
    # ax.text(-151,67, '2 ScS phases'.format(len(scs)))
    # ax.text(-151,65.75, '2 SKS phases'.format(len(sks)))
    # ax.text(-151,64.5, '3 SKKS phases'.format(len(skks)))
    # ax.set_title(r"Phases passing through domain Lower {}".format(dom_ID)))
    ax.set_title('Data used in inversions')
    ax.legend()
    grd = ax.gridlines(draw_labels=True,linewidth=0.5)
    grd.top_labels = None
    plt.show()
    #plt.savefig(f'{FIG_DIR}/Phases_in_D{dom_ID}.png',format='png',dpi=400,
                # bbox_inches='tight')
                
def map_data_tomo(data, region, draw_paths=False, fname=None):

    hqdata = data[(data.SNR >= 10.) & ((data.Q >= 0.7) | (data.Q <= -0.9))]
    scs = hqdata[hqdata.PHASE == 'ScS']
    sks = hqdata[hqdata.PHASE == 'SKS']
    skks = hqdata[hqdata.PHASE == 'SKKS']
    
    fig = pygmt.Figure()
    fig.basemap(region=region, projection='M15c', frame=['a10','WESN'])

    cpt = f'{FIG_DIR}/S40RTS/S40RTS.cpt'
    s40rts = f'{FIG_DIR}/S40RTS/S40RTS_2800km.grd'
    fig.grdimage(grid=s40rts, cmap=cpt)
    fig.colorbar(cmap=cpt,
                position='jBC+o0c/-1.5c/+w8c/0.5c+h+e+m',
                frame=['a1f0.5g0.25','x+l"dVs (%)"'])
    fig.coast(shorelines='1/0.5p,black', resolution='l')        
    # fig.plot(x=d1160_vlon,y=d1160_vlat,pen='2p,gray30,dashed')
    if draw_paths:

        for i, row in hqdata.iterrows():
            fig.plot(x=[row.EVLO, row.STLO], y=[row.EVLA, row.STLA],pen='0.75p,black')
        
        rays = pd.read_csv('~/Projects/Epac_fast_anom/HQ_data/HQ_phase_lmm_rays.txt',
                           delim_whitespace=True)        
        for i, ray in rays.iterrows():
            fig.plot(x=[ray.Lon1, ray.Lon2], y=[ray.Lat1, ray.Lat2], pen='2p,red')
    else:
        print('Not adding paths')
    fig.plot(x=scs.LOWMM_LON, y=scs.LOWMM_LAT, style='c0.2c', color='darkred', pen='1p,black')
    fig.plot(x=sks.LOWMM_LON, y=sks.LOWMM_LAT, style='c0.2c', color='green', pen='1p,black')
    fig.plot(x=skks.LOWMM_LON, y=skks.LOWMM_LAT, style='c0.2c', color='orange', pen='1p,black')

    fig.plot(x=hqdata.STLO,y=hqdata.STLA, style='i0.25c', color='red', pen='1p,black')
    #fig.text(x=hqdata.STLO + 0.5, y=hqdata.STLA + 0.5, text=hqdata.STAT)
    
    # fig.plot(x=-140, y=46, style='x0.5c', color='black', pen='2p,black')
    if fname:
        fig.savefig(f'{FIG_DIR}/{fname}.eps',dpi=400, crop=True, show=True)
    else:
        fig.show(method='external')
  
def map_data_paths(data, region=[119, 267, -12, 58], show_pp=False,
                   add_tomo=False, add_roi=True, fname=None):
    ''' 
    Maps data used for inversion. Data that fits our 'high quality' threshold will be drawn on top and 
    the paths will be highlighted 
    '''
    hqdata = data[(data.SNR >= 10.) & ((data.Q >= 0.7) | (data.Q <= -0.9))]
    lqdata = data[~data.index.isin(hqdata.index)]

    fig = pygmt.Figure()
    fig.basemap(region=region,
                projection='B-140/45/40/50/10c', frame=['a5','wesn'])
    # Zoomed out projection is 'B195/25/10/50/10c'
    # Zoomin proejction is B-140/45/40/50/10c
    if add_tomo:
        add_s40rts(fig)
    
    
    fig.coast(shorelines='1/0.5p,black', resolution='l')
    for i, row in lqdata.iterrows():
        fig.plot(x=[row.EVLO, row.STLO], y=[row.EVLA, row.STLA],pen='0.75p,gray30')
    
    for i, row in hqdata.iterrows():
        fig.plot(x=[row.EVLO, row.STLO], y=[row.EVLA, row.STLA],pen='0.75p,darkgoldenrod1')
    
    fig.plot(x=data.EVLO, y=data.EVLA, style='a0.25c', color='yellow', pen='0.5p,black')
    fig.plot(x=data.STLO, y=data.STLA, style='i0.25c', color='red', pen='0.5p,black')

    
    if show_pp:
        fig.plot(x=lqdata.LOWMM_LON, y=lqdata.LOWMM_LAT, style='c0.25c', color='gray30', pen='0.5p,black')
        fig.plot(x=hqdata.LOWMM_LON, y=hqdata.LOWMM_LAT, style='c0.25c', color='darkgoldenrod1', pen='0.5p,black')

    if add_roi:
        fig.plot(x=[-150, -150, -130, -130, -150], y=[40, 50, 50, 40, 40], pen="1p,black")

    if fname:
        fig.savefig(f'{FIG_DIR}/{fname}.png',dpi=720, crop=True, show=True, transparent=True)
    else:
        fig.show(method='external')

def fast_anom_phases(data):
    '''
    Plot phases on fast anom
    '''    
    fig = pygmt.Figure()
    fig.basemap(region=[119, 267, -12, 58],
                projection='B195/25/10/50/15c', frame=['a20','WESn'])
    add_s40rts(fig)

def add_s40rts(fig):
    ''' adds s40rts tomography and cpt to a PyGMT fig'''
    cpt = f'{FIG_DIR}/S40RTS/S40RTS.cpt'
    s40rts = f'{FIG_DIR}/S40RTS/S40RTS_2800km.grd'
    fig.grdimage(grid=s40rts, cmap=cpt)
    fig.colorbar(cmap=cpt,
            position='jBC+o0c/-1.5c/+w10c/0.4c+h+e+m',
            frame=['a1f0.5g0.25','x+l"dVs (%)"'])
    return fig
 
def map_s40rts(region=[0,360,-80,80]):
    '''
    Maps the tomography model S40RTS

    Args:
        region (array) - bounding co-ordinates of region to draw
    Using pyGMT in place of Cartopy
    '''
    fig = pygmt.Figure()
    fig.basemap(region=region, projection='M12c', frame=['a10', 'WESn'])  
    cpt = f'{FIG_DIR}/S40RTS/S40RTS.cpt'
    s40rts = f'{FIG_DIR}/S40RTS/S40RTS_2800km.grd'
    fig.grdimage(grid=s40rts, cmap=cpt)
    fig.coast(shorelines='0.5p,black', resolution='i')
    fig.colorbar(cmap=cpt,
        position='jBC+o0c/-1.5c/+w10c/0.4c+h+e+m',
        frame=['a1f0.5g0.25','x+l"dVs (%)"'])
    fig.savefig(f'{FIG_DIR}/E_pac_blank_S40RTS.eps',crop=True, show=True)
    fig.show(method='external')

def map_sks_corrections(data):
    '''
    Make a map showing the stations used in inversion along with Jacks stations average SKS results
    used as the reciever side correction
    '''
    data = data[(data.SNR >= 10.) & ((data.Q >= 0.7) | (data.Q <= -0.9))]
    corr = pd.read_csv('/Users/ja17375/SWSTomo/Inversions/Station_Corrections_from_Stacks.txt',
                       delim_whitespace=True)
    c2plot = corr[corr.STAT.isin(data.STAT)]
    vec1 = [c2plot.Gamma.values, c2plot.Strength.values*50]
    vec2 = [c2plot.Gamma.values + 180, c2plot.Strength.values*50]
    fig = pygmt.Figure()
    fig.basemap(region=[-130, -100, 25, 50], projection='M15c', frame='a5')
    fig.coast(shorelines='0.5p,black', resolution='l')
    fig.plot(x=data.STLO, y=data.STLA, direction=vec1, style='V0.75', pen='1p,black')
    fig.plot(x=data.STLO, y=data.STLA, direction=vec2, style='V0.75', pen='1p,black')
    
    fig.plot(x=data.STLO, y=data.STLA, style='i0.25c', color='red', pen='1p,black')
    for i,stat in data.iterrows():
        if stat.STAT == 'FACU':
            fig.text(x=stat.STLO, y=stat.STLA-0.5, text=stat.STAT)
        else:
            fig.text(x=stat.STLO, y=stat.STLA +0.5, text=stat.STAT)
    fig.plot(x=-125, y=27, direction=[[90], [0.5]], style='V0.75', pen='2p,black')
    fig.plot(x=-125, y=27, direction=[[270], [0.5]], style='V0.75', pen='2p,black')
    fig.text(x=-125, y=27.6, text='1 s')
    fig.show(method='external')

def global_phase_map(data, fname=None):
    
    fig = pygmt.Figure()
    fig.basemap(region="g",
                projection='ks12c', frame=['a20','WESn'])
        
   # fig.coast(shorelines='1/0.5p,black', resolution='l')
    
    fig.plot(x=data.EVLO, y=data.EVLA, style='a0.25c', color='yellow', pen='0.5p,black')
    fig.plot(x=data.STLO, y=data.STLA, style='i0.25c', color='red', pen='0.5p,black')
   
    scs = data[data.PHASE == 'ScS']
    sks = data[data.PHASE == 'SKS']
    skks = data[data.PHASE == 'SKKS']

    fig.plot(x=scs.LOWMM_LON, y=scs.LOWMM_LAT, style='c0.2c', color='darkred', pen='1p,black')
    fig.plot(x=sks.LOWMM_LON, y=sks.LOWMM_LAT, style='c0.2c', color='green', pen='1p,black')
    fig.plot(x=skks.LOWMM_LON, y=skks.LOWMM_LAT, style='c0.2c', color='orange', pen='1p,black')
    if fname:
        fig.savefig(f'{FIG_DIR}/{fname}.png')
    else:
        fig.show()

if __name__ == '__main__':
    phasefile = 'All_phases.sdb'
    date_time_conv = {'TIME': lambda x: str(x),
                      'DATE': lambda x: str(x)}
    filepath = f'/Users/ja17375/Projects/Splitting_PhD/{phasefile}'
    data = pd.read_csv(filepath, converters=date_time_conv, delim_whitespace=True)
    clat = 46
    clon = -140.
    crit = 4.
    idx = []
    for i, row in data.iterrows():
        if vincenty_dist(clat, clon, row.LOWMM_LAT, row.LOWMM_LON)[0] <= crit:
            idx.append(i)
    data2plot = data.iloc[idx]
                          
    # map_data_paths(data2plot, show_pp=True,
    #                add_tomo=True, fname='data_paths_summary')
    map_data_paths(data2plot, show_pp=True, region=[-150,-130, 40, 50],
                   add_tomo=True, fname='data_paths_zoomin', add_roi=False)
    # map_data_tomo(data2plot,[-150, -100, 25, 50], draw_paths=True, fname='paths_in_inversion_stat_labelled')
    #map_data_tomo(data2plot,[-160, -100, 25, 55], draw_paths=True, fname='colour_path_test')
    #global_phase_map(data)
    #map_s40rts([-170, -80, 10, 60])