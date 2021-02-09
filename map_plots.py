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
import cartopy.feature as cfeature 
import cartopy.crs as ccrs
from plot_models import draw_trigonal_doms
from pathset import find_phases_in_domain
from sphe_trig import vincenty_dist

FIG_DIR = '/Users/ja17375/SWSTomo/Figures'

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
                
def map_data_tomo(stats, region=[110,265,-10,55], draw_paths=False, fname=None):

    # if stations == 'all':
    #     stats = data
    # elif stations == 'data-used':
    #     s = ['DAN', 'RDM', 'C10A', 'HUMO', 'L24A', '116A', 'COR']
    #     stats = data[data.STAT.isin(s)]
    scs = stats[stats.PHASE == 'ScS']
    sks = stats[stats.PHASE == 'SKS']
    skks = stats[stats.PHASE == 'SKKS']
    
    fig = pygmt.Figure()
    fig.basemap(region=region, projection='M9i', frame='a10')

    cpt = f'{FIG_DIR}/S40RTS/S40RTS.cpt'
    s40rts = f'{FIG_DIR}/S40RTS/S40RTS_2800km.grd'
    fig.grdimage(grid=s40rts, cmap=cpt)
    fig.colorbar(cmap=cpt,
                position='jBC+o0c/-1.5c/+w8c/0.5c+h+e+m',
                frame=['a1f0.5g0.25','x+l"dVs (%)"'])
    fig.coast(shorelines='1/0.5p,black', resolution='l')        
    # fig.plot(x=d1160_vlon,y=d1160_vlat,pen='2p,gray30,dashed')
    if draw_paths:
        for i, stat in stats.iterrows():
            fig.plot(x=[stat.EVLO, stat.STLO], y=[stat.EVLA, stat.STLA],pen='1p,black')
    else:
        print('Not adding paths')
    fig.plot(x=scs.LOWMM_LON, y=scs.LOWMM_LAT, style='c0.3c', color='darkred', pen='1p,black')
    fig.plot(x=sks.LOWMM_LON, y=sks.LOWMM_LAT, style='c0.3c', color='green', pen='1p,black')
    fig.plot(x=skks.LOWMM_LON, y=skks.LOWMM_LAT, style='c0.3c', color='orange', pen='1p,black')
    fig.plot(x=stats.EVLO, y=stats.EVLA, style='a0.3c', color='black', pen ='black')
    fig.plot(x=stats.STLO,y=stats.STLA, style='i0.4c', color='red', pen='1p,black')
    
    fig.plot(x=-140, y=46, style='x0.5c', color='black', pen='2p,black')
    if fname:
        fig.savefig(f'{FIG_DIR}/{fname}.png',crop=True, show=True)
    else:
        fig.show(method='external')
    
def map_s40rts(region=[-170,-100,10,60]):
    '''
    Maps the tomography model S40RTS

    Args:
        region (array) - bounding co-ordinates of region to draw
    Using pyGMT in place of Cartopy
    '''
    d1160_vlat = [51.5789, 49.2829, 43.4879, 51.5789]   
    d1160_vlon = [-146.1495 , -134.1014, -144.1097, -146.1495]
    fig = pygmt.Figure()
    fig.basemap(region=region, projection='M9i', frame=True)  
    cpt = f'{FIG_DIR}/S40RTS/S40RTS.cpt'
    s40rts = f'{FIG_DIR}/S40RTS/S40RTS_2800km.grd'
    fig.grdimage(grid=s40rts, cmap=cpt)
    fig.coast(shorelines=True, resolution='i')
    fig.plot(
        x=d1160_vlon,
        y=d1160_vlat,
        pen='2p,gray30,dashed')
    fig.show(method='external')
    
if __name__ == '__main__':
    phasefile = 'E_pacific_all_phases.sdb'
    date_time_conv = {'TIME': lambda x: str(x),
                      'DATE': lambda x: str(x)}
    filepath = f'/Users/ja17375/SWSTomo/Inversions/{phasefile}'
    data = pd.read_csv(filepath, converters=date_time_conv, delim_whitespace=True)
    clat = 46
    clon = -140.
    crit = 4.
    idx = []
    for i, row in data.iterrows():
        if vincenty_dist(clat, clon, row.LOWMM_LAT, row.LOWMM_LON)[0] <= crit:
            idx.append(i)
    data2plot = data.iloc[idx]
    data2plot = data2plot[data2plot.SNR >= 10.]
    data2plot = data2plot[(data2plot.Q >= 0.7) | (data2plot.Q <= -0.9)]
    print(data2plot)
                          
    map_data(data2plot, add_tomo=True, draw_paths=False, region=[-170,-100,10,60])
    data2plot.to_csv('/Users/ja17375/SWSTomo/Inversions/HQ_phases_on_fast_anom.sdb', index=False, sep=' ')
    