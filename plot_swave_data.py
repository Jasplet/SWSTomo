'''
This script contains functions to plot ScS, SnKS data. Essentially taking over from some of my GMT scripts to make base maps etc.
'''

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import colorbar
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def map_input_data(extent=[-180,-70,0,70]):
    '''
    Make a map showing all our input data
    '''
    scs = pd.read_csv('~/SWSTomo/BluePebble/E_pacific/ScS_05_allQ_corr.sdb',delim_whitespace=True)
    snks = pd.read_csv('~/SWSTomo/BluePebble/E_pacific/SnkS_05_goodQ.pairs',delim_whitespace=True)
    
    lutz = pd.read_excel('../Lutz2020_data.xlsx')
    
    fig = plt.figure(figsize=(12,12))
    proj = ccrs.PlateCarree()
    ax = fig.add_subplot(111,projection=proj)
    
    ax.set_extent(extent, crs=proj)      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    crit = ((snks.Q_SKS >= 0.5) | (snks.Q_SKS <= -0.5)) & ((snks.Q_SKKS >=0.5) | (snks.Q_SKKS <= -0.5))
    ax.plot(scs.BNCLON, scs.BNCLAT, marker='o',lw=0, label='ScS data', transform=proj)
    ax.plot(snks[crit].SKS_PP_LON, snks[crit].SKS_PP_LAT, 
            marker='o',lw=0, label='SKS data', transform=proj)
    ax.plot(snks[crit].SKKS_PP_LON, snks[crit].SKKS_PP_LAT,
            marker='o', lw=0, label = 'SKKS data', transform=proj)
    ax.plot(lutz.PP_LON, lutz.PP_LAT, 'ro' , transform=proj)
    ax.legend()
    
    grd = ax.gridlines(draw_labels=True,linewidth=0)
    grd.top_labels = None
    
    