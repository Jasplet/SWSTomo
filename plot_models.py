''' 
This file contains a set of functions to plot models for Matisse. Including

- Map of ScS, SKS and SKKS data
- Global map of the Trigonal domains 
- Regional view showing phases and domains
- Global view showing Reciever Side corrections
- Regional map showing the number of ScS, SKS and SKKS paths in each domain for a given model (D`` and RSide)

'''

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import colorbar
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from pathset import find_phases_in_domain

def draw_trigonal_doms(ax, doms2plot,extent=[-180,180,-90,90],color='black', trns=ccrs.PlateCarree()):
    '''
    This function draws the trigonal domains contained in doms2plot on the axis ax
    '''
    trig = []
    for i,dom in enumerate(doms2plot):
        if (dom[2] > extent[0]) & (dom[2] < extent[1]) & (dom[1] > extent[2]) & (dom[1] < extent[3]):
            if (dom[4] ==180) & (dom[2] < 0) :
                dom[4] = -180
                if (dom[6] == 180):
                    dom[6] = -180
                elif (dom[8] == 180):
                    dom[8] = -180
            elif (dom[6] == 180) & (dom[2] < 0):
                dom[6] = -180
                if (dom[8] == 180):
                    dom[8] = -180
            elif (dom[8] == 180) & (dom[2] < 0):
                dom[8] = -180 
            vertices = np.zeros([3,2])
            vertices[:,0] = dom[[4,6,8]]
            vertices[:,1] = dom[[3,5,7]]
            trig += [Polygon(vertices,closed=True)]
    
    tps = PatchCollection(trig, alpha=0.6)
    tps.set_linestyle('--')
    tps.set_edgecolor('black')
    tps.set_facecolor('white')
    tps.set_linewidth(2.)
    ax.add_collection(tps) 


def map_path_counts(cfile,extent=[-160,-70,0,60]):
    '''
    This function makes a map showing the number of ScS, SKS, SKKS phases that pass through 
    '''
    domains = np.loadtxt('T3_global.bins',skiprows=1)
    counts = np.loadtxt(cfile)
    doms2plot = domains[np.isin(domains[:,0],counts[:,0])]
    ldom_counts = counts[:,1] + counts[:,2] + counts[:,3]
    fig = plt.figure(figsize=(10,22))
    ax = fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))
    ax.set_title(r"ScS, SKS and SKKS path density in D$''$")
   
    draw_trigonal_doms(ax, doms2plot, extent)
    for i,dom in enumerate(doms2plot):
        if ldom_counts[i] > 0:
            ax.scatter(doms2plot[i,2],doms2plot[i,1],c=ldom_counts[i],transform=ccrs.PlateCarree())
            
    ax2 = fig.add_subplot(212,projection=ccrs.PlateCarree())
    ax2.set_extent(extent, crs=ccrs.PlateCarree())      
    ax2.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))   
    udom_counts = counts[:,4] + counts[:,5] + counts[:,6]
    
    draw_trigonal_doms(ax2, doms2plot, extent)
    for i,dom in enumerate(doms2plot):
        if udom_counts[i] > 0:
            ax2.scatter(doms2plot[i,2],doms2plot[i,1],c=udom_counts[i],transform=ccrs.PlateCarree())
            
    plt.show()
    
def contour_map_counts(cfile,extent=[-170,-90,-10,70]):
    '''
    This fuction makes filled contour maps showing the coverage we have in ScS, SKS ,SKKS phases
    '''
       
    domains = np.loadtxt('T3_global.bins',skiprows=1)
    counts = np.loadtxt(cfile)
    doms2plot = domains[np.isin(domains[:,0],counts[:,0])]
    ldom_counts = counts[:,1] + counts[:,2] + counts[:,3]
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))
    draw_trigonal_doms(ax, doms2plot, extent, fmt='k:')
    
#     cntr = ax.tricontourf(doms2plot[:,2], doms2plot[:,1], ldom_counts, levels=10, vmin=1 )
    
#     fig.colorbar(cntr, ax=ax)
        
def plot_swave_model(mfile,title,save=False):
    '''
    This function plots a surface wave model, with the orientation of the anisotropy at each mesh point. THis can be used to plot a single depth slice of a model or a depth averaged model
    
    Args:
        mfile (str) - file containing the surface wave model to plot. Must be formatted with columns of 
                             lon, lat, phi, strength
        title (str) - title to give the plot              
    Returns:
        map view of the model
    '''
    model = np.loadtxt(mfile)
    if model[1,0] - model[0,0] == 1.0 :
        print('First column is Bin No. adjusting')
        lon = model[:,2]
        lat = model[:,1]
        phi = model[:,3]
        strength = model[:,4]
    else:    
        lon = model[:,0]
        lat = model[:,1]
        phi = model[:,2]
        strength = model[:,3]
    if (phi.max() <= pi/2) & (phi.min() >= -pi/2):
        print('Phi is in raidans, convert to degrees')
        phi = rad2deg(phi)
    
    fig = plt.figure(figsize=(11,11))
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    extent=[-140,-70,0,50]
    ax.set_extent(extent)      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))
    ax.quiver(lon,lat,strength,strength,angles=90-phi,
              headaxislength=0,transform=ccrs.PlateCarree(),
              pivot='mid',units='xy',scale=0.5)
#     ax.quiver(ds.lon,ds.lat,ds.dVs_P,ds.dVs_P,angles=270-ds.phi,headaxislength=0,transform=ccrs.PlateCarree())   
    ax.plot(lon,lat,'r.',transform=ccrs.PlateCarree())
    ax.set_title(title)
    grd = ax.gridlines(draw_labels=True)
    grd.top_labels = None
    if save == True:
        modelname = input('Enter name for the plot >>>')
        plt.savefig('../SchafferSurfaceWaveModels/{}'.format(modelname,dpi=400))
    
    plt.show()

def draw_tri_patch(ax, doms2plot, counts, cmin):
    trig = []
    fig = plt.gcf()
    for i, dom in enumerate(doms2plot):
        if counts[i] > cmin:
            vertices = np.zeros([3,2])
            vertices[:,0] = dom[[4,6,8]]
            vertices[:,1] = dom[[3,5,7]]
            trig += [Polygon(vertices,closed=True)]

    tps = PatchCollection(trig, alpha=0.6)
    tps.set_array(np.array(counts[counts > cmin])) # sets colors 
    tps.set_linestyle('-')
    tps.set_edgecolor('black')
    tps.set_linewidth(2.)
    ax.add_collection(tps)        
    #Add colorbar
    cax,kw =  colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.7)
    cbar = fig.colorbar(tps, cax=cax, **kw)
    cbar.set_label('No. of phases in domain')
    
    return tps
    
def counts_heatmap(cfile,extent=[-180,-70,0,70],cmin=5):
    '''
    This is a function to make a heatmap showing our coverage. Each trigonal domain is coloured based off the number of phases that pass through it
    
    Args:
        cfile [str] - path to the textfile containing the counts data (generated by routine in Pathset.py)
        
        extent [list] - geographic extent of maps [lon_min, lon_max, lat_min, lat_max]. NB current cartopy implementation does not wrap around -180/180. Pick a side!
        
        cmin [int] - the minimum number of counts to draw a polygon for.
        
    Returns:
       fig [figure] - a 2 panel figure showing ScS and SnKS coverage in D'' and the Upper Mantle (reciever side)
    '''
    
    domains = np.loadtxt('T3_global.bins',skiprows=1)
    counts = np.loadtxt(cfile)
    doms2plot = domains[np.isin(domains[:,0],counts[:,0])]
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))
    
    ldom_counts = counts[:,1] + counts[:,2] + counts[:,3]
    trigs = draw_tri_patch(doms2plot, ldom_counts, cmin)      
    tps = PatchCollection(trigs, alpha=0.6)
    tps.set_array(np.array(ldom_counts[ldom_counts > cmin]))
    ax.add_collection(tps)
    
    tps.set_clim([10,ldom_counts.max()])
    fig.colorbar(tps, ax=ax)
    ax.set_title(r"Coverage of SKS, SKKS and ScS phases in D$''$")
#   Now draw a subplot for Rside domain counts
    ax2 = fig.add_subplot(212,projection=ccrs.PlateCarree())
    ax2.set_extent(extent, crs=ccrs.PlateCarree())      
    ax2.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))
        
    rdom_counts = counts[:,4] + counts[:,5] + counts[:,6]
    print(rdom_counts[rdom_counts >1])
    tps2 = draw_tri_patch(ax, doms2plot, rdom_counts, cmin)      
    
    tps2.set_clim([10,ldom_counts.max()])
    fig.colorbar(tps, ax=ax2)
    ax2.set_title(r"Coverage of SKS, SKKS and ScS phases in Upper Mantle (Receiver Side)")
    plt.savefig('../Figures/E_pacific_phasecount_heatmap',format=png, dpi=500)

def map_T3_doms(domfile='T3_global.bins', extent=[-170,-60,0,70]):
    '''
    This is a function to draw all domains with their domain ID, for identificiation of domaisn and reference use with other maps
    
    Args:
        domfile [str] - path to the textfile containing the counts data (generated by routine in Pathset.py)
        
        extent [list] - geographic extent of maps [lon_min, lon_max, lat_min, lat_max]. NB current cartopy implementation does not wrap around -180/180. Pick a side!

    Returns:
        fig [figure] - a figure with all trigonal domains in domfile drawn with the domain ID added
    '''
    
    domains = np.loadtxt(domfile,skiprows=1)
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax1.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    draw_trigonal_doms(ax1, domains)
    for dom in domains:
        if (dom[2] > extent[0]) & (dom[2] < extent[1]) & (dom[3] > extent[2]) & (dom[4] < extent[3]):
            ax1.text(dom[2]-1,dom[1],str(dom[0]),transform=ccrs.PlateCarree())
    
    ax1.set_extent(extent, crs=ccrs.PlateCarree())
    
    plt.savefig('domain_map_w_IDs.eps',bbox_inches='tight')

    
def phase_counts_heatmap_byphase(cfile='ScS_SnKS.counts', extent=[-170,-60,0,70], cmin=0, save=False):
    '''
    This is a function to make seperate maps showing our coverage broken down by phase.
    
    Args:
        cfile [str] - path to the textfile containing the counts data (generated by routine in Pathset.py)
        
        extent [list] - geographic extent of maps [lon_min, lon_max, lat_min, lat_max]. NB current cartopy implementation does not wrap around -180/180. Pick a side!
        
        cmin [int] - the minimum number of counts to draw a polygon for.
        
        save [bool] - switch for if you want to save plots or not
        
    Returns:
        scs_fig [figure] - a 2 panel figure showing ScS coverage in D'' and the Upper Mantle (reciever side)
        
        snks_fig [figure] - a 2 panel figure showing  SnkS coverage in D'' and the Upper Mantle (reciever side)
    '''
    
    domains = np.loadtxt('T3_global.bins',skiprows=1)
    counts = np.loadtxt(cfile)
    doms2plot = domains[np.isin(domains[:,0],counts[:,0])]
    
    scs_lc = counts[:,1]
    snks_lc = counts[:,2] + counts[:,3]
    scs_rc = counts[:,4]
    snks_rc = counts[:,5] + counts[:,6]
    
    # ScS Figure
    scs_fig = plt.figure(figsize=(10,16))
    ax1 = scs_fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax1.set_extent(extent, crs=ccrs.PlateCarree())      
    ax1.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    draw_trigonal_doms(ax1, domains)
    scs = draw_tri_patch(ax1, doms2plot, scs_lc, cmin)
    ax1.set_title(r"ScS coverage in D$''$ domains")
    grd = ax1.gridlines(draw_labels=True,linewidth=0)
    grd.xlabels_top = None
    
    # Rside ScS  Panel
    ax2 = scs_fig.add_subplot(212, projection=ccrs.PlateCarree())
    ax2.set_extent(extent, crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    draw_trigonal_doms(ax2, domains)
    scs2 = draw_tri_patch(ax2, doms2plot, scs_rc, cmin)
    ax2.set_title(r"ScS coverage in upper mantle (receiver side) domains")
    grd2 = ax2.gridlines(draw_labels=True,linewidth=0)
    grd2.xlabels_top = None
    
    # SnKS Figure - D'' panel
    snks_fig = plt.figure(figsize=(10,16))
    ax3 = snks_fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax3.set_extent(extent, crs=ccrs.PlateCarree())      
    ax3.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    draw_trigonal_doms(ax3, domains)
    snks = draw_tri_patch(ax3, doms2plot, snks_lc, cmin)
    ax3.set_title(r"SKS/ SKKS coverage in D$''$ domains")
    grd3 = ax3.gridlines(draw_labels=True,linewidth=0)
    grd3.xlabels_top = None
    
    
    # Rside SnKS Panel
    ax4 = snks_fig.add_subplot(212, projection=ccrs.PlateCarree())
    ax4.set_extent(extent, crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    draw_trigonal_doms(ax4, domains)
    
    snks2 = draw_tri_patch(ax4, doms2plot, snks_rc, cmin)
 
    ax4.set_title(r"SKS/ SKKS coverage in upper mantle (receiver side) domains")
    grd4 = ax4.gridlines(draw_labels=True,linewidth=0)
    grd4.xlabels_top = None
    if save:
        scs_fig.savefig('ScS_coverage_maps_1158.png',format='png', dpi=600)
        snks_fig.savefig('SnKS_coverage_maps_1158.png',format='png', dpi=600)

def low_most_mantle_phasecounts(cfile='ScS_SnKS.counts', extent=[-170,-60,0,70], cmin=0,save=False):
    '''
    This is a function to make seperate maps showing our coverage in D`` broken down by phase.
    
    Args:
        cfile [str] - path to the textfile containing the counts data (generated by routine in Pathset.py)
        
        extent [list] - geographic extent of maps [lon_min, lon_max, lat_min, lat_max]. NB current cartopy implementation does not wrap around -180/180. Pick a side!
        
        cmin [int] - the minimum number of counts to draw a polygon for.
                
        save [bool] - switch for if you want to save plots or not
        
    Returns:
        fig [figure] - a 2 panel figure showing coverage in D'' for ScS and SnKS
    '''    
    
    domains = np.loadtxt('T3_global.bins',skiprows=1)
    counts = np.loadtxt(cfile)
    doms2plot = domains[np.isin(domains[:,0],counts[:,0])]
    
    scs_lc = counts[:,1]
    snks_lc = counts[:,2] + counts[:,3]

    # ScS panel
    fig = plt.figure(figsize=(10,16))
    ax1 = fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax1.set_extent(extent, crs=ccrs.PlateCarree())      
    ax1.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    draw_trigonal_doms(ax1, domains)
    scs = draw_tri_patch(ax1, doms2plot, scs_lc, cmin)
    ax1.set_title(r"ScS coverage in D$''$ domains")
    grd = ax1.gridlines(draw_labels=True,linewidth=0)
    grd.top_labels = None
    
    # SnkS panel
    ax2 = fig.add_subplot(212, projection=ccrs.PlateCarree())
    ax2.set_extent(extent, crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    draw_trigonal_doms(ax2, domains)
    scs2 = draw_tri_patch(ax2, doms2plot, snks_lc, cmin)
    ax2.set_title(r"SnKS coverage in D$''$ domains")
    grd2 = ax2.gridlines(draw_labels=True,linewidth=0)
    grd2.top_lables = None
    if save:
        fig.savefig('highquality_lowMM_coverage.png',format='png',dpi=400)
    
def plot_phase_data(file='E_pacific_SNR10_goodQ_allphases.sdb'):
    '''
    function to plot input phase data. currently written to ONLY expect ScS 
    '''
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection=proj)
    
    data = pd.read_csv('~/SWSTomo/BlueCrystal/{}'.format(file),delim_whitespace=True)
    scs = data[data.PHASE == 'ScS']
    sks = data[data.PHASE == 'SKS']
    skks = data[data.PHASE == 'SKKS']

    ax.set_extent([-180,-80,0,90], crs=proj)      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    
    crit = (data.LOWMM_LON > -170) & (data.LOWMM_LAT > 0)
    ax.plot(data.EVLO[crit], data.EVLA[crit], color='black', marker='*', markersize=10,
            label='ScS Events', linestyle='', transform=proj)
#     ax.plot(data.STLO[crit], data.STLA[crit], markerfacecolor='red', markersize=8,
#             markeredgecolor='black', linestyle='', marker='v', transform=proj, label='Stations')
    ax.plot(scs.LOWMM_LON[crit], scs.LOWMM_LAT[crit], color='blue', markeredgecolor='black',
            linestyle='', marker='o', transform=proj, label='ScS', markersize=8)
    ax.plot(sks.LOWMM_LON[crit], sks.LOWMM_LAT[crit], color='red', markeredgecolor='black',
            linestyle='', marker='o', transform=proj, label='SKS', markersize=8)
    ax.plot(skks.LOWMM_LON[crit], skks.LOWMM_LAT[crit], color='orange', markeredgecolor='black',
            linestyle='', marker='o', transform=proj, label='SKKS', markersize=8)
#         ax.plot(scs.ENTLON, scs.ENTLAT, color='green', markeredgecolor='black', linestyle='', marker='>', transform=proj)
#         ax.plot(scs.EXTLON, scs.EXTLON, color='green', markeredgecolor='black', linestyle='', marker='<', transform=proj)
    ax.legend()
    ax.set_title('ScS Data')
    grd = ax.gridlines(draw_labels=True,linewidth=0)
    grd.xlabels_top = None
    
def map_single_domain_phases(phasefile,dom_ID,extent=[-155,-100,30,70]):
    '''
    This function draws a map for a single domain (or potentially multiple domains in the future) along with the phase locations in D`` and the station locations
    '''
    domains = np.loadtxt('T3_global.bins',skiprows=1)
    doms2plot = domains[np.isin(domains[:,0],dom_ID)]
    data = find_phases_in_domain(phasefile,dom_ID)
    scs = data[data.PHASE == 'ScS']
    sks = data[data.PHASE == 'SKS']
    skks = data[data.PHASE == 'SKKS']
    count = len(data)
    # Draw figure
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto')) #coastline data
    draw_trigonal_doms(ax, doms2plot)
    for i, path in data.iterrows():
        ax.plot([path.LOWMM_LON, path.STLO], [path.LOWMM_LAT, path.STLA], 'k-',
               transform=ccrs.Geodetic())
    ax.plot(scs.LOWMM_LON, scs.LOWMM_LAT, markeredgecolor='black',
            linestyle='',color='cornflowerblue', marker='o', transform=proj, label='ScS', markersize=8)
    ax.plot(sks.LOWMM_LON, sks.LOWMM_LAT, markeredgecolor='black',
            linestyle='', color='orange', marker='s', transform=proj, label='SKS', markersize=8)
    ax.plot(skks.LOWMM_LON, skks.LOWMM_LAT, markeredgecolor='black',
            linestyle='', color='mediumseagreen', marker='p', transform=proj, label='SKKS', markersize=8)
    ax.plot(scs.STLO, scs.STLA, linestyle='', marker = 'v', markersize=10, 
           transform=proj, color='cornflowerblue', markeredgecolor='black', label='ScS Stations')
    ax.plot(sks.STLO, sks.STLA, linestyle='', marker = 'v', markersize=10, 
           transform=proj, color='orange', markeredgecolor='black', label='SKS Stations')
    ax.plot(skks.STLO, skks.STLA, linestyle='', marker = 'v', markersize=10, 
           transform=proj, color='mediumseagreen', markeredgecolor='black', label='SKKS Stations')
    
    ax.text(-151,67, '{} ScS phases'.format(len(scs)))
    ax.text(-151,65.75, '{} SKS phases'.format(len(sks)))
    ax.text(-151,64.5, '{} SKKS phases'.format(len(skks)))
    ax.set_title(r"Phases passing through domain Lower {}".format(dom_ID))
    ax.legend()
    grd = ax.gridlines(draw_labels=True,linewidth=0.5)
    grd.top_labels = None
    
    plt.savefig('Phases_in_domain_lower_{}.png'.format(dom_ID),format='png',dpi=400,)