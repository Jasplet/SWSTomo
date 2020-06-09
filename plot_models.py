''' 
This file contains a set of functions to plot models for Matisse. Including

- Map of ScS, SKS and SKKS data
- Global map of the Trigonal domains 
- Regional view showing phases and domains
- Global view showing Reciever Side corrections
- Regional map showing the number of ScS, SKS and SKKS paths in each domain for a given model (D`` and RSide)

'''

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import colorbar
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

def draw_trigonal_doms(ax, doms2plot,extent=[-180,180,-90,90],fmt='k-', trns=ccrs.PlateCarree()):
    '''
    This function draws the trigonal domains contained in doms2plot on the axis ax
    '''
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
   
            ax.plot([dom[4],dom[6],dom[8],dom[4]],
                    [dom[3],dom[5],dom[7],dom[3]],fmt,transform=trns)

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

def draw_tri_patch(doms2plot,counts, cmin):
    trigs = []
    for i, dom in enumerate(doms2plot):
        if counts[i] > cmin:
            vertices = np.zeros([3,2])
            vertices[:,0] = dom[[4,6,8]]
            vertices[:,1] = dom[[3,5,7]]
            trigs += [Polygon(vertices,closed=True)]

    return trigs
    
def counts_heatmap(cfile,extent=[-180,-70,0,70],cmin=5):
    '''
    This is a function to make a heatmap showing our coverage. Each trigonal domain is coloured based off the number of phases that pass through it
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
    draw_trigonal_doms(ax, doms2plot, extent, fmt='k:')
#   Now draw a subplot for Rside domain counts
    ax2 = fig.add_subplot(212,projection=ccrs.PlateCarree())
    ax2.set_extent(extent, crs=ccrs.PlateCarree())      
    ax2.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))
        
    rdom_counts = counts[:,4] + counts[:,5] + counts[:,6]
    print(rdom_counts[rdom_counts >1])
    trigs2 = draw_tri_patch(doms2plot, rdom_counts, cmin)      
    tps2 = PatchCollection(trigs2, alpha=0.6)
    tps2.set_array(np.array(rdom_counts[rdom_counts > cmin])) # sets colors 
    ax2.add_collection(tps2)
    
    tps2.set_clim([10,ldom_counts.max()])
    fig.colorbar(tps, ax=ax2)
    ax2.set_title(r"Coverage of SKS, SKKS and ScS phases in Upper Mantle (Receiver Side)")
    draw_trigonal_doms(ax2, doms2plot, extent=[-160,-60,10,60], fmt='k:')
    plt.savefig('../Figures/E_pacific_phasecount_heatmap',format=png, dpi=500)
    
def phase_counts_heatmap(cfile='ScS_SnKS.counts', extent=[-170,-60,0,70], cmin=0):
    '''
    This is a function to make seperate maps showing our coverage broken down by phase.
    
    Args:
        cfile [str] - path to the textfile containing the counts data (generated by routine in Pathset.py)
        
        extent [list] - geographic extent of maps [lon_min, lon_max, lat_min, lat_max]. NB current cartopy implementation does not wrap around -180/180. Pick a side!
        
        cmin [int] - the minimum number of counts to draw a polygon for.
        
    Returns:
        lscs_fig [figure] - a 2 panel figure showing ScS and SnKS coverage in D''
        
        rfig [figure] - a 2 pnael figure showing ScS and SnkS coverage ine the Upper Mantle (reciever side)
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
    
    scstrigs = draw_tri_patch(doms2plot, scs_lc, cmin)
    scs = PatchCollection(scstrigs, alpha=0.75)
    scs.set_array(scs_lc[scs_lc > cmin])
    ax1.add_collection(scs)
    draw_trigonal_doms(ax1, doms2plot, fmt='k:')
    cax1,kw =  colorbar.make_axes(ax1,location='bottom',pad=0.05,shrink=0.7)
    cbar_scs = scs_fig.colorbar(scs, cax=cax1, **kw)
    cbar_scs.set_label('No. of ScS phases')
    ax1.set_title(r"ScS coverage in D$''$ domains")
    grd = ax1.gridlines(draw_labels=True,linewidth=0)
    grd.xlabels_top = None
    # Rside Panel
    ax2 = scs_fig.add_subplot(212, projection=ccrs.PlateCarree())
    ax2.set_extent(extent, crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    
    scstrigs2 = draw_tri_patch(doms2plot, scs_rc, cmin)
    scs2 = PatchCollection(scstrigs2, alpha=0.75)
    scs2.set_array(scs_rc[scs_rc > cmin])
    ax2.add_collection(scs2)
    draw_trigonal_doms(ax2, doms2plot, fmt='k:')
    cax2,kw =  colorbar.make_axes(ax2,location='bottom',pad=0.05,shrink=0.7)
    cbar_scs2 = scs_fig.colorbar(scs2, cax=cax2, **kw)
    cbar_scs2.set_label('No. of phases in domain')
    ax2.set_title(r"ScS coverage in upper mantle (receiver side) domains")
    grd2 = ax2.gridlines(draw_labels=True,linewidth=0)
    grd2.xlabels_top = None
    # SnKS Figure

    snks_fig = plt.figure(figsize=(10,16))
    ax3 = snks_fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax3.set_extent(extent, crs=ccrs.PlateCarree())      
    ax3.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    
    snkstrigs = draw_tri_patch(doms2plot, snks_lc, cmin)
    snks = PatchCollection(snkstrigs, alpha=0.75)
    snks.set_array(snks_lc[snks_lc > cmin])
    ax3.add_collection(snks)
    draw_trigonal_doms(ax3, doms2plot, fmt='k:')
    cax3,kw =  colorbar.make_axes(ax3,location='bottom',pad=0.05,shrink=0.7)
    cbar_snks = snks_fig.colorbar(snks, cax=cax3, **kw)
    cbar_snks.set_label('No. of phases in domain')
    ax3.set_title(r"SKS/ SKKS coverage in D$''$ domains")
    grd3 = ax3.gridlines(draw_labels=True,linewidth=0)
    grd3.xlabels_top = None
    # Rside SnKS Panel
    ax4 = snks_fig.add_subplot(212, projection=ccrs.PlateCarree())
    ax4.set_extent(extent, crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.GSHHSFeature(levels=[1],scale='auto'))
    
    snkstrigs2 = draw_tri_patch(doms2plot, snks_rc, cmin)
    snks2 = PatchCollection(snkstrigs2, alpha=0.75)
    snks2.set_array(snks_rc[snks_rc > cmin])
    ax4.add_collection(snks2)
    draw_trigonal_doms(ax4, doms2plot, fmt='k:')
    cax4,kw =  colorbar.make_axes(ax4,location='bottom',pad=0.05,shrink=0.7)
    cbar_snks2 = scs_fig.colorbar(snks2, cax=cax4, **kw)
    cbar_snks2.set_label('No. of phases in domain')
    ax4.set_title(r"SKS/ SKKS coverage in upper mantle (receiver side) domains")
    grd4 = ax4.gridlines(draw_labels=True,linewidth=0)
    grd4.xlabels_top = None
    
    scs_fig.savefig('../Figures/ScS_coverage_maps_CorrTest_model.png',format='png', dpi=600)
    snks_fig.savefig('../Figures/SnKS_coverage_maps_CorrTest_model.png',format='png', dpi=600)