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

def map_path_counts(cfile,extent=[-160,-70,0,60]):
    '''
    This function makes a map showingn the number of ScS, SKS, SKKS phases that pass throgyuh 
    '''
    domains = np.loadtxt('T3_global.bins',skiprows=1)
    counts = np.loadtxt(cfile)
    doms2plot = domains[np.isin(domains[:,0],counts[:,0])]
    ldom_counts = counts[:,1] + counts[:,2] + counts[:,3]
    fig = plt.figure(figsize=(11,11))
    ax = fig.add_subplot(121,projection=ccrs.PlateCarree())
    ax.set_extent(extent)      
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))   
    for i,dom in enumerate(doms2plot):
        if (dom[2] > extent[0]) & (dom[2] < extent[1]) & (dom[1] > extent[2]) & (dom[1] < extent[3]):
            if (dom[4] ==180) & (dom[2] < 0 ):
                ax.plot([-180,dom[6],dom[8],-180],
                        [dom[3],dom[5],dom[7],dom[3]],'k-',transform=ccrs.PlateCarree()) 
            elif (dom[6] ==180) & (dom[2] < 0 ):  
                ax.plot([dom[4],-180,dom[8],dom[4]],
                        [dom[3],dom[5],dom[7],dom[3]],'k-',transform=ccrs.PlateCarree())
            elif (dom[8] ==180) & (dom[2] < 0 ):  
                ax.plot([dom[4],dom[6],-180,dom[4]],
                        [dom[3],dom[5],dom[7],dom[3]],'k-',transform=ccrs.PlateCarree())
            else:
                ax.plot([dom[4],dom[6],dom[8],dom[4]],
                        [dom[3],dom[5],dom[7],dom[3]],'k-',transform=ccrs.PlateCarree())
        if ldom_counts[i] > 0:
            print(ldom_counts[i])
            ax.scatter(doms2plot[i,2],doms2plot[i,1],c=ldom_counts[i],transform=ccrs.PlateCarree())
            
    ax2 = fig.add_subplot(122,projection=ccrs.PlateCarree())
    ax2.set_extent(extent)      
    ax2.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))   
    udom_counts = counts[:,4] + counts[:,5] + counts[:,6]
    for i,dom in enumerate(doms2plot):
        if (dom[2] > extent[0]) & (dom[2] < extent[1]) & (dom[1] > extent[2]) & (dom[1] < extent[3]):
            if (dom[4] ==180) & (dom[2] < 0 ):
                ax2.plot([-180,dom[6],dom[8],-180],
                        [dom[3],dom[5],dom[7],dom[3]],'k-',transform=ccrs.PlateCarree()) 
            elif (dom[6] ==180) & (dom[2] < 0 ):  
                ax2.plot([dom[4],-180,dom[8],dom[4]],
                        [dom[3],dom[5],dom[7],dom[3]],'k-',transform=ccrs.PlateCarree())
            elif (dom[8] ==180) & (dom[2] < 0 ):  
                ax2.plot([dom[4],dom[6],-180,dom[4]],
                        [dom[3],dom[5],dom[7],dom[3]],'k-',transform=ccrs.PlateCarree())
            else:
                ax2.plot([dom[4],dom[6],dom[8],dom[4]],
                        [dom[3],dom[5],dom[7],dom[3]],'k-',transform=ccrs.PlateCarree())
        if udom_counts[i] > 0:
            print(udom_counts[i])
            ax2.scatter(doms2plot[i,2],doms2plot[i,1],c=udom_counts[i],transform=ccrs.PlateCarree())