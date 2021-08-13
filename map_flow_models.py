#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 11:53:03 2021

@author: ja17375
"""
import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature as cfeature 
import cartopy.crs as ccrs
import netCDF4 as nc
import pandas as pd 

def extract_forte(la, lo):
    '''For an inou point extract the correct Ur, Utheta, Uphi and calculate orientation of flow vector'''    
    tx2008 = np.loadtxt('/Users/ja17375/SWSTomo/ForteModels/Flow_Models/TX2008/forteV2_E_Pac_250km.txt')
    shp = (15, 23)
    lat = tx2008[:,1].reshape(shp) #5 degree step in lat/lon for region [-180, -70, 0, 70]
    lon = tx2008[:,2].reshape(shp)
    Ur = tx2008[:,3].reshape(shp)
    Utheta = tx2008[:,4].reshape(shp)*-1 # theta is colat so invert 
    Uphi = tx2008[:,5].reshape(shp)
    idx_la = int((la - lat.min()) / 5 )# divide by 5 at thats the spacing
    idx_lo = int((lo - lon.min()) /5 )
    
    calc_azi_inc(Uphi[idx_la, idx_lo],
                 Utheta[idx_la, idx_lo],
                 Ur[idx_la, idx_lo])

def extract_flament(la, lo, dpath='/Users/ja17375/SWSTomo/FlamentModel'):
    nc_vx = nc.Dataset(f'{dpath}/C3-vx-000Ma-2677km.grd')
    nc_vy = nc.Dataset(f'{dpath}/C3-vy-000Ma-2677km.grd')
    nc_vz = nc.Dataset(f'{dpath}/C3-vz-000Ma-2677km.grd')
    vel_conv = 4.9e-4 # converts velocity to cm/year (from N. Flament - see model README.txt)

    Utheta = nc_vx['z'][:] * vel_conv *-1 #theta is colat so invert
    Uphi = nc_vy['z'][:] * vel_conv # longitudl velocity
    Ur = nc_vz['z'][:] * vel_conv # radial velocity
    idx_la = int(la +90) # flament models is at 1 deg spacing so indicies correspond to lat/lon
    idx_lo = int(lo)
    calc_azi_inc(Uphi[idx_la, idx_lo],
                 Utheta[idx_la, idx_lo],
                 Ur[idx_la, idx_lo])

def calc_azi_inc(vx, vy, vz):
    '''calculate azimuth and inclination(?) from flow model'''
    print(vx, vy, vz)
    azi = np.rad2deg(np.arctan2(vy, vx))
    incl = np.rad2deg(np.arctan2(vz, np.sqrt(vx**2+vy**2)))
    print(f' {azi:4.2f} {incl:4.2f}')
    return azi, incl

def plot_forte(fig, ax, proj,extent='epac'):
    if extent == 'global':
        tx2008 = np.loadtxt('/Users/ja17375/SWSTomo/ForteModels/Flow_Models/TX2008/forteV2_1deg_150km.txt')
        shp = (181, 361)
        dg = 15
    elif extent == 'epac':
        tx2008 = np.loadtxt('/Users/ja17375/SWSTomo/ForteModels/Flow_Models/TX2008/forteV2_E_Pac_250km.txt')
        shp = (15, 23)
        dg = 10
    lat = tx2008[:,1].reshape(shp) #5 degree step in lat/lon for region [-180, -70, 0, 70]
    lon = tx2008[:,2].reshape(shp)
    Ur = tx2008[:,3].reshape(shp)
    Utheta = tx2008[:,4].reshape(shp)*-1 # theta is colat so invert 
    Uphi = tx2008[:,5].reshape(shp)
    hzdeg = ((lat % dg == 0) & (lon % dg == 0))
    
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    proj = ccrs.PlateCarree()
    if extent == 'global':
        ax.set_global()
    elif extent == 'epac':
        ax.set_extent([-170, -80, 10, 60])
    ax.coastlines()
    ax.gridlines(draw_labels=True, linewidth=0)
    C = ax.contourf(lon, lat, Ur,
                    cmap='PuOr_r', levels=np.linspace(-1.8, 1.8, 37), vmax=1.8, vmin=-1.8,
                    transform=proj) # linspace(-1.4, 1.4, 15) for Epac
    Q = ax.quiver(lon[hzdeg], lat[hzdeg], Uphi[hzdeg], Utheta[hzdeg], pivot='middle', transform=proj, scale=0.2,
              scale_units='x',width=0.003, minlength=1)
    ax.set_title('TX2008 150km above CMB')
    add_sws(ax, proj)
    fig.colorbar(C, ax=ax, orientation='vertical',shrink=0.5, pad=0.08, label='vertical velocity (cm/year)')
    ax.quiverkey(Q, 0.5,-0.25, 5.0, 'velocity vector (5.0 cm/year)', labelpos='N') # 2 cm/yr fror epac
    ax.set_aspect(aspect=1)
    fig.savefig(f'/Users/ja17375/SWSTomo/Figures/TX2008V2_{extent}_splitting.png', dpi=500, transparent=True)


def plot_flament(fig, ax, proj, dpath='/Users/ja17375/SWSTomo/FlamentModel',extent='epac'):
    nc_vx = nc.Dataset(f'{dpath}/C3-vx-000Ma-2677km.grd')
    nc_vy = nc.Dataset(f'{dpath}/C3-vy-000Ma-2677km.grd')
    nc_vz = nc.Dataset(f'{dpath}/C3-vz-000Ma-2677km.grd')
    vel_conv = 4.9e-4 # converts velocity to cm/year (from N. Flament - see model README.txt)

    Utheta = nc_vx['z'][:] * vel_conv *-1 #theta is colat so invert
    Uphi = nc_vy['z'][:] * vel_conv # longitudl velocity
    Ur = nc_vz['z'][:] * vel_conv # radial velocity
    lon, lat = np.meshgrid(nc_vx['lon'][:], nc_vx['lat'][:])

    

    if extent == 'global':
        ax.set_global()
        dg = 15
    elif extent == 'epac':
        ax.set_extent([-170, -80, 10, 60])
        dg = 10
    hzdeg = ((lat % dg == 0) & (lon % dg == 0))

    ax.coastlines()
    ax.gridlines(draw_labels=True, linewidth=0)
    C = ax.contourf(lon, lat, Ur,
                    cmap='PuOr_r', levels=np.linspace(-1.2, 1.2, 25),
                    transform=proj)
    Q = ax.quiver(lon[hzdeg], lat[hzdeg], Uphi[hzdeg], Utheta[hzdeg], pivot='middle', transform=proj,
                  scale=0.2, scale_units='x',width=0.003, minlength=1)
    ax.set_title('Flament (2019) C3 depth 2677 km')
    ax.quiverkey(Q, 0.5,-0.25, 5.0, 'velocity vector (5.0 cm/year)', labelpos='N')
    fig.colorbar(C, ax=ax, orientation='vertical', shrink=0.5, pad=0.08, label='vertical velocity (cm/year)')
    ax.set_aspect(aspect=1)

def add_sws(ax,pj):
    '''Add SKS-SKKS shear-wave splitting  '''
    proj=ccrs.Geodetic()
    data = pd.read_csv('~/DiscrePy/Sheba/Results/Combined/Filt_05Hz/Combined_goodQ.pairs', delim_whitespace=True)
    for i, row in data.iterrows():
        ax.plot([row['SKS_PP_LON'], row['SKKS_PP_LON']], [row['SKS_PP_LAT'], row['SKKS_PP_LAT']],'k-', transform=proj)
        if (row['Q_SKS'] >= 0.5):
            #Plot split SKS - black circle
            ax.plot(row['SKS_PP_LON'], row['SKS_PP_LAT'], 'ko', transform=proj)
        elif (row['Q_SKS'] <= -0.5):
            ax.plot(row['SKS_PP_LON'], row['SKS_PP_LAT'], 'wo', transform=proj, mec='black')
        else:
            print('Bad Q for SKS')
        if (row['Q_SKKS'] >= 0.5):
            #Plot split SKKS - black circle
            ax.plot(row['SKKS_PP_LON'], row['SKKS_PP_LAT'], 'ko', transform=proj)
        elif (row['Q_SKS'] <= -0.5):
            ax.plot(row['SKKS_PP_LON'], row['SKKS_PP_LAT'], 'wo', transform=proj, mec='black')
        
        
if __name__ == '__main__':
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    proj = ccrs.PlateCarree()
    plot_forte(fig, ax, proj, extent='global')
    # extract_forte(45,220)
    # extract_flament(45, 220)
    # fig = plot_flament(fig, ax, proj, extent='global')
    # fig.savefig(f'/Users/ja17375/SWSTomo/Figures/Flament_C3_{extent}.png', dpi=500, transparent=True)
    
    # fig.savefig(f'/Users/ja17375/SWSTomo/Figures/TX2008V2_{extent}.png', dpi=500, transparent=True)
