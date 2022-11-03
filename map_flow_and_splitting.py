#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 11:24:01 2021

@author: ja17375
"""
import pygmt
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc

def plot_forte_gmt():

    tx2008 = np.loadtxt('/Users/ja17375/SWSTomo/ForteModels/Flow_Models/TX2008/forteV2_1deg_150km.txt')
    shp = (181, 361)
    dg = 15
    lat = tx2008[:,1].reshape(shp)
    lon = tx2008[:,2].reshape(shp)
    Ur = tx2008[:,3].reshape(shp)
    Utheta = tx2008[:,4].reshape(shp)*-1 # theta is colat so invert 
    Uphi = tx2008[:,5].reshape(shp)
    hzdeg = ((lat % dg == 0) & (lon % dg == 0))
    # Cast Ur (radial velocity) into xarry for pyGMT
    U_grid = xr.DataArray(data=np.flipud(Ur),
                          coords=[('latitude', np.linspace(-90,90,181),
                                       {'units': 'degrees_north'}),
                                       ('longitude', np.linspace(-180,180,361),
                                       {'units': 'degrees_east'})],
                          )
    fig = pygmt.Figure()
    africa_med = [-25,80,-5,60]
    easia = [60,150,10,70]
    epac = [-170, -80, 10, 65]
    proj = "M15c"
    gproj = "Ks12c"
    fig.basemap(region=africa_me, projection=proj, frame="afg",)
    # Flow model TX2008
    # pygmt.makecpt(cmap='roma', series=[-1.5, 1.5], reverse=True)
    # fig.grdimage(grid=U_grid)
    # fig.colorbar(frame=['a0.5', 'x+l"Vertical Velocity (cm/yr)"' ])
    # S40RTS 
    fig.grdimage(grid='/Users/ja17375/DiscrePy/Data/S40RTS/S40RTS_2800km.grd',
                  cmap='/Users/ja17375/DiscrePy/Data/S40RTS/S40RTS.cpt')
    fig.colorbar(frame=['a0.5', 'x+l"dVs (%)"' ], cmap='/Users/ja17375/DiscrePy/Data/S40RTS/S40RTS.cpt')
    fig.coast(shorelines=True)
    # flow_ang = np.rad2deg(np.arctan2(np.ravel(Utheta[hzdeg]), np.ravel(Uphi[hzdeg])))
    # flow_len = np.sqrt(np.ravel(Utheta[hzdeg])**2 + np.ravel(Uphi[hzdeg])**2)
    # flow_data = np.zeros((325, 4))
    # flow_data[:,0] = lon[hzdeg]
    # flow_data[:,1] = lat[hzdeg]
    # flow_data[:,2] = flow_ang
    # flow_data[:,3] = flow_len *0.5
    # fig.plot(data=flow_data, style = 'v0.2c+e', color='black', pen='1p')
    # flow_data[:,2] = flow_data[:,2] + 180
    # fig.plot(data=flow_data, style = 'v0c', color='black', pen='1p')
    fig.plot(x=130, y=20, direction = [[0], [1]], style = 'v0c', color='black', pen='1p')
  
    data = pd.read_csv('~/DiscrePy/Sheba/Results/Combined/Filt_05Hz/Combined_goodQ.pairs', delim_whitespace=True)
    for i, row in data.iterrows():
        fig.plot(x=[row['SKS_PP_LON'], row['SKKS_PP_LON']],
                  y=[row['SKS_PP_LAT'], row['SKKS_PP_LAT']],
                  pen="1p,black")

        if (row['Q_SKS'] >= 0.5):
            #Plot split SKS - black circle
            fig.plot(x=row['SKS_PP_LON'], 
                      y=row['SKS_PP_LAT'],
                      style='c0.15c', color='black', pen='black')
            vec = np.array([[row['SKS_PP_LON'], row['SKS_PP_LAT'], row['FAST_SKS'], row['TLAG_SKS']*0.5],
                            [row['SKS_PP_LON'], row['SKS_PP_LAT'], row['FAST_SKS']+180, row['TLAG_SKS']*0.5]])
            fig.plot(data=vec, style = 'v0c', color='black', pen='0.75p')
            
        elif (row['Q_SKS'] <= -0.5):
            fig.plot(x=row['SKS_PP_LON'], 
                      y=row['SKS_PP_LAT'],
                      style='c0.15c', color='white', pen='black')    
        else:
            print('Bad Q for SKS')
            
        if (row['Q_SKKS'] >= 0.5):
            #Plot split SKKS - black circle
            fig.plot(x=row['SKKS_PP_LON'], 
                      y=row['SKKS_PP_LAT'],
                      style='d0.15c', color='black', pen='black')
            vec = np.array([[row['SKKS_PP_LON'], row['SKKS_PP_LAT'], row['FAST_SKKS'], row['TLAG_SKKS']*0.5],
                    [row['SKKS_PP_LON'], row['SKKS_PP_LAT'], row['FAST_SKKS']+180, row['TLAG_SKKS']*0.5]])
            fig.plot(data=vec, style = 'v0c', color='black', pen='0.75p')
        elif (row['Q_SKKS'] <= -0.5):
            fig.plot(x=row['SKKS_PP_LON'], 
                      y=row['SKKS_PP_LAT'],
                      style='d0.15c', color='white', pen='black')
    fig.savefig('/Users/ja17375/Documents/Thesis-enclosing/Thesis/chapters/chapter02/Figs/Africa_Med_SKS_SKKS_onS40RTS.eps',
                crop=True, show=True)
    # fig.show(method='external')

def plot_flament(dpath='/Users/ja17375/SWSTomo/FlamentModel',extent='epac'):
    nc_vx = nc.Dataset(f'{dpath}/C3-vx-000Ma-2677km.grd')
    nc_vy = nc.Dataset(f'{dpath}/C3-vy-000Ma-2677km.grd')
    nc_vz = nc.Dataset(f'{dpath}/C3-vz-000Ma-2677km.grd')
    vel_conv = 4.9e-4 # converts velocity to cm/year (from N. Flament - see model README.txt)

    Utheta = nc_vx['z'][:] * vel_conv *-1 #theta is colat so invert
    Uphi = nc_vy['z'][:] * vel_conv # longitudl velocity
    Ur = nc_vz['z'][:] * vel_conv # radial velocity
    lon, lat = np.meshgrid(nc_vx['lon'][:], nc_vx['lat'][:])


    dg = 15
    hzdeg = ((lat % dg == 0) & (lon % dg == 0))
    U_grid = xr.DataArray(data=np.flipud(Ur),
                          coords=[('latitude', np.linspace(-90,90,181),
                                       {'units': 'degrees_north'}),
                                       ('longitude', np.linspace(-180,180,361),
                                       {'units': 'degrees_east'})],
                          )
    fig = pygmt.Figure()
    africa_med = [25,70,-5,50]
    fig.basemap(region=africa_med, projection="Ks12c", frame="afg",)
    fig.grdimage(grid=U_grid)
    fig.coast(shorelines=True)
    flow_ang = np.rad2deg(np.arctan2(np.ravel(Utheta[hzdeg]), np.ravel(Uphi[hzdeg])))
    flow_len = np.sqrt(np.ravel(Utheta[hzdeg])**2 + np.ravel(Uphi[hzdeg])**2)
    flow_data = np.zeros((325, 4))
    flow_data[:,0] = lon[hzdeg]
    flow_data[:,1] = lat[hzdeg]
    flow_data[:,2] = flow_ang
    flow_data[:,3] = flow_len *0.1
    fig.plot(data=flow_data, style = 'v0.2c+e', color='black', pen='1p')
    flow_data[:,2] = flow_data[:,2] + 180
    fig.plot(data=flow_data, style = 'v0c', color='black', pen='1p')

    data = pd.read_csv('~/DiscrePy/Sheba/Results/Combined/Filt_05Hz/Combined_goodQ.pairs', delim_whitespace=True)
    for i, row in data.iterrows():
        fig.plot(x=[row['SKS_PP_LON'], row['SKKS_PP_LON']],
                  y=[row['SKS_PP_LAT'], row['SKKS_PP_LAT']],
                  pen="1p,black")

        if (row['Q_SKS'] >= 0.5):
            #Plot split SKS - black circle
            fig.plot(x=row['SKS_PP_LON'], 
                      y=row['SKS_PP_LAT'],
                      style='c0.15c', color='black', pen='black')
            vec = np.array([[row['SKS_PP_LON'], row['SKS_PP_LAT'], row['FAST_SKS'], row['TLAG_SKS']*0.25],
                    [row['SKS_PP_LON'], row['SKS_PP_LAT'], row['FAST_SKS']+180, row['TLAG_SKS']*0.25]])
            fig.plot(data=vec, style = 'v0c', color='black', pen='0.75p')
            
        elif (row['Q_SKS'] <= -0.5):
            fig.plot(x=row['SKS_PP_LON'], 
                      y=row['SKS_PP_LAT'],
                      style='c0.15c', color='white', pen='black')    
        else:
            print('Bad Q for SKS')
            
        if (row['Q_SKKS'] >= 0.5):
            #Plot split SKKS - black circle
            fig.plot(x=row['SKKS_PP_LON'], 
                      y=row['SKKS_PP_LAT'],
                      style='d0.15c', color='black', pen='black')
            vec = np.array([[row['SKKS_PP_LON'], row['SKKS_PP_LAT'], row['FAST_SKKS'], row['TLAG_SKKS']*0.25],
                    [row['SKKS_PP_LON'], row['SKKS_PP_LAT'], row['FAST_SKKS']+180, row['TLAG_SKKS']*0.25]])
            fig.plot(data=vec, style = 'v0c', color='black', pen='0.75p')
        elif (row['Q_SKKS'] <= -0.5):
            fig.plot(x=row['SKKS_PP_LON'], 
                      y=row['SKKS_PP_LAT'],
                      style='d0.15c', color='white', pen='black')
    
    fig.show(method='external')
    
if __name__ == '__main__':
    plot_forte_gmt()