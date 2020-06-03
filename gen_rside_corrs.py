'''
This module contains functions needed to create reciever side corrections for our T3 level domains based off one of the Schaeffer surface-wave anisotropy models
'''

from numpy import mean,deg2rad,rad2deg,pi,sin,cos,arctan2
import numpy as np 
from scipy import interpolate
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from sphe_trig import vincenty_dist, vincenty_direct
import pandas as pd 
import argparse

def weighted_circmean(angles,weights=[]):
    '''
    A function to take a weighted arithmetic mean of a set of angles that falls within a specified range. Default is -90 to 90 degrees.
    
    Args:
        angles (array): input array. angles should be in degrees (they will be converted to radians inside the function)
    
        weights (array): the weighting to be applied to the angles. default is to set weights to 1 (i.e take an unweighted mean)
     
        high (float or int): high bound of range that angles should fall in
        
        low (float or int): low bound of range that angles should fall in 
        
    Returns:
        weighted mean 
    
    Examples:
        >>>weighted_circmean(deg2rad([80,-80]),[1,1])
            -1.5707963267948966
        
        >>>weighted_circmean([0.1,0.2,0.3],[2,1,2])
            0.20000000000000018
            
        >>>weighted_circmean([0.5,0,-np.pi],[2,1,2])
            0.19558559642933782
        
        >>>weighted_circmean([0.5,0,0],[2,1,2])
            0.19558559642933782
    '''
    high=pi/2
    low=-pi/2
    angles = np.asarray(angles)
    weights = np.asarray(weights)
#   Check there are enough weights for angles 
    if weights.size == 0:
        weights = np.ones(angles.size)
    elif weights.size != angles.size:
        raise ValueError('The number of elements in angles ({}) and weights ({}) do not match'.format(angles.size,weights.size))
    # Translate samples to a range of 0 - 2pi
    samps = (angles - low)*2.*pi / (high - low)

    s_angles = weights * sin(samps) / weights.sum()
    c_angles = weights * cos(samps) / weights.sum()
    s_sum = s_angles.sum()
    c_sum = c_angles.sum()
    tmp = arctan2(s_sum,c_sum)
    # transform result back to a range of -pi/2 - pi/2
    res = tmp*(high - low)/2/pi + low
    if res < -(pi/2):
        res = pi/2 - abs(pi/2 - abs(res))
    elif res > (pi/2):
        res = -pi/2 + abs(res - pi/2)
    
    return res

def depth_stack_point(mod_point,weighting,depth_max=330):
    '''
    This function takes a single gird (lat,lon) point in a surface wave model and returns a weighted average of the modelled fast direction across all the depth slices
    
    Args:
        mod_point (obj) - numpy array containing [depth,lat,lon,phi]
        depth_max (int) - the maximum depth to average up to. E.g. the default is to only average over the upper 400 km of the mantle. 
        weighting (string) - type of weighting to use ['depths','strengths','both']

    Returns
        phi_mean (float) - mean phi value for the given mod point [raidans]
        str_mean (float) - mean strength value for the mod point 
    '''
    n = (mod_point[:,0] <= depth_max).sum() - 1 # have to subtract 1 as indexing start at 0
    if mod_point[n,0] != depth_max:
        weights = np.zeros(n+1)
        phis = np.zeros(n+1)
        strengths = np.zeros(n+1)
        # if the last point is not the same as the specified depth max we will fix it to be so (this is a little iffy)
        if weighting == 'depths':
            weights[n] = (depth_max - mod_point[(n + 1),0]) # normalized weighting
        elif weighting == 'strengths':
            weights[n] = mod_point[(n + 1), 6 ]
        elif weighting == 'both':
            weights[n] = mod_point[(n + 1), 6 ] * (depth_max - mod_point[(n + 1),0])
        else:
            raise ValueError('weighting is not depths, strengths or both')
    else:
        weights = np.zeros(n)
        phis = np.zeros(n)
        strengths = np.zeros(n)
    
    top = 0 
    for ds in range(0,n):
        phis[ds] = mod_point[ds,3]
        if weighting == 'depths':
            weights[ds] = (mod_point[ds ,0] - top)
        elif weighting == 'strengths':
            weights[ds] = mod_point[ds ,6]
        elif weighting == 'both':
            weights[ds] = (mod_point[ds, 0] - top) * mod_point[ds,6]
        else:
            raise ValueError(f'weighting "{weighting}" is not depths, strengths or both')
        
        strengths[ds] = mod_point[ds,6]
        top = mod_point[ds,0]

    phi_mean = weighted_circmean(np.deg2rad(phis),weights)
    str_mean = np.sum(weights * strengths) / np.sum(weights)
    
    return phi_mean ,str_mean

def depth_stack_model(model,weighting='both',depth_max=330):
    '''
    This function takes one of the Schaffer models and creates a depth averaged version up to a set depth
    
    Args:
        model (array-like) - the models to stack
        weighting (string) - type of weighting to use ['depths','strengths','both']
        depth_max (int) - maximum depth to stack up to (default 400km)
        
    Returns:
        stacked_model (ndarray) - numpy array containing the depth stacked model.
    '''
    dm = model[:,0].min() # get one depth so we can find the size of the lat,lon grid
    ds = model[model[:,0] == dm] # get one depth slice 
    n = ds.shape[0] # length of ds array corresponds to number of points in lat,lon grid
    stacked_model = np.zeros((n,4))

    for i,row in enumerate(ds): # loop over one depth slice as its the easiest way to find all grid nodes
        lon,lat = row[1],row[2] 
        mod_point = model[(model[:,1] == lon) & (model[:,2] == lat)]
        phi, strength = depth_stack_point(mod_point, weighting,depth_max)
        stacked_model[i,0] = lon
        stacked_model[i,1] = lat
        stacked_model[i,2] = phi # phi is still in radians here 
        stacked_model[i,3] = strength
    
    return stacked_model

def circ_interp2d(x,y,xnew,ynew,phi,low=-pi/2,high=pi/2):
    '''
    Function to interpolate circular quantities. This is done by interpolating the sin and cos of the unrwapped phi (to 0 - 2pi) and then using atan2 to find the interpolated angle
    
    Args:
        x (array-like) - X-values of original grid
        y (array-like) - Y-values of original grid
        xnew (array-like) - X-values to interpolate for
        ynew (array-like) - Y-values to interpolate for
        phi (array-like) - angles to interpolate. Must be in radians 
        low (float) - lower bound of range of angles (defaults to -pi/2)
        high (float) -upper bound of range of angles (default to pi/2)
        
    Returns:
        phinew (array-like) 
    '''

    phi2int = (phi - low)*2.*pi / (high - low) 
    sinphi = sin(phi2int)
    cosphi = cos(phi2int)
    
#     fsin = interpolate.interp2d(x, y, sinphi, kind = 'cubic')
#     fcos = interpolate.interp2d(x, y, cosphi, kind = 'cubic')
    print(x.size,y.size,sinphi.size,cosphi.size)
    fsin = interpolate.griddata(x, y, sinphi)
    fcos = interpolate.griddata
    +(x, y, cosphi)
    print(xnew.size, ynew.size)
    sinint = fsin(xnew, ynew)
    cosint = fcos(xnew, ynew)
    phiint = arctan2(sinint, cosint)
    
    phinew = phiint * (high - low)/2/pi + low 
    
    return phinew
    
def resample_model(stacked_model,T3_grid='T3_global.bins'):
    '''
    This function takes the depth averaged model and interpolates it onto our T3 domain grid (triangular bin midpoints)
    To interpolate the phi values shift the range to 0 - 180 and use the %180 
    
    Args:
        stacked_model (nd-array) - numpy array containing a depth stacked model
        T3_grid (str) - path to the T3_grid file
    
    Returns:
        resampled_model (nd_array) - numpy array containing the surface wave model interpolated
                                     onto the T3 grid
    '''
    
    T3 = np.loadtxt(T3_grid,skiprows=1)
    stk_x = stacked_model[:,0] - 180 
    stk_y = stacked_model[:,1]
    stk_phi = stacked_model[:,2]
    T3_x = T3[:,1]
    T3_y = T3[:,2]
    
    phi_int = circ_interp2d(stk_x,stk_y,T3_x,T3_y,stk_phi)
    
    return phi_int
    
def plot_model(model,title):
    '''
    This function plots a surface wave model, with the orientation of the anisotropy at each mesh point. THis can be used to plot a single depth slice of a model or a depth averaged model
    
    Args:
        model (array-like) - 2d numpy array containing the model. Must be formatted with columns of 
                             lon, lat, phi, strength
        title (str) - title to give the plot              
    Returns:
        map view of the model
    '''
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
    
    plt.savefig('../SchafferSurfaceWaveModels/SL2016svA_n-k_depth_stacked_depth_weighted',dpi=400)
    
#     plt.show()
    
def comp_grids(grd_in):
    '''
    This is a function to comare the grid of the input surface wavel model and the grid spacing of the T3 grid
    '''
    lon = grd_in[:,0]
    lat = grd_in[:,1]
    T3 = np.loadtxt('T3_global.bins',skiprows=1)
    T3_lon = T3[:,2]
    T3_lat = T3[:,1]
    
    extent=[-140,-70,0,50]
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    ax.set_extent(extent)
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='low'))
    ax.plot(lon,lat,'b.',transform=ccrs.PlateCarree())
    ax.plot(T3_lon,T3_lat,'rx',transform=ccrs.PlateCarree())
    
    plt.show()
    
def find_points(points,qlat,qlon,tdist):
    '''
    This function finds all points within a fixed distance of the query point, assuming a WGS84 ellipsoid. Returning the points and their distance from the query point 
    
    Args:
        points (array-like) - array containing (at least) the lat,lon of the points to search
        qlat (float or int) - latitude of the query point
        qlon (float or int) - longitude of the query point 
        dist (float or int) - the search distance
        
    Returns:
        targets (array - like) - target points, with the distance between each point and the query point appended
    '''

    points = np.asarray(points) # Make sure points is an array, should be lon, lat, phi, str
    if (points[:,1].max() > 90) | (points[:,1].min() < -90):
        raise ValueError("Points array is not right. 2nd column should be latitude!")
    elif (points[:,0].max() > 180) | (points[:,0].min() < -180):
#         raise ValueError("Longitude (points[:,0]) incorrect. Range must be -180 to 180 degrees"  
        tocorr = points[points[:,0] > 180] 
        tocorr[:,0] = tocorr[:,0] - 360
        points[points[:,0] > 180] = tocorr

    print('Query lat {}, query lon {}'.format(qlat,qlon))
    dists = [] # empty list to hold distances 
    tpts = [] # list of indicies 
    for p, point in enumerate(points):
        plon = point[0]
        plat = point[1]
        dist, azi = vincenty_dist(qlat, qlon, plat, plon)
        if dist <= tdist:
            dists.append(dist)
            tpts.append(p)
    
    if len(tpts) == 0:
        raise ValueError('No points found, this should not happen. Check input params')
        
    dists = np.array(dists)
    targets = points[tpts]

    
    return targets, dists
    
def spatial_average_model(stacked_model,T3_grid='/Users/ja17375/SWSTomo/T3_global.bins'):
    '''
    This function take the depth averaged surface wave model and takes a weighted spatial average at each of the T3 mesh points
    '''
    
    T3 = np.loadtxt(T3_grid,skiprows=1)
    T3_swav_corr = np.zeros([T3.shape[0],6])
    T3_swav_corr[:,0:3] = T3[:,0:3] # Populate BIN, MidLat and MidLon
    if (stacked_model[:,0].min() == 0) & (stacked_model[:,0].max() == 360):
        # Need to map lon from 0 - 360 to -180 -180
        stacked_model[:,0] = stacked_model[:,0] - 180 
    
    for i, query_point in enumerate(T3):
        qlat = query_point[1]
        qlon = query_point[2]
        
        pts, weights = find_points(stacked_model, qlat, qlon, tdist=4.)
        npts = pts.shape[0]
        phi = pts[:,2]
        strengths = pts[:,3]
        phi_avg = weighted_circmean(phi, weights)
        str_avg = np.sum(weights * strengths) / np.sum(weights)
        
        T3_swav_corr[i,3] = rad2deg(phi_avg)
        T3_swav_corr[i,4] = str_avg
        T3_swav_corr[i,5] = npts
        
    T3_corr = pd.DataFrame(T3_swav_corr,
                           columns=['BIN','MID_LAT','MID_LON','PHI','STRENGTH','NPTS'])
    
    return T3_corr

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("mod_id",action="store",type=str,help="Filename for model to create Rside corrections from")
#     parser.add_argument("-o","--outfile",action="store",type=str,
#                         default="SL2016svA_resampled_to_T3", 
#                         help="File name for the reciever side correction")
    parser.add_argument("-w","--weightingstyle",action="store",type=str,default="both",
                        choices=['strengths','depths','both', 'all'],
                        help="Style of weighting to use in the depth stakcing. Options are 'depths', 'strengths','both', or 'all'. Default is 'both'")
    parser.add_argument("-d","--depth",action='store',type=int,default=330.,
                        help='Maximum depth (km) to stack the input model to. ')

    args = parser.parse_args()
    mod_id = args.mod_id
    depth = args.depth
    print(mod_id, depth, args.weightingstyle)
    if args.weightingstyle == 'all':
        weighting = ['strengths','depths','both']
    elif args.weightingstyle == 'strengths':
        weighting = ['strengths']
    elif args.weightingstyle == 'depths':
        weighting = ['depths']
    elif args.weightingstyle == 'both':
        weigthing = ['both']
    else:
        raise ValueError('Weighting style not recognised')
        
    model = np.loadtxt(mod_id)
    for w in weighting:
        print('Weighting by {}'.format(w))
        dep_stacked_model = depth_stack_model(model,weighting=w,depth_max=depth)
        T3_model = spatial_average_model(dep_stacked_model)
        T3_model.to_csv('SL2016svAs_T3mesh_w_by_{}.mod'.format(w),sep=' ',index=False)
        
    print('Done')