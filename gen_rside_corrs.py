'''
This module contains functions needed to create reciever side corrections for our T3 level domains based off one of the Schaeffer surface-wave anisotropy models
'''

from numpy import mean,deg2rad,rad2deg,pi,sin,cos,arctan2
import numpy as np 

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
        print(rad2deg(res))
        res = pi/2 - abs(pi/2 - abs(res))
        print(rad2deg(res))
    elif res > (pi/2):
        print(rad2deg(res))
        res = -pi/2 + abs(res - pi/2)
    
    return res

def depth_stack_point(mod_point,depth_max=330):
    '''
    This function takes a single gird (lat,lon) point in a surface wave model and returns a weighted average of the modelled fast direction across all the depth slices
    
    Args:
        mod_point (obj) - numpy array containing [depth,lat,lon,phi]
        depth_max (int) - the maximum depth to average up to. E.g. the default is to only average over the upper 400 km of the mantle. 

    '''
    n = (mod_point[:,0] <= depth_max).sum() - 1 # have to subtract 1 as indexing start at 0
    if mod_point[n,0] != depth_max:
        weights = np.zeros(n+1)
        phis = np.zeros(n+1)
        strengths = np.zeros(n+1)
        # if the last point is not the same as the specified depth max we will fix it to be so (this is a little iffy)
        weights[n] = (depth_max - mod_point[(n+1),0]) # normalized weighting
    else:
        weights = np.zeros(n)
        phis = np.zeros(n)
        strengths = np.zeros(n)
    
    top = 0 

    for ds in range(0,n):
        phis[ds] = mod_point[ds,3]
        weights[ds] = (mod_point[ds,0] - top) # normalized weighting
        strengths[ds] = mod_point[ds,6]
        top = mod_point[ds,0]

    cm = weighted_circmean(np.deg2rad(phis),np.deg2rad(weights))
    str_m = np.sum(weights * strengths) / np.sum(weights)
    
    return rad2deg(cm),str_m

def depth_stack_model(model_file,depth_max=330):
    '''
    This function takes one of the Schaffer models and creates a depth averaged version up to a set depth
    
    Args:
        model_file (str) - the models to stack
        depth_max (int) - maximum depth to stack up to (default 400km)
    '''
    model = np.loadtxt(model_file)
    dm = model[:,0].min() # get one depth so we can find the size of the lat,lon grid
    ds = model[model[:,0] == dm] # get one depth slice 
    n = ds.shape[0] # length of ds array corresponds to number of points in lat,lon grid
    stacked_model = np.zeros((n,4))

    for i,row in enumerate(ds): # loop over one depth slice as its the easiest way to find all grid nodes
        lon,lat = row[1],row[2] 
        mod_point = model[(model[:,1] == lon) & (model[:,2] == lat)]
        phi, strength = depth_stack_point(mod_point,depth_max)
        stacked_model[i,0] = lon
        stacked_model[i,1] = lat
        stacked_model[i,2] = phi
        stacked_model[i,3] = strength
    
    return stacked_model

def plot_stacked_model(stacked_model):
    '''
    
    '''
    fig = plt.figure()
    