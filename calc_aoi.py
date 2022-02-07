# convert slowness to AOI based on AK135 velocities
# Usage: slw2aoi(slw) - use depth = 0 (S-wave)
#        slw2aoi(slw,depth) - use depth specified (S-wave)
#        slw2aoi(slw,depth,'p') - use depth specified (P-wave)
# where phase = 'p' or 's' (s is default)
# slowness in s/deg
# The function slw2aoi is a ported version of a function slw2aoi.m, written by James Wookey, which I have adpated.
# Here I also calulate the slowness using TauP

import numpy as np
from obspy.taup import TauPyModel
from ak135 import ak135_original
from scipy.interpolate import pchip_interpolate

def get_rayparam(evdp,dist,phase):
    '''
    This function uses TauP (from obspy) to get the ray parameter for each phase to then get the aoi
    
    Args:
        evdp (float) - the event depth [km]
        dist (float) - the distance from source to reciever (for the whole path) [deg]
        phase (str) - the phase code (SKS or SKKS)

    Returns:
        rp (float) - the ray parameter for the given event
        
    Examples:
        >>> get_rayparam(10,100,'SKS')
            281.70596459588847
    '''
    model = TauPyModel('ak135')
    arrival = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=dist, phase_list=[phase])
    rp = arrival[0].ray_param # Rap paramter in s/rad
    return rp

def slw2aoi(depth,slw,wave='s'):
    '''
    This function calculates the angle of incidence of a seismic ray with a given ray      
    parameter (slw)
    
    Args:
        depth (float) - the depth of the ray of interest [km]
        slw (float) - the ray parameter to convert to angle of incidence
        wave (str) - [optional] where the seismic phase is a 'p' or 's' wave. Default is 's'.

    Returns:
        aoi (float) - the angle of incidence of the seismic ray of interest [deg]

    Examples:
        >>> slw2aoi(100,281.7,'s')
        11.650320444422107
    '''
    V=ak135_original()
    RE = 6371.0 # [km] earth radii
    if wave == 'p':
        vel = V[:,1] # p-wave velocity
    else:
        vel = V[:,2] # s-wave velocity

    akdepth = V[:,0]
    # v=interp1(akdepth,vel,depth,'pchip')
    v = pchip_interpolate(akdepth,vel,depth)

    # psperkm = slw / 111.16
    # aoi = np.arcsin( psperkm /  (1/v))
    # prad = slw / (np.pi/180) -- DONT need to do this as it is done in Obspy
    tmp = slw*v / (RE-depth)
    # print(tmp)
    aoi = np.degrees(np.arcsin(tmp))

    return aoi
