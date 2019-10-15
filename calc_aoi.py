# convert slowness to AOI based on AK135 velocities
# Usage: slw2aoi(slw) - use depth = 0 (S-wave)
#        slw2aoi(slw,depth) - use depth specified (S-wave)
#        slw2aoi(slw,depth,'p') - use depth specified (P-wave)
# where phase = 'p' or 's' (s is default)
# slowness in s/deg
# Ported from J. Wookeys matlab function slw2aoi.m, here I also calulate the slowness using TauP

import numpy as np
from obspy.taup import TauPyModel
from ak135 import ak135_original
from scipy.interpolate import pchip_interpolate

def get_rayparam(evdp,dist,ph):
    '''
    Use obspy TauP to get the ray parameter for each phase to then get the aoi
    ph [str] - the phase code (SKS or SKKS)
    evdp - the source depth [km]
    dist - the distancefrom source to reciever (for the whole path) [deg]
    '''
    model = TauPyModel('ak135')
    arrival = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=dist, phase_list=[ph])
    rp = arrival[0].ray_param # Rap paramter in s/rad
    return rp

def slw2aoi(depth,evdp,gcarc,phase,wave='s'):
    '''
    Ported (and simplified) version of slw2aoi.m function.
    Phase [str] - should either by 's' or 'p'
    '''
    V=ak135_original()
    slw = get_rayparam(evdp,gcarc,phase)
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
