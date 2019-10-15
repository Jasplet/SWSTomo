# convert slowness to AOI based on AK135 velocities
# Usage: slw2aoi(slw) - use depth = 0 (S-wave)
#        slw2aoi(slw,depth) - use depth specified (S-wave)
#        slw2aoi(slw,depth,'p') - use depth specified (P-wave)
# where phase = 'p' or 's' (s is default)
# slowness in s/deg
# Ported from J. Wookeys matlab function slw2aoi.m, here I also calulate the slowness using TauP

import numpy as np
from obspy.taup import TauPyModel

def get_rayparam(ph,evdp,dist):
    '''
    Use obspy TauP to get the ray parameter for each phase to then get the aoi
    ph [str] - the phase code (SKS or SKKS)
    evdp - the source depth [km]
    dist - the distancefrom source to reciever (for the whole path) [deg]
    '''
    model = TauPyModel('ak135')
    arrival = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=dist, phaselist=[ph])
    rp = arrival[0].ray_param # Rap paramter in s/deg
    return rp

def slw2aoi(slw,depth,phase):

  narg = nargin ;
  A=ak135;

  if narg < 1
     fprintf('ERROR: Need at least 1 argument for slw2aoi.\n') ;
  elseif narg == 1
     depth = 0;
     vel = A(:,3);
  elseif narg == 2
     vel = A(:,3);
  elseif narg == 3
     if strcmpi(phase,'p')
     vel = A(:,2);
     else
     vel = A(:,3);
     end
  end
  akdepth = A(:,1) ;
  RE = 6371.0 ;
  v=interp1(akdepth,vel,depth,'pchip') ;

  psperkm = slw / 111.16 ;

  aoi = asind ( psperkm / (1/v)) ;


  prad = slw / (pi/180) ;
  tmp = prad*v / (RE-depth) ;
  aoi = asind (tmp) ;
  return aoi
