'''
This module hold some useful functions for calculating spheical disntances on a WGS84 ellipsoid. 
The functions vincenty_dist and vincenty_direct have been tested and work correctly

'''

__author__ = "Joseph P.R. Asplet"

import numpy as np
from scipy import sin, cos, tan, arctan, arctan2, arccos, pi, deg2rad,rad2deg,sqrt

def vincenty_dist(lat1,lon1,lat2,lon2,t=1e-12,deg=True):
    '''
    Calculates distance between two points on the WGS84 ellipsoid using the Vincety formula
    Vincenty, T., Direct and Inverse Solutions of Geodesics on the Ellipsoid with application of nested equations, Surrey Review, XXIII (176)
    https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)
    This implementation returns the angular distance between the two points, sigma, or the distance in km

    Args: 
        lat1 (float) - latidude of the first point (deg)
        lon1 (float) - longitude of the first point (deg)
        lat2 (float) - latidude of the second point (deg)
        long2 (float) - longitude of the second point (deg)
        t (float) - Optional. Tolerance level for the inverse solution 
        deg (Boolean) - Optional. Switch what units the angular distance should be reported 
                        in degrees (True) or in kilometers (False)
        
    Returns:
        D (float) - distance between the two input points in degrees (or in km if deg=False)

    Examples:
        >>>vincenty_dist(5,5,10,10)
            7.03998, 44.64333
            
 
    '''
    # Constants for WGS84
    a = 6378137.0 #[m] length of semi-major axis (radius at equator)
    f = 1/298.257223563 # flattening of the ellipsoid
    b = (1 - f) * a  # [m] length of semi-minor axis  (radius at poles)
    tol = t # tolerance for iterations to find L
    #
    phi1, phi2 = deg2rad(lat1), deg2rad(lat2) # convert to radians and name variabes to match Vincety forumula for sanitys sake
    L1, L2 = deg2rad(lon1), deg2rad(lon2)
    L = L2 - L1
    # Calculate the reduced latitudes (latitude on the auxilliary sphere)
    U1 = arctan((1-f)*tan(phi1))
    U2 = arctan((1-f)*tan(phi2))
    lam_old  = L
    i = 0
    while True:
        i +=1
        ss = (cos(U2)*sin(lam_old))**2
        ss += (cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lam_old))**2
        sin_sigma = ss**0.5
        cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lam_old)
        sigma = arctan2(sin_sigma,cos_sigma)

        sin_alpha =  cos(U1)*cos(U2)*sin(lam_old) / sin_sigma
        cos_sq_alpha = 1 - sin_alpha**2
        cos_2_sig_m = cos_sigma - (2*sin(U1)*sin(U2)) / sin_sigma
        C = f*cos_sq_alpha*(4 + f*(4 - 3*cos_sq_alpha)) / 16

        ss = sigma + C*sin_sigma*(cos_2_sig_m + C*cos_sigma*(-1 + 2*cos_2_sig_m**2))
        lam_new = L + (1 - C)*f*sin_alpha*ss

        if abs(lam_old - lam_new) < tol:
            lam = lam_new
#             print(i)
            break
        elif i >= 50:
            # Points are nearly antipodal. As this routine is interested in points much closer together we will fix s to 180 degrees for now
            # Newton's method can be used to improve convergence
#             print('Points near antipodal')
            s = 180
            return s,None
        else:
            lam_old = lam_new # reset lam_old and repeat
    # # Once we have found lambda
    u2 = cos_sq_alpha*((a**2 - b**2)/ b**2)
    k1 = (sqrt(1 + u2) -1) / (sqrt(1 + u2) + 1)
    A = (1 + 0.25*k1**2) / (1 - k1 )
    B = k1*(1 - (3/8)*k1**2)
# A = 1 + (u2/16384)**(4096 + u2*(-768 + u2*(320 - 175*u2)))
#     B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
    ds = cos_2_sig_m + 0.25*B*(cos_sigma*(-1 + 2*cos_2_sig_m**2))
    ds += B*cos_2_sig_m*(-3 + 4*sin_sigma**2)*(-3 + 4*cos_2_sig_m**2) / 6
    delta_sigma = B*sin_sigma*ds
    s = b*A*(sigma - delta_sigma) # Distance in km
    a1 = arctan2(cos(U2)*sin(lam), cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lam))
    a2 = arctan2(cos(U1)*sin(lam), -sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(lam))
    if deg is True:
        D = np.around(rad2deg(sigma),decimals=10)
    elif deg is False:
        D = np.around(s,decimals=4)
    else:
        print("deg is not Boolean")
        
    return D,rad2deg(a1)

def vincenty_direct(lat1,lon1,azi,dist,ang_dist=True):
    '''Implementation of Vincenty's formula for the direct geodesic problem
    
    Args:
        lat1 (float) - the latitude of the starting point
        lon2 (float) - the longitude of the starting point
        azi (float) - azimuth from the start to the destination point
        dist (float) - distance from start to destination. Default is in degrees
        ang_dist (boolean) - switch to tell function if dist is in degrees (True) or km (False)
        
    Returns:
        lat2 (float) - latitude of the destination point
        lon2 (float) - longitude of the destination point 
    '''
    # Constants for WGS84
    a = 6378137.0 #[m] length of semi-major axis (radius at equator)
    f = 1/298.257223563 # flattening of the ellipsoid
    b = (1 - f) * a  # [m] length of semi-minor axis  (radius at poles)
# Start by calculating the latitude on the auxilliary sphere and the angular distance between the point and the equator
    phi1 = deg2rad(lat1)
    L1 = deg2rad(lon1)
    alpha1 = deg2rad(azi)
    U1 = arctan((1-f)*tan(phi1))
    sigma1 = arctan2(tan(U1),cos(alpha1))
    sin_alpha = cos(U1)*sin(alpha1)
    u2 = (1 - sin_alpha**2)*((a**2 - b**2)/b**2 )
    A = 1 + (u2/16384)**(4096 + u2*(-768 + u2*(320 - 175*u2)))
    B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))

    # As we are providing the angluar distance [distdeg], there is no need to iterate to find sigma
    # For completeness I have implented it anyway
    if ang_dist is True:
        sigma = deg2rad(dist)
        sigm = 2*sigma1 + sigma
    else:
        sigma_old = 0
        sigma = dist/(b*A)
        while (abs(sigma - sigma_old) > 1e-5):
            sigma_old = sigma
            sigm = 2*sigma1 + sigma
            ds = cos(sigma)*(-1 + 2*cos(sigm)**2)
            ds -= (B/6) * cos(sigm) * (-3 + 4*sin(sigma)**2) * (-3 + 4*cos(sigm)**2)
            delta_sigma = B * sin(sigma) * (cos(sigm) + 0.25 * B * ds)
            sigma += delta_sigma
    # Now we can solve for phi2, L2
    sq = sqrt(sin_alpha**2 + (sin(U1)*sin(sigma) - cos(U1)*cos(sigma)*cos(alpha1))**2)
    phi2 = arctan2(sin(U1)*cos(sigma) + cos(U1)*sin(sigma)*cos(alpha1),(1 - f)* sq)
    lam = arctan2(sin(sigma)*sin(alpha1),cos(U1)*cos(sigma) - sin(U1)*sin(sigma)*cos(alpha1))
    cos_sq_alpha = (1 - sin_alpha**2)
    C = (f/16) * cos_sq_alpha*(4 + f*(4 - 3*cos_sq_alpha))
    Lb = (sigma + C*sin(sigma)*(cos(sigm) + C*cos(sigma)*(-1 + 2*cos(sigm)**2)))
    L = lam - (1 - C)*f*sin_alpha*Lb
    L2 = L + L1
    alpha2 = arctan2(sin_alpha,cos(U1)*cos(sigma)*cos(alpha1) - sin(U1)*sin(sigma))
    lat2 = np.around(rad2deg(phi2),decimals=4)
    lon2 = np.around(rad2deg(L2),decimals=4)
    return lat2,lon2