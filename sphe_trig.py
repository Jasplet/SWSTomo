#! /usr/bin/env python
##############################
# sphe_trig.py
##############################
# Author: Joseph P.R. Asplet
##############################
# Holds some useful spherical trig functions
# vincenty_dist actually works. The other functions are a mish mash of other attempts and
# ports of J.Wookeys STADIS2 (which also uses the vincenty formula but my port doesnt work... go figure)
#
###
import numpy as np
from scipy import sin, cos, tan, arctan, arctan2, arccos, pi, deg2rad,rad2deg,sqrt

def vincenty_dist(lat1,lon1,lat2,lon2,t=1e-12,deg=True):
    '''
    Calculates distance between two points on the WGS84 ellipsoid using the Vincety formula
    Vincenty, T., Direct and Inverse Solutions of Geodesics on the Ellipsoid with application of nested equations, Surrey Review, XXIII (176)
    https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)

    This implementation returns the angular distance between the two points, sigma, or the distance in km
    '''
    # Constants for WGS84
    a = 6378137.0 #[m] length of semi-major axis (radius at equator)
    f = 1/298.257223563 # flattening of the ellipsoid
    b = (1 - f) * a  # [m] length of semi-minor axis  (radius at poles)
    tol = t # tolerance for iterations to find L
    #
    phi1, phi2 = deg2rad(lat1), deg2rad(lat2) # convert to radians and name variabes to match Vincety forumula for sanitys sake
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
        ss += (cos(U2)*sin(U2) - sin(U1)*cos(U2)*cos(lam_old))**2
        sin_sigma = ss**0.5
        cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lam_old)
        sigma = arctan2(sin_sigma,cos_sigma)

        sin_alpha =  cos(U1)*cos(U2)*sin(lam_old) / sin_sigma
        cos_sq_alpha = 1 - sin_alpha**2
        cos_2_sig_m = cos_sigma - (2*sin(U1)*sin(U2)) / sin_sigma
        C = f*cos_sq_alpha*(4 + f*(4 - 3*cos_sq_alpha)) / 16

        s = sigma + C*sin_sigma*(cos_2_sig_m + C*cos_sigma*(-1 + 2*cos_2_sig_m**2))
        lam_new = L + (1 - C)*f*sin_alpha*s

        if abs(lam_old - lam_new) < tol:
            lam = lam_new

            break
        elif i >= 50:
            # Points are nearly antipodal. As this routine is interested in points much closer together we will fix s to 180 degrees for now
            # Newton's method can be used to improve convergence
            s = 180
        else:
            lam_old = lam_new # reset lam_old and repeat
    # # Once we have found lambda
    u2 = cos_sq_alpha*((a**2 - b**2)/ b**2)
    A = 1 + (u2/16384)**(4096 + u2*(-768 + u2*(320 - 175*u2)))
    B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
    ds = cos_2_sig_m + 0.25*B*(cos_sigma*(-1 + 2*cos_2_sig_m**2))
    ds += B*cos_2_sig_m*(-3 + 4*sin_sigma**2)*(-3 + 4*cos_2_sig_m**2) / 6
    delta_sigma = B*sin_sigma*ds
    s = b*A*(sigma - delta_sigma) # Distance in km

    if deg is True:
        return np.around(rad2deg(sigma),decimals=4)
    elif deg is False:
        return np.around(s,decimals=4)
    else:
        print("deg is not Boolean")

def vincenty_direct(lat1,lon1,azi,dist,ang_dist=True):
    '''Implementation of Vincenty's formula for the direct geodesic problem
       We take an initial point (mlat,mlon), azimuth (azi) and distance (distdeg) and find the end point
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

def dist_deg(qlat,qlon,slat,slon):
    '''
    Calculates spherical distance between a pair of points

    Follows Example from the book:
        Plane and Spherical Trigonometry by Kells, Lyman.M; Kern, Willis. F; Bland, James. R (1940). Published by McGraw Hill Book Company
    Available online at archive.org. Identifier: planeandspherica031803mbp

    '''
    lat1 = geocen_lat(qlat)
    lat2 = geocen_lat(slat)
    dlon = math.radians(qlon - slon)
    dlat = math.radians(lat1 - lat2)
    h = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2 # d = archav(h)
    drad = 2 * math.atan2(math.sqrt(h), math.sqrt(1-h)) # arcsin(X) = arctan(X / sqrt( 1 - X^2))
    ddeg = math.degrees(drad)
    return (np.around(ddeg,decimals=3))


def geocen_colat(arg):
    '''
    input:
    arg    = geographic colatitude (radians)
  output:
    geocen = geocentric colatitude (radians)
    N.B
    fac = (1-f)**2
    conversion gc = pi - atan(cos(colat)/f*sin(colat))) is for colatitude (so dont spaz out when it is different
    from the coversion for latitude you'll find online!)
    '''
    fac = 0.993305621334896
    geocen = np.pi - np.arctan(fac*np.cos(arg)/np.max([1e-30,np.sin(arg)]))

    return geocen

def stadispy(qlat,qlon,slat,slon):
    '''
     Computes the epicentral distance and azimuth from source to receiver.
     Latitudes are converted to geocentric latitudes prior to performing
     the computations (it is assumed that input latitudes are geographic).
     Inputs:   qlat  =  quake latitude (degrees)
               qlon  =  quake longitude (degrees)
               slat  =  station latitude (degrees)
               slon  =  station longitude (degrees)
     Returns:  delta    =  epicentral distance (degrees)
               az      =  azimuth at quake to receiver, from North (degrees)

     [Ported version of STADIS2 (fortran) written by James Wookey ]
     Distance is calculated using an implementation of the Vicenty Formula
    '''

    qcolat = (90 - qlat)
    if (qlon < 0 ):
        qlon=qlon+360

    scolat = 90 - slat
    if (slon < 0 ):
        slon=slon+360

    print(qlat,qlon,slat,slon)
    print(qcolat,qlon,scolat,slon)
    t1deg=scolat
    p1deg=slon
    colat=qcolat
    lon=qlon
    t1=np.radians(t1deg)
    p1=np.radians(p1deg)
    ## For Point A
    t0 = geocen(np.radians(colat))
    p0 = np.radians(lon)

    c0 = np.cos(t0)
    s0 = np.sin(t0)
    ## For point B
    t2 = geocen(t1)
    c1 = np.cos(t2)
    s1 = np.sin(t2)
    ## Calculated distance
    dp = p1-p0
    co = c0*c1 + s0*s1*np.cos(dp)
    print(co)
    si = np.sqrt(1 - co*co)
    delr = np.arctan2(si,co)
    delta = np.degrees(delr)
    #Calculate Azimuth
    print(c1,c0,co,si,s0)
    caz=(c1-(c0*co))/(si*s0)
    if dp != 0:
        dp2=-dp
    else:
        dp2=dp
    saz=-s1*np.sin(dp2)/si
    azr=np.arctan2(saz,caz)
    az = np.degrees(azr)
    if az < 0 :
        az = az + 360

    return delta,az
