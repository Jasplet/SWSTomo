#! /usr/bin/env python
##############################
# sphe_trig.py
##############################
# Author: Joseph P.R. Asplet
##############################
# Holds some useful spherical trig functions
###
import numpy as np
import math

def geocen_lat(lat):
    '''

    '''
    e2 = 6.69437999014E-3
    gc_lat = math.atan((1-e2)* math.tan(lat))

    return gc_lat

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
