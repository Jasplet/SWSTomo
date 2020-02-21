#! /usr/bin/env python
##############################
# sphe_trig.py
##############################
# Author: Joseph P.R. Asplet
##############################
# Holds some useful spherical trig functions
###
import numpy as np


def dist_deg(latA,lonA,latB,lonB):
    '''
    Calculates spherical distance between a pair of points

    Follows Example from the book:
        Plane and Spherical Trigonometry by Kells, Lyman.M; Kern, Willis. F; Bland, James. R (1940). Published by McGraw Hill Book Company
    Available online at archive.org. Identifier: planeandspherica031803mbp

    '''
    dlon = np.abs(lonA - lonB)
    sins = np.sin(np.radians(latA))*np.sin(np.radians(latB))
    coss = np.cos(np.radians(latA))*np.cos(np.radians(latB))*np.cos(np.radians(dlon))
    drad = np.arccos(sins + coss)
    ddeg = np.rad2deg(drad)

    return (np.around(ddeg,decimals=3))
