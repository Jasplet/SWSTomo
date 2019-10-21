#! /usr/bin/env python
##############################
#   Program: plot_2d_mppd.py
#
##############################
#   (C) J. Asplet
##############################
#   This program is designed to gather to parse and them plot 2d mppd's (mutli-parameter probability density) output from Matisse.
#   We also plot the 1-D marginals (created by summing across the axes of the 2-D PDF) along each axes
#   The general idea of these plots is taken from James Wookey's matlab scripts which are packaged with Matisse. However here I
#   make some design and functionality (e.g. adding comments) tweaks to make them a little more user friendly (plus they are not in Matlab and are therefore better)

######## Imports ##############

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
###############################

def _plot_1d_ppm(ax,V,P,dompam,pthresh=0.95):
    '''
    Plot a 1-D marginal on the given axes object and returns the axes object (because I am not a bastard)
    ==== Input ====
    ax - [obj].   A matplotlib pyplot axes object (ideally the one you want the 1d PPM to be plotted on)
    V - [array].  Values of the given parameter to plot the 1-D marginal for (e.g. Gamma for Domain A)
    P - [array].  Likelihood for the given set of values
    pthresh - Probability threshold (defualt of 0.95 taking from Wookey's matlab function)
    '''

def _plot_1d_mppd(ax):
    '''
    Copied from JW's matlab scripts, iterates over variables and plots 1d_ppm's
    '''

def plot_2d_mppd(i):
    '''
    plot the 2-D mppd on
    i - [str]. A 3 digits code ('001') that corresponds to the MPPD that we want to plot
    '''
    #Whilst testing
    i = '001'
    # Remove line above when devel is done
    P = np.loadtxt('MTS_2D_MPPD.{}.p'.format(i),skiprows=1) # Skip header line
    XY = np.loadtxt('MTS_2D_MPPD.{}.xy'.format(i)) # Two-row file containing X and Y parameters
    x = XY[0] # array of X-parameters
    y = XY[1] # array of Y-parameters
    # Get Header of .p file for domain labels
    with open(''MTS_2D_MPPD.{}.p'.format(i)','r') as reader:
        h = reader.readine().strip('%') # read header line and get rid of that pesky % symbol
        dps = h.split('-')
        dp_x = doms[0] #Domain & Parameter (should be Domain:Par) in the X dir
        dp_y = doms[1]

    (X,Y) = np.meshgrid()

if __name == "__main__":
    # If this script is being run from the command line

    # Do stuff (plotting mainly)

# EOF
