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
    V - [array].  Values of the given parameter to plot the 1-D marginal for (e.g. Gamma for Domain A)
    P - [array]. Likelihood for the given set of values
    pthresh - Probability threshold (defualt of 0.95 taking from Wookey's matlab function)
    '''

def _plot_1d_mppd(ax):
    '''
    Plot 
    '''
def _plot_2d_mppd(ax):
