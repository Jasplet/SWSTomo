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
from matplotlib.gridspec import GridSpec as gs
import sys
from glob import glob
import argparse
###############################

def _plot_1d_ppd(ax,V,P,dompam,pthresh=0.95,orientation='vertical'):
    '''
    Plot a 1-D marginal on the given axes object and returns the axes object (because I am not a bastard)
    ==== Input ====
    ax - [obj].   A matplotlib pyplot axes object (ideally the one you want the 1d PPM to be plotted on)
    V - [array].  Values of the given parameter to plot the 1-D marginal for (e.g. Gamma for Domain A)
    P - [array].  Likelihood for the given set of values
    pthresh - Probability threshold (defualt of 0.95 taking from Wookey's matlab function)
    '''

    nbins = (len(P)-1)
    # Plot the whole PDF
    dx2 = (V[1] - V[0])/2
    v_line = np.array([])
    p_line = np.array([])
    # Find the start / end of each bin and add it with the midpoint (V[i])
    # The points correspond to the top line of the hist
    # Def a lambda function to flatten lists
    flatten = lambda l: [item for sublist in l for item in sublist] # flattens nested lists
    v_line = flatten([ v-dx2, v, v+dx2] for v in V)
    p_line = flatten([pv,pv,pv] for pv in P )

    # Make a copy of the the top lines of the bins so we can turn them into polygons
    v_poly,p_poly = v_line.copy(),p_line.copy()
    v_poly.insert(0,V[0]-dx2)
    p_poly.insert(0,0)
    v_poly.extend([V[-1]+dx2,V[0]])
    p_poly.extend([0,0])
    # Now we have arrays that form closed polygons
    # print(p_poly[0])
    # print(v_poly[0])

    #Calculate and highlight the specified threshold
    psum = 0
    imax = np.argmax(P,axis=None)
    iR, iL = imax,imax # Counters moving right (iR) and left (iL) from the most likely bin
    # print(len(P),np.sum(P))
    while psum<=pthresh:
        # print(psum,iL,iR)
        if ((P[iR] >= P[iL]) or (iL == 0) )and (iR<nbins) :
            psum += P[iR]
            iR += 1 # move one bin to the right (up to nbins)
        elif iL >0:
            psum += P[iL]
            iL += -1 # move one bin the the left (down to 0)
    # Now make polygons for region below threshold
    v_polyt = flatten([ V[i]-dx2, V[i], V[i]+dx2] for i in np.arange(iL,iR))
    p_polyt = flatten([P[i],P[i],P[i]] for i in np.arange(iL,iR))
    v_polyt.insert(0,v_polyt[0])
    p_polyt.insert(0,0)
    v_polyt.extend([v_polyt[-1],v_polyt[0]])
    p_polyt.extend([0,0])
    # Highlight Likliest values
    v_l = [V[imax]-dx2,V[imax],V[imax]+dx2,V[imax]+dx2,V[imax]-dx2]
    p_l = [P[imax],P[imax],P[imax],0,0]



    # Now do the plotting
    if orientation is 'vertical':
        ax.fill(v_poly,p_poly,linestyle='-',color=(0.7,0.7,0.7)) # Polygon for full region
        ax.fill(v_polyt,p_polyt,linestyle='-',color=(0.7,0.7,1.0))
        ax.fill(v_l,p_l,color=(0.85,0.85,1.0))
        # Draw a line atop the shaded hist.
        ax.plot(v_line,p_line,'k')
        ax.set_ylim([0,np.around(1.2*P[imax],decimals=3)])
    elif orientation is 'horizontal':
        ax.fill(p_poly,v_poly,linestyle='-',color=(0.7,0.7,0.7)) # Polygon for full region
        ax.fill(p_polyt,v_polyt,linestyle='-',color=(0.7,0.7,1.0))
        ax.fill(p_l,v_l,color=(0.85,0.85,1.0))
        # Draw a line atop the shaded hist.
        ax.plot(p_line,v_line,'k')
        ax.set_xlim([0,np.around(1.2*P[imax],decimals=3)])
    else:
        print('Orientation {} not recognised'.format(orientation))
    return ax

def plot_2d_mppd(i,save=False,f_uid=None):
    '''
    plot the 2-D mppd on
    i - [str]. A 3 digits code ('001') that corresponds to the MPPD that we want to plot
    '''
    # Remove line above when devel is done
    P = np.loadtxt('MTS_2D_MPPD.{}.p'.format(i),skiprows=1) # Skip header line
    XY = np.loadtxt('MTS_2D_MPPD.{}.xy'.format(i)) # Two-row file containing X and Y parameters
    x = XY[0] # array of X-parameters
    y = XY[1] # array of Y-parameters
    # Get Header of .p file for domain labels
    with open('MTS_2D_MPPD.{}.p'.format(i),'r') as reader:
        h = reader.readline().strip('%') # read header line and get rid of that pesky % symbol
        dps = h.split('-')
        dp_x = dps[0] #Domain & Parameter (should be Domain:Par) in the X dir
        dp_y = dps[1]

    (X,Y) = np.meshgrid(x,y)

    fig = plt.figure(figsize = (8,8))
    axs = gs(4,4, hspace =0.3, wspace=0.3)
    ax_y = fig.add_subplot(axs[:-1,0])
    ax_x = fig.add_subplot(axs[-1,1:])
    ax_main = fig.add_subplot(axs[:-1,1:],sharey=ax_y,sharex=ax_x)

    # Calc 1D Marginals
    P_x = P.sum(axis=0)
    P_y = P.sum(axis=1)
    # Find indicies of most likely Region
    (irow,icol) = np.unravel_index(np.argmax(P,axis=None),P.shape)
    px_cppd = P[icol,:] / P_y[icol]
    py_cppd = P[:,irow] / P_x[irow]
    _plot_1d_ppd(ax_x,x,px_cppd,dp_x)
    # ax_x.hist(x,x,weights=px_cppd,orientation='vertical',histtype='stepfilled')

    ax_x.invert_yaxis()
    # ax_y.hist(y,y,weights=py_cppd,orientation='horizontal',histtype='stepfilled')
    _plot_1d_ppd(ax_y,y,py_cppd,dp_y,orientation='horizontal')
    # ax_y.set_ylim([-90,90])
    ax_y.invert_xaxis()
    ax_y.set_ylabel('{}'.format(dp_y))
    ax_y.set_xlabel('p({})'.format(dp_y))
    ax_x.set_xlabel('{}'.format(dp_x))
    ax_x.set_ylabel('p({})'.format(dp_x))

    # plt.colorbar(C)
    print(x.max())
    if x.max() == 0.0495: # This is the max value of unrestricted domain:
        ax_extent= (0,0.05,-90,90)
        ax_main.imshow(P,extent = ax_extent,aspect='auto',origin='lower')
        ax_x.set_xticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05])
        ax_y.set_yticks([-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
        ax_x.set_xlim([0,0.05])
        ax_y.set_ylim([-90,90])
    else:
        print('Resricted Domain, adjusting scale')
        ax_extent = (0,np.around(x.max(),decimals=2),np.around(y.min(),decimals=0),np.around(y.max(),decimals=0))
        ax_main.imshow(P,extent = ax_extent,aspect='auto',origin='lower')
        print(ax_extent)
        ax_x.set_xlim([np.around(x.min(),decimals=2),np.around(x.max(),decimals=2)])
        ax_y.set_ylim([np.around(y.min(),decimals=0),np.around(y.max(),decimals=0)])

    ax_main.plot(x[icol],y[irow],'xb',markersize=15)
    print(r'Most likely (maxima) points is $\gamma = ${}, $s = ${}'.format(y[irow],x[icol]))
    plt.setp(ax_main.get_xticklabels(),visible=False)
    plt.setp(ax_main.get_yticklabels(),visible=False)

    plt.title('MPPD for Domain {}'.format(dp_x.split(':')[0]))
    # fig,ax = plt.subplots(1,1)
    # _plot_1d_ppd(ax,x,px_cppd,dp_x)
    if sv is True:
        plt.savefig('MPPD_{}_{}.png'.format(dp_x.split(':')[0].strip(' '),f_uid),format='png',dpi=400)

if __name__ == "__main__":
    # If this script is being run from the command line
    # idx = sys.argv[0] # expected usage is plot_2d_mppd.py 001
    ## Parse argueemnts
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--save",action="store_true",help="saves MPPD plots")
    parser.add_argument("-f","--filename",default="plot",action="store",type="str",help="file name that will be appended to MPPD_xxx")
    args = parser.parse_args()

    if args.save:
        print('Plots will be saved')
        # f_uid = input('Enter Unique Identifier for the MPPD plots: ')
        sv = True
    else:
        print('Not Saving')
        # f_uid='placeholder'
        sv = False

    n = glob('MTS_2D_MPPD*.xy')
    idx = [f.split('.')[1] for f in n]
    idx.sort() # sort idx in ascending order (for tidyness sake)
    for i in idx:
        print(i)
        plot_2d_mppd(i,sv,args.filename)
    # Do stuff (plotting mainly)
    plt.show()
# EOF
