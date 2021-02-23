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

def _plot_1d_ppd(ax,V,P,pthresh=0.95,orientation='vertical'):
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
    if orientation == 'vertical':
        ax.fill(v_poly,p_poly,linestyle='-',color=(0.7,0.7,0.7)) # Polygon for full region
        ax.fill(v_polyt,p_polyt,linestyle='-',color=(0.7,0.7,1.0))
        ax.fill(v_l,p_l,color=(0.85,0.85,1.0))
        # Draw a line atop the shaded hist.
        ax.plot(v_line,p_line,'k')
        ax.set_ylim([0,np.around(1.2*P[imax],decimals=3)])
    elif orientation == 'horizontal':
        ax.fill(p_poly,v_poly,linestyle='-',color=(0.7,0.7,0.7)) # Polygon for full region
        ax.fill(p_polyt,v_polyt,linestyle='-',color=(0.7,0.7,1.0))
        ax.fill(p_l,v_l,color=(0.85,0.85,1.0))
        # Draw a line atop the shaded hist.
        ax.plot(p_line,v_line,'k')
        ax.set_xlim([0,np.around(1.2*P[imax],decimals=3)])
    else:
        print('Orientation {} not recognised'.format(orientation))
    return ax

def plot_2d_mppd(P,XY,dp_x,dp_y,save=False,f_uid=None):
    '''
    plot the 2-D mppd on
    i - [str]. A 3 digits code ('001') that corresponds to the MPPD that we want to plot
    '''
    # Remove line above when devel is done
    x = XY[0] # array of X-parameters
    y = XY[1] # array of Y-parameters
    # Get Header of .p file for domain labels
    (X,Y) = np.meshgrid(x,y)
    
    fig = plt.figure(figsize = (10,10))
    axs = gs(ncols = 5,nrows = 4, hspace =0.3, wspace=0.3)
    ax_y = fig.add_subplot(axs[:-1,0])
    ax_x = fig.add_subplot(axs[-1,1:-1])
    ax_main = fig.add_subplot(axs[:-1,1:-1],sharey=ax_y,sharex=ax_x)
    ax_c = fig.add_subplot(axs[:-1,-1])

    # Calc 1D Marginals
    P_x = P.sum(axis=0)
    P_y = P.sum(axis=1)
    # Find indicies of most likely Region
    (irow,icol) = np.unravel_index(np.argmax(P,axis=None),P.shape)
    px_cppd = P[icol,:] / P_y[icol]
    py_cppd = P[:,irow] / P_x[irow]
    _plot_1d_ppd(ax_x,x,P_x)
    # ax_x.hist(x,x,weights=px_cppd,orientation='vertical',histtype='stepfilled')

    ax_x.invert_yaxis()
    # ax_y.hist(y,y,weights=py_cppd,orientation='horizontal',histtype='stepfilled')
    _plot_1d_ppd(ax_y,y,P_y,orientation='horizontal')
    # ax_y.set_ylim([-90,90])
    ax_y.invert_xaxis()
    ax_y.set_ylabel('{}'.format(dp_y))
    ax_y.set_xlabel('p({})'.format(dp_y))
    ax_x.set_xlabel('{}'.format(dp_x))
    ax_x.set_ylabel('p({})'.format(dp_x))

    xdiff = abs(x[1] - x[0])
    ydiff = abs(y[1] - y[0])
    print(xdiff)
    xmax = x.max() + (xdiff/2)
    xmin = x.min() - (xdiff/2)
    ymax = y.max() + (ydiff/2)
    ymin = y.min() - (ydiff/2)

    ax_extent= (xmin, xmax, ymin, ymax)
    xlim = [xmin, xmax]
    ylim = [ymin, ymax]
    C = ax_main.imshow(P,extent = ax_extent,aspect='auto',
                   origin='lower',interpolation='none')
    ax_x.set_xticks(np.linspace(xmin, xmax, 9))
    ax_y.set_yticks(np.linspace(ymin, ymax, 13))
    ax_x.set_xlim(xlim)
    ax_y.set_ylim(ylim)
    plt.colorbar(C,cax=ax_c)
        
    #ax_main.plot(x[icol],y[irow],'xb',markersize=15)
    print(r'Most likely (maxima) points is $\gamma = ${}, $s = ${}'.format(y[irow],x[icol]))
    plt.setp(ax_main.get_xticklabels(),visible=False)
    plt.setp(ax_main.get_yticklabels(),visible=False)

    ax_main.set_title('MPPD {}-{} for Domain {}'.format(dp_y.split(':')[1].strip('\n'), 
                                                dp_x.split(':')[1].strip('\n'),
                                                dp_x.split(':')[0]))
    # fig,ax = plt.subplots(1,1)
    # _plot_1d_ppd(ax,x,px_cppd,dp_x)
    if sv == True:
        print(dp_y.split(':')[1].strip('\n'))
        print(dp_x.split(':')[1].strip('\n'))
        plt.savefig('MPPD_{}_{}_{}_{}.png'.format(dp_y.split(':')[1].strip('\n'), 
                                                  dp_x.split(':')[1].strip('\n'),
                                                  dp_x.split(':')[0].strip(' '),
                                                  f_uid),format='png',dpi=400)
    else:
        plt.show()


def write_out_most_likely(idx,g,s):
    '''Finds (again) the most likely solution from the MPPD and writes it to a textfile '''


if __name__ == "__main__":
    # If this script is being run from the command line
    # idx = sys.argv[0] # expected usage is plot_2d_mppd.py 001
    ## Parse argueemnts
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--save",action="store_true",help="saves MPPD plots")
    parser.add_argument("-f","--filename",default="plot",action="store",type=str,help="file name that will be appended to MPPD_xxx")
    parser.add_argument("--fixed",action="store_true",help="Flag for if Domains has been fixed to only have one parameter")
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
    
    if args.fixed:

        P = np.loadtxt('MTS_1D_MPPD.p',comments='%')
        df1d = pd.read_csv('MTS_1D_MPPD.p',delimiter='%',names=['data','dom_par'])
        dom_par = df1d.dom_par.values
        X = np.loadtxt('MTS_1D_MPPD.x',comments='%')
        n = P.size / 50 
        for i in range(0,int(n)):
            d_par_i = dom_par[i]
            if d_par_i.split(':')[-1] == 's':
                # Is this dom-par pair for strength?
                dom = d_par_i.split(':')[0].strip(" ")
                print(dom)
                fig = plt.figure()
                ax = fig.add_subplot(111)
                if n == 1:
                    _plot_1d_ppd(ax, X, P)
                else:
                    _plot_1d_ppd(ax, x, p)
                gamma = input('Enter fixed gamma >')
                ax.set_title('1D PPD for {}. Gamma = {}'.format(dom,gamma))
                ax.set_xlabel('Strength Parameter')
                ax.set_ylabel('p')
                ax.set_xlim([0,0.05])
                if sv == True:
                    plt.savefig('{}_s_1D_PPD.png'.format(dom),format='png',dpi=400)
                else:
                    plt.show()
    else:

        with open('MTS_most_likely.results','w') as writer:
            writer.write('Layer Domain Gamma Strength\n')
            for i,dom in enumerate(idx):
                print(i,dom)
                P = np.loadtxt('MTS_2D_MPPD.{}.p'.format(dom),skiprows=1) # Skip header line
                # Read the headers
                with open('MTS_2D_MPPD.{}.p'.format(dom),'r') as reader:
                    h = reader.readline().strip('%') # read header line and get rid of that pesky % symbol
                    dps = h.split('-')
                    layer = dps[0].split('_')[0]  
                    print(layer)      
                    try:
                        domain = dps[0].split('_')[1].split(':')[0]
                        dp_x = dps[0]
                        dp_y = dps[1]
                    except IndexError:
                        domain = dps[0].split(':')[0]
                    
                        dp_x = dps[0]
                        dp_y = dps[1].strip('\n')

                    XY = np.loadtxt('MTS_2D_MPPD.{}.xy'.format(dom)) # Two-row file containing X and Y parameters
                    # For each MPPD find most likely solution
                    (irow,icol) = np.unravel_index(np.argmax(P,axis=None),P.shape)
                    gamma = XY[0][irow]
                    s = XY[1][icol]
                    writer.write('{} {} {:5.3f} {:5.3f}\n'.format(layer,domain,gamma,s))
                    # if layer == ' RSide':
                    #     print('Skip plotting for Rside corr domains')
                    #     continue
                    # else:
                    plot_2d_mppd(P,XY,dp_x,dp_y,sv,args.filename)
        # Do stuff (plotting mainly)

    # plt.show()
# EOF
