#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 16:18:00 2021

@author: ja17375
"""
import matplotlib.pyplot as plt
import mplstereonet
FIG_DIR = '/Users/ja17375/Projects/Epac_fast_anom/Figures/Stereonets'
def plot_stereonet(model, add_flow=False, add_secondaries=False):
    
    fig, ax = mplstereonet.subplots(figsize=(3,3))
    if add_secondaries:
        for i in range(0,len(model['a2'])):
            ax.plane(model['g2'][i]+180, model['a2'][i], linestyle='--', color='darkgrey')
            ax.line(model['a2'][i], model['g2'][i] - 90, marker='o', color='darkgrey')
            
    ax.plane(model['gamma']+180, model['alpha'],'k')
    ax.line(model['alpha'], model['gamma']-90, 'ko')
    ax.set_title(model['title']) 
    
    if add_flow:
        #add forte flow model
        ax.line(5.66, 108.65, 'ro')
        #add flament flwo model
        ax.line(3.88, 141.766, 'bx')
        fout = f'{FIG_DIR}/{model["title"]}_stereonet_w_flow.png'
    else:
        fout = f'{FIG_DIR}/{model["title"]}_stereonet.png'
     # ax.grid()
    fig.savefig(fout, dpi = 400, transparent=False)
if __name__ == '__main__':
    
    # ellip = {'alpha':25.541, 'gamma':28.197, 'title':'Elliptical TI', 'file':'ellipTI'}
    # br = {'alpha':42.859, 'gamma':22.268, 'title':'Bridgmanite', 'file':'bridgmanite'}
    # ppv001 = {'alpha':22.802, 'gamma':-84.611, 'title':'Post-perovskite [100](001)','file':'ppv001_100'}
    # ppv010 = {'alpha':12.193, 'gamma':93.999, 'title':'Post-perovskite [100](010)','file':'ppv_010_100'}
    
    ellip = {'alpha':81.6, 'gamma':-140.6, 'title':'Elliptical TI', 'file':'ellipTI', 'flag':True,
             'a2':[90], 'g2':[40]}
    br = {'alpha':84.4, 'gamma':-144.1, 'title':'Bridgmanite', 'file':'bridgmanite', 'flag':True,
          'a2':[45], 'g2':[60]}
    ppv001 = {'alpha':27.2, 'gamma':-31.0, 'title':'Post-perovskite [100](001)','file':'ppv001_100',
              'flag':True,'a2':[60, 22, 25], 'g2':[160, 100, 150]}
    ppv010 = {'alpha':28.7, 'gamma':-14.6, 'title':'Post-perovskite [100](010)','file':'ppv_010_100',
              'flag':True, 'a2':[15], 'g2':[120]}
    for model in [ellip, br, ppv001, ppv010]:
        plot_stereonet(model, add_secondaries=model['flag'], add_flow=False)
    # plot_stereonet(ellip)
    