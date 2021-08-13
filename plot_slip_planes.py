#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 16:18:00 2021

@author: ja17375
"""
import matplotlib.pyplot as plt
import mplstereonet
FIG_DIR = '/Users/ja17375/SWSTomo/Figures/Stereonets'
def plot_stereonet(model):
    
    fig, ax = mplstereonet.subplots(figsize=(3,3))
    ax.plane(model['gamma']+180, model['alpha'],'k')
    ax.line(model['alpha'], model['gamma']-90, 'kx')
    #forte flow model
    ax.line(5.66, 108.65, 'ro')
    # flament flwo model
    ax.line(3.88, 141.766, 'bx')
    ax.set_title(model['title'])   
     # ax.grid()
    fig.savefig('{}/{}_stereonet_w_flow.png'.format(FIG_DIR,model['title']), dpi = 400, transparent=False)
if __name__ == '__main__':
    
    # ellip = {'alpha':25.541, 'gamma':28.197, 'title':'Elliptical TI', 'file':'ellipTI'}
    # br = {'alpha':42.859, 'gamma':22.268, 'title':'Bridgmanite', 'file':'bridgmanite'}
    # ppv001 = {'alpha':22.802, 'gamma':-84.611, 'title':'Post-perovskite [100](001)','file':'ppv001_100'}
    # ppv010 = {'alpha':12.193, 'gamma':93.999, 'title':'Post-perovskite [100](010)','file':'ppv_010_100'}
    
    ellip = {'alpha':74.491, 'gamma':-140.557, 'title':'Elliptical TI', 'file':'ellipTI'}
    br = {'alpha':71.502, 'gamma':-150.188, 'title':'Bridgmanite', 'file':'bridgmanite'}
    ppv001 = {'alpha':60.673, 'gamma':38.264, 'title':'Post-perovskite [100](001)','file':'ppv001_100'}
    ppv010 = {'alpha':34.404, 'gamma':-12.267, 'title':'Post-perovskite [100](010)','file':'ppv_010_100'}
    for model in [ellip, br, ppv001, ppv010]:
        plot_stereonet(model)
    # plot_stereonet(ellip)
    