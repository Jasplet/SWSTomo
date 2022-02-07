#! /bin/python
# make_synthetic_pathset.py
# Script to make synthetic pathset.XML for a set of synthetics (SWAVXX) 
# Use this for ''inspiration'' as to how the main Pathset setup can be cleaned up and improved
# Author: Joseph Asplet, joseph.asplet@bristol.ac.uk


from xml.etree import ElementTree 
# ElementTree is a standard (as of python 2.5) library which we can use to parse XML.
# N.B ET is NOT secure against "malicous XML" however as we are only intersted in very simple XML this shouldnt be an issue
from xml.dom import minidom
import obspy
import argparse
import numpy as np
from xml.dom import minidom
import os
# Local code (i.e other code i've written)
from wrapper import Wrapper
from sactools import get_sac,get_mts


def set_synth_pathset(idx, synth_dir, azis, incs, dists, pathset_out):
    
    for i in idx:
        n = i + 1 # Becuase SWAVs are numbered from 1 not 0!
        data = obspy.read(f'{synth_dir}/data/SWAV{n:02}.BH?')
        Sheba = Wrapper(data, 'Synth', rundir=f'{synth_dir}/run')
        # Preprocess method filters and trims data
        Sheba.preprocess()
        result = Sheba.measure_splitting(f'SWAV{n:02}', sheba_exec_path='/Users/ja17375/Ext_programs/bin', window=False)
     # Now make XML for this Path
        path = ElementTree.SubElement(pathset,'path')
        pathname = f'Path SWAV {n}'          
        path_uid = ElementTree.SubElement(path,'path_uid')
        path_uid.text = pathname
        # Add Data (from .mts)
        data = get_mts(f'{synth_dir}/run/SWAV{n:02}', syn=True)
        path.append(data)
        stat_uid = ElementTree.SubElement(path,'station_uid')
        stat_uid.text = f'SWAV{n:02}'
        evt_uid = ElementTree.SubElement(path,'event_uid')
        evt_uid.text = 'Synthetic'       
        # Add D`` Operator
        operator = ElementTree.SubElement(path, 'operator')
        dom_uid = ElementTree.SubElement(operator,'domain_uid')
        dom_uid.text = 'Synthetic Test'
        azimuth = ElementTree.SubElement(operator,'azi')
        azimuth.text = str(azis[i])
        # azimuth.text = '0'
        inclination = ElementTree.SubElement(operator,'inc')
        inclination.text = str(incs[i])
        d = ElementTree.SubElement(operator,'dist')
        d.text = str(dists[i])
    
    # write out
    tree = ElementTree.ElementTree(root)
    ElementTree.indent(tree, space="\t", level=0)
    tree.write(pathset_out,
               encoding="utf-8")

if __name__ == '__main__':
  
    parser = argparse.ArgumentParser()
    parser.add_argument("-n","--noise",action="store",required=True, help='Name of noise dir. e.g. a noise level of 0.1 is named Noise01')
    parser.add_argument("-t","--type", action="store",required=True, help='Type of synthetics (real or ideal)')
    parser.add_argument("-b", "--bootstrap", action="store_true", default=False)
    #parser.add_argument("-p","--npaths", action="store",required=True, help='Number of synthetic paths')
    args = parser.parse_args()
    
    if args.type not in ['ideal', 'real']:
        raise ValueError(f'Unexpected synthetic type {args.type}')
    
    # Fix the azimuth, inclination and distances to use
    if args.type == 'ideal':
        azi = np.array([  0.,   40.,   80.,  120.,  160.,  200.,  240.,  280.,  320., 90., 270])
    elif args.type == 'real':
        azi = np.array([96.593, 95.705, 94.830, 88.108, 84.497, 110.886, 114.547, 84.342, 88.017, 112.022, 113.370])
    inc = np.array([57.422, 58.620, 58.102, 60.010, 58.334, 32.828, 34.370, 34.736, 35.814, 34.983, 35.197])
    dists = np.array([296.679, 292.833, 294.468, 288.646, 293.729, 461.181, 442.840, 438.753, 427.231, 872.102, 867.477])
    # dist = 400
    # Set variables
    mtsml = 'http://www1.gly.bris.ac.uk/cetsei/xml/MatisseML/'
    synth_dir = f'/Users/ja17375/Projects/Matisse_Synthetics/ppv1/{args.type}/Noise{args.noise}'
    
    #Start making XML structure
    root = ElementTree.Element('MatisseML')
    root.set("xmlns",mtsml)
    pathset = ElementTree.SubElement(root,'pathset')
    psuid = f'{args.type} BAZ test'
    pathset_uid = ElementTree.SubElement(pathset,'pathset_uid')
    pathset_uid.text = psuid
    # Use sheba to measure splitting for synthetics. We need to do this to get the .mts XML for each path
    if args.bootstrap:
        for i in range(0,500):
            sample_dir = f'{synth_dir}/Bootstrapping/Sample_{i}'
            os.mkdir(sample_dir)
            idx = np.randon.choice(11, 11)
            outfile = f'{sample_dir}/Synthetic_Pathset.xml'
            set_synth_pathset(idx, synth_dir, azi, inc, dists, outfile)
    else:
        idx = np.arange(0, 11)
        outfile = f'/Users/ja17375/Projects/Matisse_Synthetics/ppv1/{args.type}/{args.noise}/Synthetic_{args.type}_Pathset.xml'    
        set_synth_pathset(idx, synth_dir, azi, inc, dists, outfile)
        
    
        
