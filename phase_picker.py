#! /usr/bin/env python
'''
This module contains functions to read raw SAC data and window it for ScS (or any desired phase). Picked phases are output in a format ready to be used by SHEBA to measure shear-wave splitting
'''
import numpy as np
import obspy
import Windower
from obspy.taup import TauPyModel as TauP
import os
from glob import glob

def preproc_sac(st_in,c1 = 0.01, c2 = 0.5):
    '''
    This function reads the raw sac files and pre-processes them by:
        - rotating horizontal components to N,E if needed
        - rtrend, rmean, rtaper
        - applying a bandpass filter 
    
    Args:
        st (Stream) - obspy Stream object containg the data to preprocess
        c1 (float) - lower corner frequency of bandpass filter
        c2 (float) - upper corner frequency of nandpass filter
    '''
    st = st_in.copy()
    st.detrend(type='simple') # de trend data
    st.detrend(type='demean')
#     st.taper(type = 'hann') # hanning window taper
#   filter data
    st.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)

    return st 

def predict_arrivals(sachdr,phases):
    '''
    function to calculate predicted phase arrival times using TauP
    
    Args:
        sachdr (dict) - the SAC header dictionary from the Stream object
        phases (list) - list of phases to use
    
    Returns:
        ttimes (Arrivals object) - obspy object containing phases and travel times
    '''
    model = TauP(model="iasp91")
    if sachdr.evdp >= 1000:
        ttimes = model.get_travel_times(source_depth_in_km=sachdr.evdp/1000, distance_in_degree=sachdr.gcarc,
                               phase_list=phases)
    else:
        ttimes = model.get_travel_times(source_depth_in_km=sachdr.evdp/1000, distance_in_degree=sachdr.gcarc,
                               phase_list=phases)
    
    return ttimes

def make_station_list():
    '''
    function to make the station list for the bash script to rotate traces
    '''
    files = glob('*[00,10]*.SAC')
    with open('station_list','w+') as writer:
        for f in files:
            print(f)
            s1 = f.split('.')[0:3]
            s2 = f.split('.')[4:-1]
            ch = f.split('.')[3]
            name = s1 + s2 + [ch]
            fname = '.'.join(name)
            os.rename(f, fname)
            if ch == 'BHZ':
                stem = s1 + s2
                fstem = '.'.join(stem)
                writer.write('{}\n'.format(fstem))
    
def set_dir():
    '''
    This function sets up the cwd for the rest of Picker to be used. Creates ZNE directory, file list for rotating/ reading files
    '''
    cwd = os.getcwd()
    print("Working in {}".format(cwd))
    os.mkdir('ZNE')
    make_station_list()


if __name__ == "__main__":
    '''
    Script to run picker
    '''
    phases = ['ScS','SS','SKS','SKKS']
    if os.path.isfile('station_list'):
        print('List of files exists')
    else:
        set_dir()
        
    with open('station_list','r+') as reader:
        for line in reader.readlines():
            fstem = line.strip('\n')
            print(fstem)
            st_raw = obspy.read('./ZNE/{}.BH?'.format(fstem))
            st_proc = preproc_sac(st_raw)
            hdr = st_proc[0].stats.sac
            origintime = obspy.UTCDateTime(year = self.BHN[0].stats.sac.nzyear, 
                            julday = self.BHN[0].stats.sac.nzjday,hour=self.BHN[0].stats.sac.nzhour,   
                            minute=self.BHN[0].stats.sac.nzmin,second=self.BHN[0].stats.sac.nzsec,
                            microsecond=self.BHN[0].stats.sac.nzmsec)
            ttimes = predict_arrivals(hdr, phases)
            for tt in ttimes:
                if tt.phase.name == 'ScS'
                    scs_tt = 
    print('Done')