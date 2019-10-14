#! /usr/bin/env python
##############################
#   Program: pathset.py
#
##############################
#   Author: J. Asplet
##############################
#   This program is designed to gather to generate the required Pathset.xml file.
#   Primarily by collecting the .mts files that are output by sheba when measuring SWS, for each phase a seperate <data> tag is required
#   Stating the paths to the requisit SAC files (which are copied to the data directory in the current model path [for now at least])
#

######## Imports ##############
import xml.etree.ElementTree as ET # ElementTree is a standard (as of python 2.5) library which we can use to parse XML.
                                   # N.B ET is NOT secure against "malicous XML" however as we are only intersted in very simple XML this shouldnt be an issue
import pandas as pd
import os
from pathlib import Path
from shutil import copy
import glob
###############################

class Setter:
    """A class to hold the metadata for the run (rdir, station? [for now], outdir etc. ) and fucntions
    to parse/generate the XML needed for the PathSet file
    Inputs: df [obj] - a DataFrame of the relevent .pairs file (Setter will select rows from the relevent station)
            station [str] - the station code for the station[s] we want to include data for (starting with a single station)
            ddir [str] - the data directory. Directory where sheba has output the .mts (and SAC) files to
            odir [str] - the output directory. Path to where we want our output. If none, we use the current working directory
    """
    def __init__(self,df_in,station,ddir,odir=None):

        self.df = df_in[df_in.STAT == station]
        if odir == None:
            self.opath = os.getcwd() # Gets current working directory and stores is as a Path
        self.station = station
        self.ddir = ddir # Data directory (hopefully)

    def get_mts(self,phase,fileID):
        '''Function to get the .mts file for a phase and read in the xml.
           Phase [str] - SKS or SKKS (as my sheba runs are split by phase)
           fileID [str] - file name (.mts is appended here)
        '''
        mts = '{}/{}/{}/{}.mts'.format(self.ddir,self.station,phase,fileID)
        xml = ET.parse(mts) # Parse the xml file (output from sheba as .mts)
        print(xml) # Just print this for now as a test that it works

    def get_sac(self,phase,fileID):
        '''Function to copy data to the local data directory "data" if it is not already there (mainly useful for test/first runs).
           Phase [str] - SKS or SKKS (as my sheba runs are split by phase)
           fileID [str] - file name (.mts is appended here)
        '''
        for comp in ['E','N','Z']:
            f = Path('{}/data/{}.BH{}'.format(self.opath,fileID,comp))
            if f.is_file():
                print('./data/{}.BH{} exists, not copying'.format(fileID,comp))
            else:
                print('File not found, copying from Sheba Run Dir E_pacific if possible')
                file = '{}/{}/{}/{}.BH{}'.format(self.ddir,self.station,phase,fileID,comp)
                dst = '{}/data/{}.BH{}'.format(self.opath,fileID,comp)
                p = copy(file, dst)

    def iter_files(self,phases=['SKS','SKKS']):
        '''
        Function to iterate over the DataFrame and construct the XML we need
        Phases [list] - list of phase codes to iterate over (assuming each row in the DataFame corresponds to all phases)
        '''
        for i, row in self.df.iterrows():
            for ph in phases:
                f = '{}/{}/{}/{}_{}_{}??_{}.mts'.format(self.ddir,row.STAT,ph,row.STAT,row.DATE,row.TIME,ph)
                fileID = glob.glob(f)[0].strip('.mts').split('/')[-1] # Strip out .mts and split by '/', select end to get filestem
                self.get_sac(ph,fileID)
                self.get_mts(ph,fileID)
