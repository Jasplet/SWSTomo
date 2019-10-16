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
from xml.etree import ElementTree # ElementTree is a standard (as of python 2.5) library which we can use to parse XML.
                                   # N.B ET is NOT secure against "malicous XML" however as we are only intersted in very simple XML this shouldnt be an issue
from xml.etree.ElementTree import Element,SubElement
from xml.dom import minidom
import pandas as pd
import os
from pathlib import Path
from shutil import copy
import glob
from calc_aoi import slw2aoi
import numpy as np
###############################

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

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
            self.odir = self.opath.split('/')[-1]

        self.station = station
        self.ddir = ddir # Data directory (hopefully)
        self.dom_h = 250 # [km] height of the domains (fixed for UM and D'' for now!)
        self.xmlns = {'mtsML':'http://www1.gly.bris.ac.uk/cetsei/xml/MatisseML/'} # Dict with the Matisse XML namespace (for easy ref)
        ## Read Model.xml (assumed to be in working directory because why wouldnt it?)
        root = ElementTree.parse('Model.xml').getroot()
        self.model = root.find('mtsML:model',self.xmlns)
        model_name = self.model[0].text
        print('Using model named ... {}'.format(model_name))


    def get_mts(self,phase):
        '''Function to get the .mts file for a phase and read in the xml.
           Phase [str] - SKS or SKKS (as my sheba runs are split by phase)
           fileID [str] - file name (.mts is appended here) [ now an attribute of Setter class]
        '''
        mts = '{}/{}/{}/{}.mts'.format(self.ddir,self.station,phase,self.fileID)
        xml = ElementTree.parse(mts) # Parse the xml file (output from sheba as .mts)
        data = xml.getroot() # Gets the root element of the XML. In this case (the .mts) this is the tag <data> which we want to inject into the
                                     # the bigger Pathset XML file
        for file in ['file1','file2','file3']:# Loop over file tags to add in pointer to data dir
            f = data.find(file).text
            print(f)
            f_new = 'data/{}'.format(f)
            data.find(file).text = f_new

        return data

    def get_sac(self,phase):
        '''Function to copy data to the local data directory "data" if it is not already there (mainly useful for test/first runs).
           Phase [str] - SKS or SKKS (as my sheba runs are split by phase)
           fileID [str] - file name (.mts is appended here) (now an attribute)
        '''
        for comp in ['E','N','Z']:
            f = Path('{}/data/{}.BH{}'.format(self.opath,self.fileID,comp))
            if f.is_file():
                print('./data/{}.BH{} exists, not copying'.format(self.fileID,comp))
            else:
                print('File not found, copying from Sheba Run Dir E_pacific if possible')
                file = '{}/{}/{}/{}.BH{}'.format(self.ddir,self.station,phase,self.fileID,comp)
                dst = '{}/data/{}.BH{}'.format(self.opath,self.fileID,comp)
                p = copy(file, dst)

    def domain2operator(self,domain,phase):
        '''Function to read model.xml and extract the required domain
        Input -------
        domain [str] - the domain name (tag <domain_uid>) for the domain we want to find and "cast" as the operator
        '''
        domains = self.model.findall('mtsML:domain',self.xmlns)
        for dom in domains:
            print(dom)
            uid_tmp = dom.find('mtsML:domain_uid',self.xmlns).text
            print(uid_tmp)
            if uid_tmp == domain:
                print('Domain Found')
                uid = uid_tmp
                if uid.split('_')[0] == 'Lower':
                    # Domains starting with Lower are at CMB. So depth == 2890 km
                    depth = 2890. # Approx depth of CMB [km]
                elif uid.split('_')[0] == 'Upper':
                    depth = 250. # Depth of upper domain (keeping domain the same thickness for now)
                else:
                    raise NameError('Domain is incorreclty named. Should be either "Upper" or "Lower".')

                aoi = slw2aoi(depth,self.evdp,self.gcarc,phase) # Calculate ray param and then incidence angle
                dist = self.dom_h / np.cos(np.radians(aoi)) # Calculate distance travelled through domain

        ## Now return what we need to make the opertor (will do this above?)
        operator = ElementTree.Element('operator')
        dom_uid = ElementTree.SubElement(operator,'domain_uid')
        dom_uid.text = uid
        azimuth = ElementTree.SubElement(operator,'azi')
        azimuth.text = str(self.az)
        inclination = ElementTree.SubElement(operator,'inc')
        inclination.text = str(aoi)
        l = ElementTree.SubElement(operator,'dist')
        l.text = str(dist)

        return operator

    def gen_XML(self,phases=['SKS','SKKS']):
        '''
        Function to iterate over the DataFrame and construct the XML we need
        Phases [list] - list of phase codes to iterate over (assuming each row in the DataFame corresponds to all phases)
        '''
        root = ElementTree.Element('MatisseML')
        tree = ElementTree.ElementTree(root)
        root.set("xmlns",self.xmlns['mtsML'])
        pathset = ElementTree.SubElement(root,'pathset')
        ps_uid = 'Path for run in dir {} .'.format(self.odir)
        for i, row in self.df.iterrows():
            # All XML generation must sit within this loop (function calls) so that we make sure that Az, EVDP etc. are right for the current phase
            self.evdp = row.EVDP # read out event depth [km]. Attrib as needed for path length calcs.
            self.az = row.AZI
            self.gcarc = row.DIST
            evla = row.EVLA
            evlo = row.EVLO
            stla = row.STLA
            stlo = row.STLO
            stat = row.STAT
            for ph in phases:
                #Each iteration of this loop is a seperate path (1 per phase and event)
                f = '{}/{}/{}/{}_{}_{}??_{}.mts'.format(self.ddir,stat,ph,stat,row.DATE,row.TIME,ph)
                self.fileID = glob.glob(f)[0].strip('.mts').split('/')[-1] # Strip out .mts and split by '/', select end to get filestem
                self.get_sac(ph)
                # Now make XML for this Path
                path = ElementTree.SubElement(pathset,'path')
                pathname = 'Path {} {}'.format(i,ph)
                path_uid = ElementTree.SubElement(path,'path_uid')
                path_uid.text = pathname
                # Add Data (from .mts)
                data = self.get_mts(ph)
                path.append(data)
                stat_uid = ElementTree.SubElement(path,'station_uid')
                stat_uid.text = stat
                evt_uid = ElementTree.SubElement(path,'event_uid')
                evt_uid.text = '{}_{}'.format(row.DATE,row.TIME)
                # Now we need to select the correct domains and in the right order (order of operators).
                # As the model get more complex these tests will have to get more "clever"
                # Hardcoded for now, need to get a function to read Model.xml (possibly as part of __init__)
                op_UM = self.domain2operator('Upper_01',ph)
                if ph == 'SKS':
                    op_LM = self.domain2operator('Lower_01',ph) # This will need to be a dictionary of domain UIDs
                elif ph == 'SKKS':
                    op_LM = self.domain2operator('Lower_02',ph)
                path.append(op_LM)
                path.append(op_UM)

        return tree
