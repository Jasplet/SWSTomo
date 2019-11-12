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



class PathSetter:
    """A class to hold the metadata for the run (rdir, station? [for now], outdir etc. ) and fucntions
    to parse/generate the XML needed for the PathSet file
    Inputs: df [obj] - a DataFrame of the relevent .pairs file (Setter will select rows from the relevent station)
            station [str] - the station code for the station[s] we want to include data for (starting with a single station)
            ddir [str] - the data directory. Directory where sheba has output the .mts (and SAC) files to
            odir [str] - the output directory. Path to where we want our output. If none, we use the current working directory
    """
    def __init__(self,df_in,ddir,station=None,model=None,odir=None,config_uid='Test Run'):

        if odir == None:
            self.opath = os.getcwd() # Gets current working directory and stores is as a Path
            self.odir = self.opath.split('/')[-1]
        else:
            self.opath = '/Users/ja17375/SWSTomo/BluePebble/{}'.format(odir)

        if model == None:
            print('Using Model.xml from MTS_Setup, copying to cwd {}'.format(self.opath))
            self.modelxml = '/Users/ja17375/SWSTomo/MTS_Setup/Model.xml'
        else:
            self.modelxml = '/Users/ja17375/SWSTomo/MTS_Setup/{}'.format(model)

        if station == None:
            self.stations = df_in.STAT
            self.df = df_in
            print(self.stations)
        else:
            self.stations = station
            self.df = df_in[df_in['STAT'].isin(station)]

        self.ddir = ddir # Data directory (hopefully)
        self.dom_h = 250 # [km] height of the domains (fixed for UM and D'' for now!)
        self.xmlns = {'mtsML':'http://www1.gly.bris.ac.uk/cetsei/xml/MatisseML/'} # Dict with the Matisse XML namespace (for easy ref)
        ## Read Model.xml (assumed to be in working directory because why wouldnt it?)
        root = ElementTree.parse('{}'.format(self.modelxml)).getroot()
        self.model = root.find('mtsML:model',self.xmlns)
        model_name = self.model[0].text
        print('Using model named ... {}'.format(model_name))
        # Set config uID
        self.config_uid = config_uid

    def get_mts(self,phase):
        '''Function to get the .mts file for a phase and read in the xml.
           Phase [str] - SKS or SKKS (as my sheba runs are split by phase)
           fileID [str] - file name (.mts is appended here) [ now an attribute of Setter class]
        '''
        mts = '{}/{}/{}/{}.mts'.format(self.ddir,self.stat,phase,self.fileID)
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
                file = '{}/{}/{}/{}.BH{}'.format(self.ddir,self.stat,phase,self.fileID,comp)
                dst = '{}/data/{}.BH{}'.format(self.opath,self.fileID,comp)
                p = copy(file, dst)

    def domain2operator(self,domain):
        '''Function to read model.xml and extract the required domain
        Input -------
        domain [str] - the domain name (tag <domain_uid>) for the domain we want to find and "cast" as the operator
        '''
        domains = self.model.findall('mtsML:domain',self.xmlns)

        for dom in domains:
            uid_tmp = dom.find('mtsML:domain_uid',self.xmlns).text
            if uid_tmp == domain:
                print('Domain {} Found'.format(uid_tmp))
                uid = uid_tmp
                if uid.split('_')[0] == 'Lower':
                    # Domains starting with Lower are at CMB. So depth == 2890 km
                    depth = 2890. # Approx depth of CMB [km]
                    # azi = '0'
                elif uid.split('_')[0] == 'Upper':
                    depth = 250. # Depth of upper domain (keeping domain the same thickness for now)
                    # azi = str(self.az)
                else:
                    raise NameError('Domain is incorreclty named. Should be either "Upper" or "Lower".')

        # aoi = slw2aoi(depth,self.evdp,self.gcarc,phase) # Calculate ray param and then incidence angle
        aoi = 0 # Assume rays are vertical, not true but using this for testing.
        dist = self.dom_h / np.cos(np.radians(aoi)) # Calculate distance travelled through domain

        ## Now return what we need to make the opertor (will do this above?)
        operator = ElementTree.Element('operator')
        dom_uid = ElementTree.SubElement(operator,'domain_uid')
        dom_uid.text = domain
        azimuth = ElementTree.SubElement(operator,'azi')
        # azimuth.text = str(self.az)
        azimuth.text = '0'
        inclination = ElementTree.SubElement(operator,'inc')
        inclination.text = str(90 - aoi) # change from aoi to inclination
        l = ElementTree.SubElement(operator,'dist')
        l.text = str(dist)

        return operator

    def gen_PathSet_XML(self,phases=['SKS','SKKS'],fname=None):
        '''
        Function to iterate over the DataFrame and construct the XML we need
        Phases [list] - list of phase codes to iterate over (assuming each row in the DataFame corresponds to all phases)
        '''
        # Some quick i/o management
        if fname == None:
            self.pathset_xml = 'Pathset' # Set default Pathset XML file name to 'Pathset.xml' (the .xml is added later)
        else:
            self.pathset_xml = fname
        #Start of main function
        self.pathset_root = ElementTree.Element('MatisseML')
        tree = ElementTree.ElementTree(self.pathset_root)
        self.pathset_root.set("xmlns",self.xmlns['mtsML'])
        pathset = ElementTree.SubElement(self.pathset_root,'pathset')
        psuid = 'Paths for run in dir {} .'.format(self.odir)
        pathset_uid = ElementTree.SubElement(pathset,'pathset_uid')
        pathset_uid.text = psuid
        dom_used = []

        for stat in self.stations:
        # Loop over station list
            print(stat)

            self.stat = stat
            sdf = self.df[self.df.STAT ==stat]
            for i, row in sdf.iterrows():
                # All XML generation must sit within this loop (function calls) so that we make sure that Az, EVDP etc. are right for the current phase
                self.evdp = row.EVDP # read out event depth [km]. Attrib as needed for path length calcs.
                self.az = row.AZI
                self.gcarc = row.DIST
                evla = row.EVLA
                evlo = row.EVLO
                stla = row.STLA
                stlo = row.STLO

                for ph in phases:
                    #Each iteration of this loop is a seperate path (1 per phase and event)
                    f = '{}/{}/{}/{}_{}_{}??_{}.mts'.format(self.ddir,stat,ph,stat,row.DATE,row.TIME,ph)

                    try:
                        self.fileID = glob.glob(f)[0].strip('.mts').split('/')[-1] # Strip out .mts and split by '/', select end to get filestem
                    except IndexError:
                        print('File {} Not found'.format(f))
                        continue

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
                    # Assign Upper Domain
                    if row.STLA >=40. :
                        if row.STLO <= -110.0:
                            op_UM = self.domain2operator('Upper_01') # This will need to be a dictionary of domain UIDs
                        elif (row.STLO > -110.0) and (row.STLO <= -95.0):
                            op_UM = self.domain2operator('Upper_02')
                        elif row.STLO > -95.0:
                            op_UM = self.domain2operator('Upper_03')
                        else:
                            print("Error No Domain for this!")
                    elif row.STLA <= 40. :
                        if row.STLO <= -110.0:
                            op_UM = self.domain2operator('Upper_04') # This will need to be a dictionary of domain UIDs
                        elif (row.STLO > -110.0) and (row.STLO <= -95.0):
                            op_UM = self.domain2operator('Upper_05')
                        elif row.STLO > -95.0:
                            op_UM = self.domain2operator('Upper_06')
                        else:
                            print("Error No Domain for this STLO!")
                    else:
                        print("Error No Upper Domain for this STLA!")
                    # Assign Lower Domain
                    if ph == 'SKS':
                        # Now parameterised for E_Pacifc, test of lat/lon and assign domain accordingly
                        if row.SKS_PP_LAT >=35. :
                            if row.SKS_PP_LON <= -130.6:
                                op_LM = self.domain2operator('Lower_01') # This will need to be a dictionary of domain UIDs
                            elif (row.SKS_PP_LON > -130.6) and (row.SKS_PP_LON <= -110.0):
                                op_LM = self.domain2operator('Lower_02')
                            elif row.SKS_PP_LON > -110.0:
                                op_LM = self.domain2operator('Lower_03')
                            else:
                                print("Error No Domain for this!")
                        else:
                            op_LM = self.domain2operator('Lower_04')
                    # if 'Lower_{}'.format(dom_ext) not in dom_used:
                    #     dom_used.append('Lower_{}'.format(dom_ext))
                    elif ph == 'SKKS':
                        if row.SKKS_PP_LAT >=35. :
                            if row.SKKS_PP_LON <= -130.6:
                                op_LM = self.domain2operator('Lower_01') # This will need to be a dictionary of domain UIDs
                            elif (row.SKKS_PP_LON > -130.6) and (row.SKKS_PP_LON <= -110.0):
                                op_LM = self.domain2operator('Lower_02')
                            elif row.SKKS_PP_LON > -110.0:
                                op_LM = self.domain2operator('Lower_03')
                            else:
                                print("Error No Domain for this!")
                        else:
                            op_LM = self.domain2operator('Lower_04')
                    path.append(op_LM)
                    path.append(op_UM)
        print(dom_used)
        self._write_pretty_xml(self.pathset_root,file='{}/{}.xml'.format(self.opath,self.pathset_xml))

    def baz_test(self,baz):
        '''
        Function to test which 30 deg backzimuth bin each phase sits in
        '''
        for i,v in enumerate([30,60,90,120,150,180,210,240,270,300,330,360]):
            if (baz > v-30) and (baz <= (v)):
                return "{:03d}".format(v)

    def gen_MTS_Info(self):
        '''
        Function to create the simple XML file "MTS_Info"
        '''
        root = ElementTree.Element('MatisseML')
        tree = ElementTree.ElementTree(root)
        root.set("xmlns",self.xmlns['mtsML'])
        root.append(ElementTree.Comment('Generated by gen_MTS_Info from pathset.py'))
        samplerinfo = ElementTree.SubElement(root,'SamplerInfo')
        c_uid = ElementTree.SubElement(samplerinfo,'config_uid')
        c_uid.text = self.config_uid
        # Now write the XML to a file
        self._write_pretty_xml(root,file='{}/MTS_Info.xml'.format(self.opath))

    def gen_MTS_Config(self,options=None,model='Model.xml'):
        '''
        Function to create the MTS config file
        Inputs ===============
            Model   - The name of the Model file
            Option  - A dictionary of the operational options that one can pass to MTS.
                      If no dict is provided that the defaults will be set. calc_2d_mppad has to be provided as it depends on domain names
        ======================
        Pathset XML file name is taken from self.pathset_xml and is set by gen_PathSet_XML

        '''
        root = ElementTree.Element('MatisseML')
        tree = ElementTree.ElementTree(root)
        root.set("xmlns",self.xmlns['mtsML'])
        root.append(ElementTree.Comment('Generated by gen_MTS_Config from pathset.py'))
        config = ElementTree.SubElement(root,'config')
        c_uid = ElementTree.SubElement(config,'config_uid')
        c_uid.text = self.config_uid
        config.append(ElementTree.Comment(' input files'))
        mf = ElementTree.SubElement(config,'model_file') # Creat model file tag
        mf.text = model
        ps = ElementTree.SubElement(config,'pathset_file')
        ps.text = '{}.xml'.format(self.pathset_xml)
        config.append(ElementTree.Comment('operational options')) # Break before operational options ( for readablility)
        if options is None:
            print('Config not specified, using Defaults. Calc_2d_mmpad have to be added manually!')
            options = {'verbose': '0', 'max_ensemble_size': '10000000', 'tlag_domain':'time', 'nmodels':'2000000','nmod_burn': '10000',
                       'nmod_skip':'1', 'step_size': '0.05', 'term_mode':'convergence', 'conv_limit': '0.0005', 'nmod_stable': '1000',
                       'save_ensemble': '1', 'save_expvals': '1', 'nbin_1d_mppd': '50', 'nbin_2d_mppd': '50' }
        for opt in options:
            conf_SE = ElementTree.SubElement(config,opt)
            conf_SE.text=options[opt]
            if opt == 'tlag_domain':
                config.append(ElementTree.Comment(' MH Control Options '))
            elif opt == 'stepsize':
                config.append(ElementTree.Comment(' Termination Options '))
            elif opt == 'nmod_stable':
                config.append(ElementTree.Comment(' Output and Storage Options '))

        # tree.write('MTSConfig.xml')
        self._write_pretty_xml(root,file='{}/MTSConfig.xml'.format(self.opath))

    def _write_pretty_xml(self,root,file):
        """Return a pretty-printed XML string for the Element.
        """
        rough_string = ElementTree.tostring(root)
        reparsed = minidom.parseString(rough_string)
        print(reparsed)
        with open('{}'.format(file),'w') as writer :
            reparsed.writexml(writer,indent="    ",addindent="    ",newl="\n",encoding="utf-8")

        print('XML written to {}'.format(file))
