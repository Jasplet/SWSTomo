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
from sphe_trig import dist_deg
###############################
__author__ = "Joseph Asplet"
__license__ = "MIT"
__version__ = "0.1"
__email__ = "joseph.asplet@bristol.ac.uk"
__status__ = "Development"


class PathSetter:
    """A class to hold the metadata for the run (rdir, station? [for now], outdir etc. ) and fucntions
    to parse/generate the XML needed for the PathSet file
    Inputs:
            Required:
            #########
            df [obj]   -  filename for the relevent .pairs file, with T3 bins for SKS, SKKS in lower domains and in the upper mantle
                         Titled: SKS_BIN, SKKS_BIN, UPPER_BIN (Setter will select rows from the relevent station)
            ddir [str] - The data directory. Directory where sheba has output the .mts (and SAC) files to
            domains    - A domains file (e.g. E_pac.T3.counts.doms). This file must contain the trigonal domain midpoints (having Vertices is
                         nice, as it the Upper and Lower domains counts for each bin (i.e. how many paths were assigned to each bin by geogeom)
            Optional:
            station [str] - the station code for the station[s] we want to include data for (starting with a single station).
                            Now that we are looking at all stations in the East_Pacific we can read them in from a textfile
            odir [str]    - The output directory. Path to where we want our output. If none, we use the current working directory
            model [str]   - Name of the model file (assumed to be within MTS_Setup) that is to be used. If none is provided Model.xml is used
            config_uid [str] - An optional unique identifier that will be added to the MTSConfig file
    """
    def __init__(self,df_file,ddir,domains,station=None,model=None,odir=None,config_uid='Test Run',use_setup=False):

        print(df_file)
        print('Reading df')
        date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
        df_in = pd.read_csv(df_file,converters=date_time_convert,delim_whitespace=True)
        print('Read Trigonal Domians file')
        self.doms = pd.read_csv(domains,delim_whitespace=True)
        print('Set Outdir')
        if odir == None:
            self.opath = os.getcwd() # Gets current working directory and stores is as a Path
            self.odir = self.opath.split('/')[-1]
        else:
            self.opath = '/Users/ja17375/SWSTomo/BluePebble/{}'.format(odir)
        # Set XML namespace
        self.xmlns = {'mtsML':'http://www1.gly.bris.ac.uk/cetsei/xml/MatisseML/'}

        if model == None:
            print('Model will need to be generated before Pathsetting can be done')
            # self.modelxml = '/Users/ja17375/SWSTomo/MTS_Setup/Model.xml'
        elif use_setup is True:
            print("Reading Model file from MTS_Setup")
            self.modelxml = '/Users/ja17375/SWSTomo/MTS_Setup/{}'.format(model)
            ## Read Model.xml (assumed to be in working directory because why wouldnt it?)
            root = ElementTree.parse('{}'.format(self.modelxml)).getroot()
            self.model = root.find('mtsML:model',self.xmlns)
            model_name = self.model[0].text
            print('Using model named ... {}'.format(model_name))
        else:
            print("Reading Model file from {}".format(os.getcwd()))
            self.modelxml = '{}/{}'.format(self.opath,model)
            ## Read Model.xml (assumed to be in working directory because why wouldnt it?)
            self.modelroot = ElementTree.parse('{}'.format(self.modelxml)).getroot()
            self.model = self.modelroot.find('mtsML:model',self.xmlns)
            model_name = self.model[0].text
            print('Using model named ... {}'.format(model_name))

        if station == None:
            stat_check = pd.read_csv('/Users/ja17375/SWSTomo/Jacks_stats_in_Epac.txt',delim_whitespace=True)
            self.stations = stat_check.STA
            self.df = df_in
            # print(self.stations)
        else:
            self.stations = station
            self.df = df_in[df_in['STAT'].isin(station)]

        self.ddir = ddir # Data directory (hopefully)
        self.dom_h = 250 # [km] height of the domains (fixed for UM and D'' for now!)

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

    def domain2operator(self,domain,phase):
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
            # else:
            #     print('Domain {} is not in Model!'.format(uid_tmp))

        aoi = slw2aoi(depth,self.evdp,self.gcarc,phase) # Calculate ray param and then incidence angle
        # aoi = 0 # Assume rays are vertical, not true but using this for testing.
        dist = self.dom_h / np.cos(np.radians(aoi)) # Calculate distance travelled through domain

        ## Now return what we need to make the opertor (will do this above?)
        operator = ElementTree.Element('operator')
        dom_uid = ElementTree.SubElement(operator,'domain_uid')
        dom_uid.text = domain
        azimuth = ElementTree.SubElement(operator,'azi')
        azimuth.text = str(self.az)
        # azimuth.text = '0'
        inclination = ElementTree.SubElement(operator,'inc')
        inclination.text = str(90 - aoi) # change from aoi to inclination
        l = ElementTree.SubElement(operator,'dist')
        l.text = str(dist)

        return operator

    def parsedoms(self):
        '''

        '''

        dml = '{{{}}}domain'.format(self.xmlns['mtsML']) # need the triple curly braces to get a get of braces surrounding the mtsML
        udoms = []
        ldoms = []
        for dom in self.model.iter(dml):
            d_id = dom[0].text
            if d_id.split('_')[0] == 'Upper':
                udoms.append(int(d_id.split('_')[1]))
            elif d_id.split('_')[0] == 'Lower':
                ldoms.append(int(d_id.split('_')[1]))

        return (udoms,ldoms)


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

        k = 0
        nf = 0 # counter for files not found
        q_fail = 0 # counter for pahses that fail Q test
        ex = 0
        # Before Pathsetting, parse the upper/lower domains from the model file
        udoms,ldoms = self.parsedoms()

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
                    k = k + 1
                    #Each iteration of this loop is a seperate path (1 per phase and event)
                    f = '{}/{}/{}/{}_{}_{}??_{}.mts'.format(self.ddir,stat,ph,stat,row.DATE,row.TIME,ph)
                    if (ph == 'SKS') and (row.Q_SKS > 0.5) or (row.Q_SKS < -0.7):
                        print(' SKS Pass')
                        # pass # no operation required other than to keep the loop iterating
                    elif (ph == 'SKKS') and (row.Q_SKKS > 0.5) or (row.Q_SKKS < -0.7) :
                        print('SKKS Pass')
                        # pass # Phase is SKKS and event is a clear split or null, so we can keep going
                    else:
                        print('Phase {} fails Q tests, continuing to next'.format(ph))
                        q_fail += 1
                        continue  # SKS,SKKS phase is NOT a clear split or null, so we don't want to use it. continue to next iteration of loop
                    try:
                        self.fileID = glob.glob(f)[0].strip('.mts').split('/')[-1] # Strip out .mts and split by '/', select end to get filestem
                    except IndexError:
                        print('File {} Not found'.format(f))
                        nf += 1
                        continue
                    # Now we need to select the correct domains and in the right order (order of operators).
                    # As the model get more complex these tests will have to get more "clever"
                    # Hardcoded for now, need to get a function to read Model.xml (possibly as part of __init__)


                    # Now add Path to XML file
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
                    path.append(op_LM)
                    path.append(op_UM)

        in_u_dom = u1+u2+u3+u4+u5+u6
        print(dom_used)
        print('Total Phases:', k)
        print('Upper_01: {}, Upper_02: {}, Upper_03: {}, Upper_04: {}, Upper_05: {}, Upper_06: {}'.format(u1,u2,u3,u4,u5,u6))
        print('Lower_01: {}, Lower_02: {}, Lower_03: {}, Lower_04: {}'.format(l1,l2,l3,l4))
        print('Phases assigned an upper domain:', in_u_dom)
        print('Phases not assigned a lower domain ', ex)
        print('Fail Q tests:', q_fail)
        print('Not Found:', nf)
        print('Sanity test [0 if sane]: {}'.format(k-(in_u_dom+q_fail+nf)))
        self._write_pretty_xml(self.pathset_root,file='{}/{}.xml'.format(self.opath,self.pathset_xml))

    def baz_test(self,baz):
        '''
        Function to test which 30 deg backzimuth bin each phase sits in
        '''
        for i,v in enumerate([30,60,90,120,150,180,210,240,270,300,330,360]):
            if (baz > v-30) and (baz <= (v)):
                return "{:03d}".format(v)

    def bin2domain(self,dom_name):
        '''
        Takes a single bin and creates the requisit domain XML (for Model.xml)
        '''
        domain = ElementTree.Element('domain')
        domain.append(ElementTree.Comment(' Identifier '))
        dom_uid = ElementTree.SubElement(domain,'domain_uid')
        dom_uid.text = dom_name
        domain.append(ElementTree.Comment(' Elastic medium '))
        medium = ElementTree.SubElement(domain,'medium')
        medium.text = 'elliptical:2000.,1000.,2000.'
        domain.append(ElementTree.Comment(' Inversion Parameters '))
        alpha = ElementTree.SubElement(domain, 'alpha', type="fixed",value="90")
        beta = ElementTree.SubElement(domain, 'beta', type="fixed",value="0")
        gamma = ElementTree.SubElement(domain, 'gamma', type="periodic",min="-90",max="90",init="0")
        domain.append(ElementTree.Comment('<gamma  type="fixed" value="-30"/>'))
        s = ElementTree.SubElement(domain, 'strength', type="fixed",min="0.00",max="0.05",init="0.0125")
        return domain

    def gen_Model_XML(self,mod_name=None,Low_Domains=None):
        '''
        Generate a default Model XML file, containing all requested Lower Domains (default is all in .domains file)
        and all corresponding Upper domains.
        '''
        if Low_Domains is None:
            print('Using all domains in .domains file')
            lowb = self.df.SKS_BIN.append(self.df.SKKS_BIN)
            Low_Domains = lowb.drop_duplicates().values

        if mod_name is None:
            mod_name = input('No model name provided. Enter one now :')

        #find corresponding Upper Domains
        # first test outside of loop so we can append the rest into 1 Series more easily
        u_df = self.df[(self.df['SKS_BIN'] == Low_Domains[0]) | (self.df['SKKS_BIN'] == Low_Domains[0])]['UPPER_BIN']
        for i in range(1,len(Low_Domains)):
            tmp = self.df[(self.df['SKS_BIN'] == Low_Domains[i]) | (self.df['SKKS_BIN'] == Low_Domains[i])]['UPPER_BIN']
            u_df = u_df.append(tmp)
        Up_Domains = u_df.drop_duplicates()


        # Write the xml
        root = ElementTree.Element('MatisseML')
        tree = ElementTree.ElementTree(root)
        root.set("xmlns",self.xmlns['mtsML'])
        root.append(ElementTree.Comment('Generated by gen_Model from pathset.py'))
        m = ElementTree.SubElement(root,'model')
        m_uid = ElementTree.SubElement(m,'model_uid')
        m_uid.text = 'Trigonal Domains'

        u,l = 0,0
        for udom in Up_Domains:
            dom = self.bin2domain('Upper_{}'.format(udom))
            m.append(dom)
            u += 1

        for ldom in Low_Domains:
            # Loop over requested lower domains
            dom = self.bin2domain('Lower_{}'.format(ldom))
            m.append(dom)
            l+=1

        print('Model Generated')
        self.model = m
        print('{} Upper domains, {} Lower domains'.format(u,l))
        self._write_pretty_xml(root,file='{}/{}.xml'.format(self.opath,mod_name))


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

    def gen_MTS_Config(self,options=None):
        '''
        Function to create the MTS config file
        Inputs ===============
            Model   - The name of the Model file
            Option  - A dictionary of the operational options that one can pass to MTS.
                      If no dict is provided that the defaults will be set. calc_2d_mppad has to be provided as it depends on domain names
        ======================
        Pathset XML file name is taken from self.pathset_xml and is set by gen_PathSet_XML

        '''
        model = self.modelxml.split('/')[0]

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
