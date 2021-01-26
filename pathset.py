'''
This program is designed to gather to generate the required Pathset.xml file.
Primarily by collecting the .mts files that are output by sheba when measuring SWS, for each phase a seperate <data> tag is required
Stating the paths to the requisit SAC files (which are copied to the data directory in the current model path [for now at least])
#
'''
######## Imports ##############
from xml.etree import ElementTree 
# ElementTree is a standard (as of python 2.5) library which we can use to parse XML.
# N.B ET is NOT secure against "malicous XML" however as we are only intersted in very simple XML this shouldnt be an issue
from xml.dom import minidom
import pandas as pd
import os
import glob
import numpy as np
from sphe_trig import vincenty_dist
from calc_aoi import slw2aoi,get_rayparam
from sactools import get_sac,get_mts
from corrections import add_sside_correction, add_rside_correction,bin2domain
###############################
__author__ = "Joseph Asplet"
__license__ = "MIT"
__version__ = "0.1"
__email__ = "joseph.asplet@bristol.ac.uk"
__status__ = "Development"

class PathSetter:
    """A class to hold the metadata for the run (rdir, station? [for now], outdir etc. ) and fucntions
    to parse/generate the XML needed for the PathSet file
    """
    def __init__(self,phasefile,phases,domains,model=None,odir=None,config_uid='Test Run'):
        '''
        Initialise Pathsetter with some essential data, metadata. 

        Args:
            phasefile (str) - file that contains splitting results for all phases
            phases (list) - list of phases that are to be used
            domains (str) - file that contains information on the size/ shape of the domains
            model (str) - model XML file.If the model file does not exists Pathsetter can create one
            odir (str) - the output directory that XML files are written to
            config_uid (str) - a unique identifier written to MTSConfig.xml, default is 'Test Run'

        Returns:
            Set (obj) - a configured Pathsetter Object.

        Examples:
            >>> Set = Pathsetter('E_pacific_05_binned.pairs',['SKS','SKKS'],
                                'ScS_w_bouncepoints_binned_05.sdb','ScS','TestModel.xml')                       
        '''
        
        print('Reading df')
        self.date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
        self.df = pd.read_csv(phasefile,converters=self.date_time_convert,delim_whitespace=True)
        self.phases = phases
        print('Phases specified are: ', self.phases)
        print('Read Trigonal Domians file')
        self.doms = pd.read_csv(domains,delim_whitespace=True).sort_values(by='BIN')
        print('Set Outdir')
        if odir == None:
            self.opath = os.getcwd() # Gets current working directory and stores is as a Path
            self.odir = self.opath.split('/')[-1]
        elif len(odir.split('/')) > 1:
            self.opath = odir
            self.odir = odir.split('/')[-1]
        else:
            self.opath = '/Users/ja17375/SWSTomo/BlueCrystal/{}'.format(odir)
            self.odir = odir
        # Set XML namespace
        self.xmlns = {'mtsML':'http://www1.gly.bris.ac.uk/cetsei/xml/MatisseML/'}
        if model == None:
            print('Model will need to be generated before Pathsetting can be done')
        else:
            print("Reading Model file from SWSTomo/Models")
            self.modelxml = '/Users/ja17375/SWSTomo/Models/{}'.format(model)
            self.modelroot = ElementTree.parse('{}'.format(self.modelxml)).getroot()
            self.model = self.modelroot.find('mtsML:model',self.xmlns)
            model_name = self.model[0].text
            print('Using model named ... {}'.format(model_name))
        # Set config uID
        self.config_uid = config_uid

    def domain2operator(self,domain,phase,evdp,gcarc,azi):
        '''Function to read model.xml and extract the required domain and generate the operator XML needed for the Pathset file
        
        Args:
            domain (str) - the domain uid for the domain to make the operator for
            phase (str) - the phase passing through the domain
            evdp (str) - source event depth for phase
            gcarc (float) - event-station distance [deg]
            azi (float) - event station azimuth
        Returns:
            operator (obj) - An ElementTree Strucutre conainting the operator XML
        '''
        domains = self.model.findall('mtsML:domain',self.xmlns)
        uid = None
        for dom in domains:
            # Loop through domains in the model file 
#             uid_tmp = dom.find('mtsML:domain_uid',self.xmlns).text
            uid_tmp = dom[0].text
            if domain == uid_tmp:
                uid = uid_tmp
                break
        
        if uid:
            print('Domain {} Found'.format(uid_tmp))
        else:
            raise ValueError('Domain {} Not Found'.format(domain))
           
                
        if uid.split('_')[0] == 'Lower':
            # Domains starting with Lower are at CMB. So depth == 2890 km
            depth = 2890. # Approx depth of CMB [km] 

            slw = get_rayparam(evdp,gcarc,phase)
            aoi = slw2aoi(depth,slw) # Calculate ray param and then incidence angle
            dom_h = 250
            dist = dom_h / np.cos(np.deg2rad(aoi))
        elif (uid.split('_')[0] in ['Upper','RSide','Station']):
            depth = 100. # 100km so we can use conversion dt/100 = strength

            slw = get_rayparam(evdp,gcarc,phase)
            aoi = slw2aoi(depth,slw) # Calculate ray param and then incidence angle
            dist = 100
        elif uid.split('_')[0] == 'SSide':
            # If domain ID is SSide then we know this domain is for a source side correction. 
            azi = 0 # For the splitting corrections we fix azi to 0
            aoi = 0
            dist = 100
            print('SSide')
        else:
            raise ValueError(f'Domain "{uid}"name does not match expected values')
        
        print(dist)
        operator = ElementTree.Element('operator')
        dom_uid = ElementTree.SubElement(operator,'domain_uid')
        dom_uid.text = domain
        azimuth = ElementTree.SubElement(operator,'azi')
        azimuth.text = str(azi)
        # azimuth.text = '0'
        inclination = ElementTree.SubElement(operator,'inc')
        inclination.text = str(90 - aoi) # change from aoi to inclination
        l = ElementTree.SubElement(operator,'dist')
        l.text = str(np.around(dist,decimals=3))
        return operator

    def parsedoms(self):
        '''
        Function to parse the input model XML into lists of domains for use in other functions
        
        Args:
            None - model XML is read in __init__ and is an attribute of the Setter Class (self)
        Returns:
            rdoms (list) - a list of all upper mantle reciever side domains in the model
            ldoms (list) - a list of all D`` domains in the model
            sdoms (list) - a list of all upper mantle source side domains in the model
        Examples:
            >>>(rdoms, ldoms, sdoms) = Set.parsedoms()
        '''
        dml = '{{{}}}domain'.format(self.xmlns['mtsML']) 
        # need the triple curly braces to get a set of braces surrounding the mtsML
        u = []
        l = []
        sdoms = []
        self.udom_c = {}
        self.ldom_c = {}
        self.sdom_c = {}
        for dom in self.model.iter(dml):
            d_id = dom[0].text
            if d_id.split('_')[0] == 'Station':
                u.append(d_id)
                self.udom_c.update({d_id : 0})
            elif d_id.split('_')[0] == 'Lower':
                l.append(d_id)
                self.ldom_c.update({d_id : 0})
            elif d_id.split('_')[0] == 'SSide':
                # If the uid is SSide we know it is a source side correction domain. We just need to 
                # match the correct domain to the correct phase 
                sdoms.append(d_id)
                self.sdom_c.update({d_id : 0 })
                
        ldoms = self.doms[self.doms['BIN'].isin([i.split('_')[1] for i in l ])]
        rdoms = self.doms[self.doms['BIN'].isin([i.split('_')[1] for i in u ])]        
        return (rdoms,ldoms,sdoms)

    def is_phase_good(self,row,split_c=0.7, null_c=-0.9):
        '''
        Function to test if a phase is not split or null based on its Q value (Wuestefeld et al., 2010)
        
        Args:
            row (df) - a row from a pandas DataFrame containing data for multiple phases
            ph (str) - the phase to do the Q test for
        
        Returns:
            res (bool) - The result of the Q test. True if phase passed, False if it fails.
        '''
        if (row['Q'] < 0.7) and (row['Q'] > -0.9):
            res = False
        else:
            res = True
        return res

    def loop_thru_domains(self,doms,raylat,raylon,uid,ph,crit,evdp,gcarc,azi):
        ''' 
        This function loops through a set of input domains and identifies domains whose midpoints are within a set distance from a raypath (the point the raypath is at at a given depth in the mantle). For these domains it then creates the required operator XML sequence
        
        Args:
            doms (list) - list of domains to loop over
            raylat (float) - latitude of the raypath
            raylon (float) - longitude of the raypath
            uid (str) - domain uid, for creating the XML (e.g. 'Upper')
            ph (str) - the phase ID ('SKS','SKKS','ScS')
            crit (float) - the distance ciriterea we want to find domains for
            evdp (str) - source event depth for phase
            gcarc (float) - event-station distance [deg]
            azi (float) - event station azimuth
            
        Returns:
            operators (list) - a list containing operator objects from domain2operator for the domains that were within the test distance from the raypath. 
'''
        operators = [ ]
        for d,row in doms.iterrows():
            domlat = row.MID_LAT
            domlon = row.MID_LON
            dom_id = "{}_{:.0f}".format(uid,row.BIN)
            dist, azi = vincenty_dist(raylat,raylon,domlat,domlon)
            if dist <= crit:
                op = self.domain2operator(dom_id,ph,evdp,gcarc,azi)
                operators.append(op)
                                                
        return operators                                      
                                                
    def gen_PathSet_XML(self,phases=['ScS','SKS','SKKS'],fname=None,crit=5.0):
        '''
        Function to iterate over the shear-wave phase data provided to Setter and to generate the correct XML to from the Pathset XML file required for Matisse. This includes adding source/ reciever side corrections where required (and where the neccesary corrections have been provided as needed by functions called by gen_PathSet_XML)
        
        Args:
            phases [list] - list of phase codes to iterate over. Default is ['ScS','SKS','SKKS']
            fname [str] - File name to write the Pathset XML to. If no name is provided then this defaults to 'Pathset'
            
        Returns:
            Output XML is written to a textfile
        '''
        # Some quick i/o management
        if fname == None:
            self.pathset_xml = 'Pathset' 
            # Set default Pathset XML file name to 'Pathset.xml' (the .xml is added later)
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
   
        p = 1
        nf = 0 # counter for files not found
        q_fail = 0 # counter for pahses that fail Q test
        # Before Pathsetting, parse the upper/lower domains from the model file
        rdoms,ldoms,sdoms = self.parsedoms()
        # [deg] - the distance criterea for including phases in a domain.
                   #         designed to give some overlap in neighbouring domians for reduce edge effects...
        ## Wrangle the SnKS and ScS dataframes to produce one combined dataframe we can then iterate over

        for i, row in self.df.iterrows():
            # All XML generation must sit within this loop (function calls) so that we make sure that Az, EVDP etc. are right for the current phase
            attribs = {'evdp':row.EVDP,'azi':row.AZI,'gcarc':row.DIST,'evla':row.EVLA,
                            'evlo':row.EVLO,'stla':row.STLA,'stlo':row.STLO,'stat':row.STAT,
                            'date':row.DATE,'time':row.TIME
                            }
            for ph in phases:
                if ph == 'ScS':
                    f = '/Users/ja17375/DiscrePy/Sheba/Runs/ScS/{}/{}/{}_{}_{}??_{}.mts'.format(
                        attribs['stat'],ph,attribs['stat'],attribs['date'],attribs['time'],ph)
                else:
                    f = '/Users/ja17375/DiscrePy/Sheba/Runs/E_pacific/{}/{}/{}_{}_{}??_{}.mts'.format(
                        attribs['stat'],ph,attribs['stat'],attribs['date'],attribs['time'],ph)

                ph_good = self.is_phase_good(row,ph) 
                # Test if phases are "good" splits or nulls
                try:
                    fileID = glob.glob(f)[0].strip('.mts').split('/')[-1]
                    # Strip out .mts and split by '/', select end to get filestem
                except IndexError:
                    nf += 1
                    # print('No file')
                    continue
                if ph_good is True:     
                    # print('Pass')
                    llat = row['LOWMM_LAT']
                    llon = row['LOWMM_LON']
                    lowmm_ops = self.loop_thru_domains(ldoms,llat,llon,"Lower",ph,crit,
                                                       attribs['evdp'],attribs['gcarc'],attribs['azi'])
                    rside_ops = self.loop_thru_domains(rdoms,attribs['stla'],attribs['stlo'],"RSide",
                                                       ph,crit,attribs['evdp'],attribs['gcarc'],attribs['azi'])
                    sside_ops = [0] # placeholder, if not ScS we still need this so we iterate once over the doms
                    if ph == 'ScS': # Only ScS needs source side domains
#                   We have source side correction domians for each ScS phase. 
#                   Now we need to find the right correction for the current phase in the loop.
#                   As we have precalculated the corrections, a domain should exist for each ScS phase. 
#                   These corrections domains are identified by Station Code and event Date/Time
                        dom_id = 'SSide_{}_{}_{}'.format(attribs['stat'],attribs['date'],attribs['time'])
                        if np.isin(sdoms, dom_id):
                            sside_ops[0] = self.domain2operator(dom_id,'ScS',attribs['evdp'],attribs['gcarc'],
                                                 attribs['azi'])
                        else:
                            print(f'{dom_id} not in Domain, skip')
                            continue # move loop on as this ScS phase does not pass through the model domains   

                    for s,s_op in enumerate(sside_ops):
                        for m,l_op in enumerate(lowmm_ops):
                            for k,r_op in enumerate(rside_ops):
                                get_sac(fileID,attribs['stat'],ph)
                                # Now make XML for this Path
                                path = ElementTree.SubElement(pathset,'path')
                                pathname = 'Path {} {}'.format(p,ph)               
                                p +=1 
                                path_uid = ElementTree.SubElement(path,'path_uid')
                                path_uid.text = pathname
                                # Add Data (from .mts)
                                data = get_mts(fileID,attribs['stat'],ph)
                                path.append(data)
                                stat_uid = ElementTree.SubElement(path,'station_uid')
                                stat_uid.text = attribs['stat']
                                evt_uid = ElementTree.SubElement(path,'event_uid')
                                evt_uid.text = '{}_{}'.format(attribs['date'],attribs['time'])
                                if ph == 'ScS':
                                    print('Add SSide Operator to path')
                                    print(ph, s)
                                    path.append(s_op)
                                print('Rside Op')
                                path.append(l_op) 
                                path.append(r_op)

                else:
                    q_fail +=1
        #Now write out the Pathset XML
        self._write_pretty_xml(self.pathset_root,file='{}/{}.xml'.format(self.opath,self.pathset_xml))

    def gen_Model_XML(self,mod_name=None,Low_Domains=None,Rside_Domains=None,Sside_Domains=None,Rside_corr=True):
        '''
        Generates the model XML file required by matisse, based off the input domains. Can configure reciever side domains with corrections for phi based off surface wave models (Schaffer et al.). It is assumed that the domains used are T3 trigonal domains defined from the geogeom routine (written by J.Wookey). Domains should be assigned an identifying number, which in theory should work for other domain set-ups but this has not bee tested (as of 15/6/20)
        
        Args:
            mod_name (str) - name of the model. The model XML will be written to a file mod_name.xml
            Low_Domains (array-like) - [optional] array of D`` domains to use in model. If none, all domains that have some phases in the input data to Setter are used (IF you have used geogeom to bin the data already)
            Rside_Domains (array-like) - [optional] array of all reciever-side domains to use in model. If none, reciever side for all phases (from the input data to the Setter class) assigned to a D`` domain will be added
            Sside_Domains (array-like) - [optiona] array of all source-side domains to use in the model. These domaisn are only needed for ScS phases. If none, it is assumed that there are sourceside corrections to be added for any ScS phases. These take the form of one domain per phase, with the source side correction fixed in gamma and strength. 
            Rside_corr (bool) - [default = True]. Switch for whether receiver side corrections should be added or not. Reciever side corrections should be pre-calculated with one correction per domain. Our preferred source for these corrections is the surface-wave anisotropy model of Schaffaer et al., 2016 SL2016svA models.
            
        Returns:
            model XML written to the file mod_name.xml
        '''
        if mod_name is None:
            mod_name = input('No model name provided. Enter one now :')
      
        s,l,r = 0,0,0
    
        # Write the initial xml
        root = ElementTree.Element('MatisseML')
        tree = ElementTree.ElementTree(root)
        root.set("xmlns",self.xmlns['mtsML'])
        root.append(ElementTree.Comment('Generated by gen_Model from pathset.py'))
        m = ElementTree.SubElement(root,'model')
        m_uid = ElementTree.SubElement(m,'model_uid')
        m_uid.text = 'Trigonal Domains'
        
        if Low_Domains is None:
            print('Using all domains in .domains file')
            Low_Domains= self.df.LOWER_BIN.drop_duplicates().values

        if Rside_Domains is None:
            print('Finding all Reciever Side Domains')
            rside = rdf.drop_duplicates().values
        else:
            rside = Rside_Domains

        for i,row in self.df.iterrows():
            if (row.LOWER_BIN in Low_Domains) & (row.PHASE == 'ScS'):
                print(row.LOWER_BIN)
                date, time, stat = row.DATE, row.TIME, row.STAT
                sdom = add_sside_correction(date, time, stat)
                m.append(sdom)
                s += 1
                    
        if Rside_corr == True:
            for rdom in rside:
                rcorr = add_rside_correction(rdom)
                m.append(rcorr)
                r += 1
        # elif sside:
        #     ud = np.concatenate((rside,sside))
        #     udoms = np.unique(ud)
        #     udoms.sort()
        #     for udom in udoms:
        #         dom = bin2domain('Upper_{}'.format(int(udom)))
        #         m.append(dom)
        #         r += 1
        else:
            for udom in rside:
                dom = bin2domain('RSide_{}'.format(int(udom)))
                m.append(dom)
                r += 1

        for ldom in Low_Domains:
            # Loop over requested lower domains
            dom = bin2domain('Lower_{}'.format(int(ldom)))
            m.append(dom)
            l+=1

        print('Model Generated')
        self.model = m
        print('{} Upper (Sside) domains, {} Lower domains, {} Upper (Rside) domains'.format(s,l,r))
        self._write_pretty_xml(root,file='{}/{}.xml'.format(self.opath,mod_name))

    def gen_stat_PathSet_XML(self,stations, lower_dom_no=None, phases=['ScS','SKS','SKKS'],fname=None,crit=5.0,stat_only=True):
        '''
        Function to iterate over the shear-wave phase data provided to Setter and to generate the correct XML to from the Pathset XML file required for Matisse. This includes adding source/ reciever side corrections where required (and where the neccesary corrections have been provided as needed by functions called by gen_PathSet_XML)
        
        Args:
            stations [list] - list of stations to include data for
            lower_dom_no [int or list] - lower domain of interest
            phases [list] - list of phase codes to iterate over. Default is ['ScS','SKS','SKKS']
            fname [str] - File name to write the Pathset XML to. If no name is provided then this defaults to 'Pathset'
            crit (float) - the distance ciriterea we want to find domains for
            stat_only [bool] - flag for whether we want to add any D'' domains or not
        Returns:
            Output XML is written to a textfile
        '''
        # Some quick i/o management
        if fname == None:
            self.pathset_xml = 'Pathset' 
            # Set default Pathset XML file name to 'Pathset.xml' (the .xml is added later)
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
        p = 0
        # Before Pathsetting, parse the upper/lower domains from the model file
        rdoms,ldoms,sdoms = self.parsedoms()
        print(ldoms)
        df = self.df[self.df.STAT.isin(stations)]
        
        for i, row in df.iterrows():
            # All XML generation must sit within this loop (function calls) so that we make sure that Az, EVDP etc. are right for the current phase
            attribs = {'evdp':row.EVDP,'azi':row.AZI,'gcarc':row.DIST,'evla':row.EVLA,
                            'evlo':row.EVLO,'stla':row.STLA,'stlo':row.STLO,'stat':row.STAT,
                            'date':row.DATE,'time':row.TIME,'stat':row.STAT
                            }
            ph = row.PHASE
            if ph == 'ScS':
                f = '/Users/ja17375/DiscrePy/Sheba/Runs/ScS/{}/{}/{}_{}_{}??_{}.mts'.format(
                    attribs['stat'],ph,attribs['stat'],attribs['date'],attribs['time'],ph)
            else:
                f = '/Users/ja17375/DiscrePy/Sheba/Runs/E_pacific/{}/{}/{}_{}_{}??_{}.mts'.format(
                    attribs['stat'],ph,attribs['stat'],attribs['date'],attribs['time'],ph)

#             ph_good = self.is_phase_good(row,ph) 
            # Test if phases are "good" splits or nulls
            try:
                fileID = glob.glob(f)[0].strip('.mts').split('/')[-1]
                # Strip out .mts and split by '/', select end to get filestem
            except IndexError:
                continue
#             if ph_good is True:     
            statdom = 'Station_{}'.format(attribs['stat'])
            stat_op = self.domain2operator(statdom,ph,attribs['evdp'],
                                               attribs['gcarc'],attribs['azi'])
            if ph == 'ScS': # Only ScS needs source side domains
                dom_id = 'SSide_{}_{}_{}'.format(attribs['stat'],
                                                 attribs['date'],attribs['time'])
                print(dom_id)
                sside_op = self.domain2operator(dom_id,'ScS',attribs['evdp'],attribs['gcarc'],
                                     attribs['azi'])    
                print(sside_op)
            lowmm_ops = [0]
            if stat_only:
                lowmm_ops = [False]
            else:
                ldom_id = 'Lower_{}'.format(lower_dom_no)
                dlat = ldoms.MID_LAT.values[0]
                dlon = ldoms.MID_LON.values[0]
                print(dlon)
                dist = vincenty_dist(row['LOWMM_LAT'], row['LOWMM_LON'],dlat, dlon )[0]
                if dist < crit:
                    print('{} {}'.format(attribs['stat'],ph))
                    lowmm_ops[0] = self.domain2operator(ldom_id,ph,attribs['evdp'],
                                           attribs['gcarc'],attribs['azi'])
                else:
                    print('Path for {} {} too far from domain {} deg'.format(
                    attribs['stat'],ph, dist))
                    continue

            for m,l_op in enumerate(lowmm_ops):
                get_sac(fileID,attribs['stat'],ph)
                # Now make XML for this Path
                path = ElementTree.SubElement(pathset,'path')
                pathname = 'Path {} {}'.format(p,ph)          
                path_uid = ElementTree.SubElement(path,'path_uid')
                path_uid.text = pathname
                # Add Data (from .mts)
                data = get_mts(fileID,attribs['stat'],ph)
                path.append(data)
                stat_uid = ElementTree.SubElement(path,'station_uid')
                stat_uid.text = attribs['stat']
                evt_uid = ElementTree.SubElement(path,'event_uid')
                evt_uid.text = '{}_{}'.format(attribs['date'],attribs['time'])
                if ph == 'ScS':
                    print('Add SSide Operator to path')
                    path.append(sside_op)
                if l_op:
                    print("Add D`` Operator")
                    path.append(l_op) 
                path.append(stat_op)
                p+=1
        #Now write out the Pathset XML
        self._write_pretty_xml(self.pathset_root,file='{}/{}.xml'.format(self.opath,self.pathset_xml))

        
    def gen_stat_Model_XML(self, stations, mod_name=None, Low_Domains=None, ):
        '''
        Generates the model XML file required by matisse, with lowermost mantle domains and seperate upper mantle domains for every station requested. Source Side corrections for ScS are auto-added
        
        Args:
            stations (lsit) - list of all Stations to create domains for. 
            mod_name (str) - name of the model. The model XML will be written to a file mod_name.xml
            Low_Domains (array-like) - [optional] array of D`` domains to use in model. If none we assume this is a single domain run 
            
            
        Returns:
            model XML written to the file mod_name.xml
        '''
        if mod_name is None:
            mod_name = input('No model name provided. Enter one now :')
        modpath = '/Users/ja17375/SWSTomo/Models'
        # Write the initial xml
        root = ElementTree.Element('MatisseML')
        tree = ElementTree.ElementTree(root)
        root.set("xmlns",self.xmlns['mtsML'])
        root.append(ElementTree.Comment('Generated by gen_Model from pathset.py'))
        m = ElementTree.SubElement(root,'model')
        m_uid = ElementTree.SubElement(m,'model_uid')
        m_uid.text = 'Trigonal Domains'
        
        if Low_Domains is None:
            print('No Lowermost Mantle Domains')
        
        for stat in stations:
            dom = bin2domain('Station_{}'.format(stat))
            m.append(dom)
            
        for i,row in self.df.iterrows():
            if Low_Domains is None:
                if (row.STAT in stations) & (row.PHASE == 'ScS'):
                    date, time, stat = row.DATE, row.TIME, row.STAT
                    sdom = add_sside_correction(date, time, stat)
                    m.append(sdom)
            else:
                print(row.LOWER_BIN)
                if (row.STAT in stations) & (row.PHASE == 'ScS'):
                    date, time, stat = row.DATE, row.TIME, row.STAT
                    sdom = add_sside_correction(date, time, stat)
                    m.append(sdom) 
                
        if Low_Domains:
            for ldom in Low_Domains:
                # Loop over requested lower domains
                print(ldom)
                dom = bin2domain('Lower_{}'.format(int(ldom)))
                m.append(dom)

        print('Model Generated')
        self.model = m
        self._write_pretty_xml(root,file='{}/{}.xml'.format(modpath,mod_name))        
        
    def count_paths_in_model(self,Low_Domains=None,Rside_Domains=None):
        '''
        This function takes a test set of D'' and Reciever Side UM domains (default is all available domains) and counts the number of ScS, SKS and SKKS phases that are within the criteria distance from each of them. This is so that we can get an idea of how many paths are sampling each domain before we generate all the XML. This function repeats some things done in gen_Pathset and gen_Model, but here we do not have to worry about any XML. All we want to output is an array of domains and the number of paths for each phase that pass through the domain
        
        Args:
            Low_Domains (array-like) - [optional] array of D`` domains to use in model. If none, all domains that have some phases in the input data to Setter are used (IF you have used geogeom to bin the data already)
            Rside_Domains (array-like) - [optional] array of all reciever-side domains to use in model. If none, reciever side for all phases (from the input data to the Setter class) assigned to a D`` domain will be added.
            
        Returns:
            counts (array-like) - a 2d numpy array (shape(1280,7)) that contains the counts for the number of phases that pass through each D`` and reciever side domain requested. 
        '''
        date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
        df = pd.read_csv('E_pacific_SNR10_goodQ.sdb',converters=date_time_convert,delim_whitespace=True)
        crit = 5.0 # [deg] the criteria distance 
        if Low_Domains is None:
            ldoms = self.doms.BIN.values
        else:
            ldoms = np.asarray(Low_Domains)

        if Rside_Domains is None:
            rdoms = self.doms.BIN.values
        else:
            rdoms = np.asarray(Rside_Domains)
     
        doms = self.doms[self.doms.BIN.isin(ldoms) | self.doms.BIN.isin(rdoms)]
        print('Counting paths in domains {}'.format(doms.BIN.values))
        
        counts = np.zeros([1280,4])
#         countsM = np.zeros([1280,1280])
        # initialise counts array for all 1280 domains. Cols [bin, ScS_Low, SKS_Low, SKKS_Low, ScS_R, SKS_R, SKKS_R]

        for i,dom in doms.iterrows():
            for i, row in df.iterrows():
                dlat, dlon = dom.MID_LAT, dom.MID_LON
                dom_id = int(dom.BIN)
                dist_lowmm = vincenty_dist(row.LOWMM_LAT, row.LOWMM_LON, dlat, dlon)[0]
#                 dist_skks = vincenty_dist(row.SKKS_LM_LAT, row.SKKS_LM_LON, dlat, dlon)[0]
                dist_rside = vincenty_dist(row.STLA, row.STLO, dlat, dlon)[0]
                if (dist_lowmm <= crit) & (dom.BIN in ldoms):
                    if row.PHASE == 'ScS':
                        counts[dom_id-1,0] = dom_id
                        counts[dom_id-1,1] += 1
                    elif (row.PHASE == 'SKS') | (row.PHASE == 'SKKS'):
                        counts[dom_id-1,0] = dom_id
                        counts[dom_id-1,2] += 1
                    else:
                        raise ValueError('Phase code {} incorrect'.format(row.phase))                    
#                 if (dist_rside <= crit) & (dom.BIN in rdoms):
#                     print('RSIDE DOM: ',dom.BIN)
#                     if row.PHASE == 'ScS':
#                         counts[dom_id-1,0] = dom_id
#                         counts[dom_id-1,4] += 1
#                     elif (row.PHASE == 'SKS') | (row.PHASE == 'SKKS'):
#                         counts[dom_id-1,0] = dom_id
#                         counts[dom_id-1,5] += 1

        mask = np.all(np.equal(counts, 0), axis=1)                
        c_out = counts[~mask]
        print('Done')
        return c_out
                        
                        
    
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
        Function to create the MTS config file. N.B this function will just generate a default config file. It is easier to edit/copy an existing file as only a few changes are needed to tweak a run.
        
        Args:
            options (dict) - options to be set in the Config file. N.B all options must be specified

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

    def _write_domain_counts(self,file):
        """
        Adds the domains phase counts to the domains DF and then writes it out
        """
        self.doms['Lower_Count'] = self.doms['BIN'].map(self.ldom_c)
        self.doms['Upper_Count'] = self.doms['BIN'].map(self.udom_c)
        self.doms = self.doms.fillna(0)
        self.doms.to_csv(file,sep=' ',index=False)

def phases2sdb():
    '''
    Function to take data from multiple phases and make in into a single file to be used in Pathsetting
    '''
    date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
    df = pd.read_csv('E_pacific_05.pairs',converters=date_time_convert,delim_whitespace=True)
    df_snks = pd.read_csv('E_pacific_05_binned.pairs',converters=date_time_convert,delim_whitespace=True)
    cols = ['TLAG_SKS','DTLAG_SKS','FAST_SKS','DFAST_SKS','TLAG_SKKS','DTLAG_SKKS','FAST_SKKS',
            'DFAST_SKKS','SNR_SKS','SNR_SKKS']
    df_snks[cols] = df[cols]
    df_scs = pd.read_csv('Jacks_Final_ScS_overlap.rdb',converters=date_time_convert,delim_whitespace=True)
    with open('E_pacific_all_data_binned.sdb','w+') as writer:  
        writer.write('STAT DATE TIME PHASE EVLA EVLO EVDP STLA STLO DIST BAZ AZI LOWMM_LAT LOWMM_LON FAST DFAST TLAG DTLAG Q SNR LOWER_BIN RSIDE_BIN \n')  
        for ph in ['SKS','SKKS']:         
            pplon = '{}_PP_LON'.format(ph) 
            pplat = '{}_PP_LAT'.format(ph) 
            q = 'Q_{}'.format(ph) 
            snr = 'SNR_{}'.format(ph) 
            bn = '{}_BIN'.format(ph) 
            fast = 'FAST_{}'.format(ph) 
            dfast = 'DFAST_{}'.format(ph) 
            tlag = 'TLAG_{}'.format(ph)  
            dtlag = 'DTLAG_{}'.format(ph)  
            for i,row in df_snks.iterrows():
                writer.write('{} {} {} {} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {} {} \n'.format(
                row['STAT'], row['DATE'], row['TIME'], ph,  row['EVLA'], row['EVLO'], row['EVDP'], 
                row['STLA'],row['STLO'], row['DIST'], row['BAZ'], row['AZI'], row[pplat], row[pplon],
                row[fast], row[dfast], row[tlag], row[dtlag], row[q], row[snr], row[bn], row['RSIDE_BIN'])) 
        for i, row in df_scs.iterrows():
            writer.write('{} {} {} {} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {:4.3f} {} {} \n'.format(
            row['STAT'], row['DATE'], row['TIME'], 'ScS',  row['EVLA'], row['EVLO'], row['EVDP'], 
            row['STLA'],row['STLO'], row['DIST'], row['BAZ'], row['AZI'], row['BNCLAT'], row['BNCLON'],
            row['FAST'], row['DFAST'], row['TLAG'], row['DTLAG'], row['Q'], row['SNR'], row['LOWMM_BIN'], row['RSIDE_BIN']))
    print('Done')

def find_phases_in_domain(phasefile,dom_ID,layer='Lower',crit=5.0,save=False):
    '''
    This function searches the input phase file for all paths that go through the given domain (intended use is D'' domains). This subset of phases is then written out for future use.
    '''
    date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
    df = pd.read_csv(phasefile,converters=date_time_convert,delim_whitespace=True)
    domains = np.loadtxt('/Users/ja17375/SWSTomo/Inversions/T3_global.bins',skiprows=1)
#     dom = domains[domains[:,0] == dom_ID]
    dom = domains[np.isin(domains[:,0],dom_ID)]
    dlat = dom[0,1]
    dlon = dom[0,2]
    idx = [ ] # empty list to hold indexs of phases we want
    for i, phase in df.iterrows():
        if layer == 'Lower':
            plat = phase.LOWMM_LAT
            plon = phase.LOWMM_LON
        elif phase == 'RSide':
            plat = phase.STLA
            plon = phase.STLO
        
        dist = vincenty_dist(dlat, dlon, plat, plon)[0]
        
        if dist <= crit:
            #print(phase.STAT, phase.PHASE, dist)
            #Is the phase within a fixed distance from domain center
            idx.append(i)
    df_out = df.iloc[idx]
    if save:
        df_out.to_csv('Domain_{}_SNR10_goodQ_allphases.sdb'.format(dom_ID),index=False,sep=' ')
    else:
        return df_out
    
def split_by_stats(phasefile, stats):
    '''
    function to breakup phase list by station and write out these subsections
    
    Args: 
        phasefile (str) - path to phase file to use
        stats (array-like) - array of stations to iterate over
    
    Returns:
    
        writes out textfiles of phases lists for each station
    '''
    date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
    df = pd.read_csv(phasefile,converters=date_time_convert,delim_whitespace=True)
    
    for stat in stats:
        stat_df = df[df.STAT == stat]
        print(stat, len(stat_df))
        stat_df.to_csv('./StatCorr/{}_goodQ.sdb'.format(stat),sep=' ',index=False)
        
    print('Done')

if __name__ == '__main__':
    Set = PathSetter(phasefile='/Users/ja17375/SWSTomo/BlueCrystal/LDom_1160_phases_SNR10_goodQ.sdb',
    phases=['SKS','SKKS','ScS'], domains='T3_global.bins',
    odir='/Users/ja17375/SWSTomo/BlueCrystal/Dom1160/StackedCorrs/',model='D1160_Stacked_StatCorr_model.xml')
    Set.gen_PathSet_XML(phases=['SKS','SKKS','ScS'],
                        fname='D1160_Stacked_SC_Paths',crit=5.0)
    