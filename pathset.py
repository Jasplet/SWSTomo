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
from xml.etree.ElementTree import Element,SubElement
from xml.dom import minidom
import pandas as pd
import os
import glob
import numpy as np
from obspy.clients import iris

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
    Inputs:
            Required:
            #########
            df [obj]   -  filename for the relevent .pairs file, with T3 bins for SKS, SKKS in lower domains and in the upper mantle
                         Titled: SKS_BIN, SKKS_BIN, UPPER_BIN (Setter will select rows from the relevent station)
            domains    - A domains file (e.g. E_pac.T3.counts.doms). This file must contain the trigonal domain midpoints (having Vertices is
                         nice, as it the Upper and Lower domains counts for each bin (i.e. how many paths were assigned to each bin by geogeom)
            Optional:
            station [str] - the station code for the station[s] we want to include data for (starting with a single station).
                            Now that we are looking at all stations in the East_Pacific we can read them in from a textfile
            odir [str]    - The output directory. Path to where we want our output. If none, we use the current working directory
            model [str]   - Name of the model file (assumed to be within MTS_Setup) that is to be used. If none is provided Model.xml is used
            config_uid [str] - An optional unique identifier that will be added to the MTSConfig file
    """
    def __init__(self,file1,phase1,domains,file2=None,phase2=None,model=None,odir=None,config_uid='Test Run'):
        '''
        Initialise Pathsetter with some essential data, metadata. We have to use file1, file2 to handle that SKS-SKKS and ScS splitting result files contain slightly different information

        Args:
            file1 (str) - file that contains BINNED splitting results 
            phase1 (str) - phase (or phases) contained in file1
            domains (str) - file that contains information on the size/ shape of the domains
            file2 (str) - file that contains BINNED splitting results for a second set of phases
            phase2 (str) - phase (or phases) contained in file2
            model (str) - model XML file.If the model file does not exists Pathsetter can create one
            odir (str) - the output directory that XML files are written to
            config_uid (str) - a unique identifier written to MTSConfig.xml, default is 'Test Run'

        Returns:
            Set (obj) - a configured Pathsetter Object.

        Examples:
            >>> Set = Pathsetter('E_pacific_05_binned.pairs',['SKS','SKKS'],
                                'ScS_w_bouncepoints_binned_05.sdb','ScS','TestModel.xml')                       
        '''
        self.phases = []
        print('Reading df')
        self.date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
        self.df_snks = pd.read_csv(file1,converters=self.date_time_convert,delim_whitespace=True)
        stat1 = self.df_snks[['STAT','STLA','STLO']]
        [self.phases.append(p) for p in phase1] # do list comp to handle a list of phases being specified
        if file2 is not None:
            print('Second data file added, reading')
            self.df_scs = pd.read_csv(file2,converters=self.date_time_convert,delim_whitespace=True)
            stat2 = self.df_scs[['STAT','STLA','STLO']]
            stations = stat1.append(stat2)
            [self.phases.append(p) for p in phase2]
        else:
            self.df_scs = None
        print(self.phases)
        self.stations = stations.drop_duplicates()
        print('Phases specified are: ', self.phases)
        print('Read Trigonal Domians file')
        self.doms = pd.read_csv(domains,delim_whitespace=True)
        print('Set Outdir')
        if odir == None:
            self.opath = os.getcwd() # Gets current working directory and stores is as a Path
            self.odir = self.opath.split('/')[-1]
        elif len(odir.split('/')) > 1:
            self.opath = odir
            self.odir = odir.split('/')[-1]
        else:
            self.opath = '/Users/ja17375/SWSTomo/BluePebble/{}'.format(odir)
        # Set XML namespace
        self.xmlns = {'mtsML':'http://www1.gly.bris.ac.uk/cetsei/xml/MatisseML/'}

        if model == None:
            print('Model will need to be generated before Pathsetting can be done')
            # self.modelxml = '/Users/ja17375/SWSTomo/MTS_Setup/Model.xml'
        else:
            print("Reading Model file from SWSTomo/Models".format(os.getcwd()))
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

        for dom in domains:
#             uid_tmp = dom.find('mtsML:domain_uid',self.xmlns).text
            uid_tmp = dom[0].text
            if uid_tmp == domain:
                print('Domain {} Found'.format(uid_tmp))
                uid = uid_tmp
                if uid.split('_')[0] == 'Lower':
                    # Domains starting with Lower are at CMB. So depth == 2890 km
                    depth = 2890. # Approx depth of CMB [km] 
                
                    slw = get_rayparam(evdp,gcarc,phase)
                    aoi = slw2aoi(depth,slw) # Calculate ray param and then incidence angle
                    dom_h = 250
                elif uid.split('_')[0] == 'Upper':
                    depth = 250. # Depth of upper domain (keeping domain the same thickness for now)
                   
                    slw = get_rayparam(evdp,gcarc,phase)
                    aoi = slw2aoi(depth,slw) # Calculate ray param and then incidence angle
                    dom_h = 250
                elif uid.split('_')[0] == 'SSide':
                    # If domain ID is SSide then we know this domain is for a source side correction. 
                    azi = 0 # For the splitting corrections we fix azi to 0
                    aoi = 0
                    dom_h = 100
                    print('SSide')
                else:
                    raise ValueError(f'Domain \"{uid_tmp}\" not found ')
                
        dist = dom_h / np.cos(np.deg2rad(aoi))
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

        '''
        dml = '{{{}}}domain'.format(self.xmlns['mtsML']) # need the triple curly braces to get a set of braces surrounding the mtsML
        u = []
        l = []
        sdoms = []
        self.udom_c = {}
        self.ldom_c = {}
        self.sdom_c = {}
        print(dml)
        for dom in self.model.iter(dml):
            d_id = dom[0].text

            if d_id.split('_')[0] == 'Upper':
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
        udoms = self.doms[self.doms['BIN'].isin([i.split('_')[1] for i in u ])]        
        return (udoms,ldoms,sdoms)

    def is_phase_good(self,row,ph):
        '''Function to test if a phase is not split or null based on its Q '''
        Q = 'Q_{}'.format(ph)
        if (row[Q] < 0.5) and (row[Q] > -0.7):
            # print('Phase {} fails Q tests, continuing to next'.format(ph))
            sn = False
        else:
            sn = True
        return sn

    def loop_thru_domains(self,doms,raylat,raylon,uid,ph,crit,evdp,gcarc,azi):
        ''' 
        This function loops through a set of input domains and identifies domains whose midpoints are within a set distance from a raypath (the point the raypath is at at a given depth in the mantle)
        
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
                                                
    def gen_PathSet_XML(self,phases=['ScS','SKS','SKKS'],fname=None):
        '''
        Function to iterate over the DataFrame and construct the XML we need
        Phases [list] - list of phase codes to iterate over (assuming each row in the DataFame corresponds to all phases)
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
        dom_used = []
        ##
        p = 1
        nf = 0 # counter for files not found
        q_fail = 0 # counter for pahses that fail Q test
        ex = 0
        # Before Pathsetting, parse the upper/lower domains from the model file
        rdoms,ldoms,sdoms = self.parsedoms()
        crit = 5.5 # [deg] - the distance criterea for including phases in a domain.
                   #         designed to give some overlap in neighbouring domians for reduce edge effects...
        ## Wrangle the SnKS and ScS dataframes to produce one combined dataframe we can then iterate over
        rows2comp = ['DATE','TIME','STAT','EVLA','EVLO','EVDP','STLA','STLO','AZI','BAZ','DIST']
        ScSrows = ['Q','SNR','BNCLAT','BNCLON']
        SnKSrows = ['Q_SKS','Q_SKKS','SKS_PP_LAT','SKS_PP_LON','SKKS_PP_LAT','SKKS_PP_LON']
        ScS_t = self.df_scs[rows2comp + ScSrows]
        SnKS_t = self.df_snks[rows2comp + SnKSrows]
        pdf = pd.concat([ScS_t,SnKS_t],sort=False)
        pdf.rename(columns={'Q':'Q_ScS','SNR':'SNR_ScS','BNCLAT':'ScS_LM_LAT','BNCLON':'ScS_LM_LON',
                            'SKS_PP_LAT':'SKS_LM_LAT','SKS_PP_LON':'SKS_LM_LON',
                            'SKKS_PP_LAT':'SKKS_LM_LAT','SKKS_PP_LON':'SKKS_LM_LON'},inplace=True)
        for i, row in pdf.iterrows():
            # All XML generation must sit within this loop (function calls) so that we make sure that Az, EVDP etc. are right for the current phase
            attribs = {'evdp':row.EVDP,'azi':row.AZI,'gcarc':row.DIST,'evla':row.EVLA,
                            'evlo':row.EVLO,'stla':row.STLA,'stlo':row.STLO,'stat':row.STAT,
                            'date':row.DATE,'time':row.TIME,'stat':row.STAT
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
                    continue
                if ph_good is True:                  
                    llat = row['{}_LM_LAT'.format(ph)]
                    llon = row['{}_LM_LON'.format(ph)]
                    lowmm_ops = self.loop_thru_domains(ldoms,llat,llon,"Lower",ph,crit,
                                                       attribs['evdp'],attribs['gcarc'],attribs['azi'])
                    rside_ops = self.loop_thru_domains(rdoms,attribs['stla'],attribs['stlo'],"Upper",
                                                       ph,crit,attribs['evdp'],attribs['gcarc'],attribs['azi'])
                    sside_ops = [0] # This is so the Pathsetter will iterate at least once (so for SnKS phases it will still append the right domains)
                    if ph == 'ScS': # Only ScS needs source side domains
                        if len(sdoms) > 0:
#                            We have source side correction domians for each ScS phase. 
#                            Now we need to find the right correction for the current phase in the loop.
#                            As we have precalculated the corrections, a domain should exist for each ScS phase. 
#                            These corrections domains are identified by Station Code and event Date/Time
                                dom_id = 'SSide_{}_{}_{}'.format(attribs['stat'],attribs['date'],attribs['time'])
                                
                                sside_ops[0] = self.domain2operator(dom_id,'ScS',attribs['evdp'],attribs['gcarc'],
                                                     attribs['azi'])
                        else:
#                           Overwrite sside_ops with a list of operators from loop_thru_domains
                            sside_ops = self.loop_thru_domains(udoms,attribs['evla'],attribs['evlo'], "Upper",
                                                               ph,crit,attribs['evdp'],
                                                               attribs['gcarc'],attribs['azi'])                   
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
                                    path.append(s_op)

                                path.append(l_op) 
                                path.append(r_op)

                else:
                    q_fail +=1
        #Now write out the Pathset XML
        self._write_pretty_xml(self.pathset_root,file='{}/{}.xml'.format(self.opath,self.pathset_xml))
        # self._write_domain_counts(file='{}/E_pac.T3.goodphases.domains.counts'.format(self.opath))

    def baz_test(self,baz):
        '''
        Function to test which 30 deg backzimuth bin each phase sits in
        '''
        for i,v in enumerate([30,60,90,120,150,180,210,240,270,300,330,360]):
            if (baz > v-30) and (baz <= (v)):
                return "{:03d}".format(v)




    def gen_Model_XML(self,mod_name=None,Low_Domains=None,Rside_Domains=None,Sside_Domains=None):
        '''
        Generate a default Model XML file, containing all requested Lower Domains (default is all in .domains file)
        and all corresponding Upper domains.
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
            lowb = self.df_snks.SKS_BIN.append(self.df_snks.SKKS_BIN)
            if self.df_scs is not None:
                lowb2 = lowb.append(self.df_scs.LOWMM_BIN)
            Low_Domains = lowb2.drop_duplicates().values

        if Rside_Domains is None:
            print('Finding all Reciever Side Domains')
            rside_sks = self.df_snks[self.df_snks['SKS_BIN'].isin(Low_Domains)]['RSIDE_BIN']
            rside_skks = self.df_snks[self.df_snks['SKKS_BIN'].isin(Low_Domains)]['RSIDE_BIN']
            rdf = rside_sks.append(rside_skks)
            if self.df_scs is not None:
                rside_scs = self.df_scs[self.df_scs['LOWMM_BIN'].isin(Low_Domains)]['RSIDE_BIN']
                rdf = rdf.append(rside_scs)
            rside = rdf.drop_duplicates().values
        else:
            rside = Rside_Domains

        if Sside_Domains:
#             If no Sside Domains requsted, we assume each ScS phase will have its own fixed "sourceside 
#             domain" with a correction. 
            sside = Sside_Domains
        else:
            sside = None
            for i,row in self.df_scs.iterrows():
                date, time, stat = row.DATE, row.TIME, row.STAT
                sdom = add_sside_correction(date, time, stat)
                m.append(sdom)
                s += 1
                
        if sside:
            ud = np.concatenate((rside,sside))
            udoms = np.unique(ud)
        else:
            udoms = rside
        
        for udom in udoms:
            dom = bin2domain('Upper_{}'.format(int(udom)))
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

    def _write_domain_counts(self,file):
        """
        Adds the domains phase counts to the domains DF and then writes it out
        """
        self.doms['Lower_Count'] = self.doms['BIN'].map(self.ldom_c)
        self.doms['Upper_Count'] = self.doms['BIN'].map(self.udom_c)
        self.doms = self.doms.fillna(0)
        self.doms.to_csv(file,sep=' ',index=False)
