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
from corrections import add_sside_correction, add_station_correction,bin2domain
###############################
__author__ = "Joseph Asplet"
__license__ = "MIT"
__version__ = "0.1"
__email__ = "joseph.asplet@bristol.ac.uk"
__status__ = "Development"

# Set XML namespace as global varaible
xmlns = {'mtsML':'http://www1.gly.bris.ac.uk/cetsei/xml/MatisseML/'}
DATADIR = '/Users/ja17375/Projects/Splitting_PhD/Runs/E_pacific/'

class PathSetter:
    """A class to hold the metadata for the run (rdir, station? [for now], outdir etc. ) and fucntions
    to parse/generate the XML needed for the PathSet file
    """
    def __init__(self,phasefile, model,odir=None,config_uid='Test Run'):
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
        print('Read Trigonal Domians file')
        print('Set Outdir')
        if odir == None:
            self.opath = os.getcwd() # Gets current working directory and stores is as a Path
            self.odir = self.opath.split('/')[-1]
        elif len(odir.split('/')) > 1:
            self.opath = odir
            self.odir = odir.split('/')[-1]
        else:
            raise ValueError('Output directory/path required')      
        print(f"Reading Model file {model}")
        self.modelxml = '{}'.format(model)
        self.modelroot = ElementTree.parse('{}'.format(self.modelxml)).getroot()
        self.model = self.modelroot.find('mtsML:model',xmlns)
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
        domains = self.model.findall('mtsML:domain',xmlns)
        uid = None
        for dom in domains:
            # Loop through domains in the model file 
            uid_tmp = dom[0].text
            if domain == uid_tmp:
                uid = uid_tmp
                break
        
        if uid is None:
            raise ValueError('Domain {} Not Found'.format(domain))
        
        if phase not in ['SKS', 'SKKS']:
            if (phase =='ScS') and (uid.split('_')[0] == 'SSide'):
                #Bodge for now. Fix when refectoaring domain2operator.
                pass
            else:
                raise ValueError(f'Unsupported phase > {phase}')
        
        if uid.split('_')[0] == 'Lower':
            # Domains starting with Lower are at CMB. So depth == 2890 km. 
            # We assume SKS, SKKS piercepoints are at this depth.
            depth = 2890. # Approx depth of CMB [km] 
            slw = get_rayparam(evdp,gcarc,phase)
            aoi = slw2aoi(depth,slw) # Calculate ray param and then incidence angle
            dom_h = 250
            dist = dom_h / np.cos(np.deg2rad(aoi))
        elif (uid.split('_')[0] in ['Upper','RSide','Station','SSide']):
            # These prefixes should be used for correction domians.
            aoi = 0 
            dist = 100 # 100km so we can use conversion dt/100 = strength
        else:
            raise ValueError(f'Domain "{uid}"name does not match expected values')
            
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
    
    def configure_stat_corr_operator(self, row):
        station_corr = ElementTree.Element('operator')
        dom_uid = ElementTree.SubElement(station_corr,'domain_uid')
        dom_uid.text = f'Station_{row.STAT}'
        azimuth = ElementTree.SubElement(station_corr,'azi')
        azimuth.text = '0'
        inclination = ElementTree.SubElement(station_corr,'inc')
        inclination.text = "90" # Inclination is 90 for corrections
        length = ElementTree.SubElement(station_corr,'dist')
        length.text = "100"
        return station_corr
    
    def configure_scs_operators(self, row, direction, dom_id):
        '''Handles the (slightly more involved) process of setting up ScS. 
        Downgoing leg should be applied first with a negative aoi. This is because
        a positive aoi corresponds to a ray propagating up through the domain.
        Input:
            row (Series) - row of a DataFrame
        '''
        # Calculate the ray params for ScS
        # Now there is a bit of user choice here in what angle to use. 
        # I am going to to chose to use the incindence angle at the CMB bouncepoint.
        depth = 2890. # Approx depth of CMB [km] 
        slw = get_rayparam(row.EVDP,row.GCARC,phase='ScS')
        aoi = slw2aoi(depth,slw) # Calculate ray param and then incidence angle
        dom_h = 250
        dist = dom_h / np.cos(np.deg2rad(aoi))
        _ , lmm_azi = vincenty_dist(row.LOWMM_LAT, row.LOWMM_LON, row.STLA, row.STLO)
        #  Downgoing ScS Limb. 
        scs_op = ElementTree.Element('operator')
        dom_uid = ElementTree.SubElement(scs_op,'domain_uid')
        dom_uid.text = dom_id
        azimuth = ElementTree.SubElement(scs_op,'azi')
        azimuth.text = str(lmm_azi)
        inclination = ElementTree.SubElement(scs_op,'inc')
        if direction == 'Upgoing':
            inc = (90 - aoi)
        elif direction == 'Downgoing':
            inc = (90 - aoi) * -1
        inclination.text = str(inc) # Inclination is 90 for corrections
        length = ElementTree.SubElement(scs_op,'dist')
        dist = dom_h / np.cos(np.deg2rad(aoi))
        length.text = str(dist)
        
        return scs_op
        
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
        dml = '{{{}}}domain'.format(xmlns['mtsML']) 
        # need the triple curly braces to get a set of braces surrounding the mtsML
        rdoms = []
        ldoms = []
        sdoms = []
        self.udom_c = {}
        self.ldom_c = {}
        self.sdom_c = {}
        for dom in self.model.iter(dml):
            d_id = dom[0].text
            if d_id.split('_')[0] == 'Station':
                rdoms.append(d_id)
                self.udom_c.update({d_id : 0})
            elif d_id.split('_')[0] == 'Lower':
                ldoms.append(d_id)
                self.ldom_c.update({d_id : 0})
            elif d_id.split('_')[0] == 'SSide':
                # If the uid is SSide we know it is a source side correction domain. We just need to 
                # match the correct domain to the correct phase 
                sdoms.append(d_id)
                self.sdom_c.update({d_id : 0 })
                      
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
 
    def gen_PathSet_XML(self,stations, lower_dom_no=None, phases=['ScS','SKS','SKKS'],fname=None, f=1):
        '''
        Function to iterate over the shear-wave phase data provided to Setter and to generate the correct XML
        to from the Pathset XML file required for Matisse. This includes adding source/ reciever side corrections
        where required (and where the neccesary corrections have been provided as needed by functions called
                        by gen_PathSet_XML)
        This funciton is written assuming there is a singular D'' domain to invert for 
        
        Args:
            stations [list] - list of stations to include data for
            lower_dom_no [int or list] - lower domain of interest
            phases [list] - list of phase codes to iterate over. Default is ['ScS','SKS','SKKS']
            fname [str] - File name to write the Pathset XML to. If no name is provided then this defaults to 'Pathset'
        Returns:
            Output XML is written to fname
        '''
        # Some quick i/o management
        if fname == None:
            self.pathset_xml = 'Pathset' 
        else:
            self.pathset_xml = fname
        
        #Start of main function
        self.pathset_root = ElementTree.Element('MatisseML')
        self.pathset_root.set("xmlns",xmlns['mtsML'])
        pathset = ElementTree.SubElement(self.pathset_root,'pathset')
        psuid = 'Paths for run in dir {} .'.format(self.odir)
        pathset_uid = ElementTree.SubElement(pathset,'pathset_uid')
        pathset_uid.text = psuid
        p = 0
        # Before Pathsetting, parse the upper/lower domains from the model file
        rdoms,ldoms,sdoms = self.parsedoms()
        if stations == 'all':
            df = self.df
        elif stations == 'random_sample':
            df = self.df.sample(frac=f, replace=True)
        else:   
            df = self.df[self.df.STAT.isin(stations)]
        
        for i, row in df.iterrows():
            # All XML generation must sit within this loop (function calls) so that we make sure that Az, EVDP etc. are right for the current phase
            attribs = {'evdp':row.EVDP,'azi':row.AZI,'gcarc':row.GCARC,'evla':row.EVLA,
                            'evlo':row.EVLO,'stla':row.STLA,'stlo':row.STLO,
                            'date':row.DATE,'time':row.TIME,'stat':row.STAT,
                            'phase' :row.PHASE}
            filepath = f'{DATADIR}/{row.STAT}/{row.PHASE}/{row.STAT}_{row.DATE}_{row.TIME}??_{row.PHASE}.mts'
            try:
                mts_file = glob.glob(filepath)[0]
                fileID = mts_file.strip('.mts').split('/')[-1]
                # Strip out .mts and split by '/', select end to get filestem
            except IndexError:
                print(f'WARNING: Glob has failed to find anything for {filepath}')
                continue

            ldom_id = 'Lower_{}'.format(lower_dom_no)
            get_sac(fileID, self.opath, DATADIR, attribs['stat'],attribs['phase'])
            # Now make XML for this Path
            path = ElementTree.SubElement(pathset,'path')
            pathname = 'Path {} {}'.format(p,attribs['phase'])          
            path_uid = ElementTree.SubElement(path,'path_uid')
            path_uid.text = pathname
            # Add Data (from .mts)
            data = get_mts(mts_file,attribs['stat'],attribs['phase'])
            path.append(data)
            stat_uid = ElementTree.SubElement(path,'station_uid')
            stat_uid.text = attribs['stat']
            evt_uid = ElementTree.SubElement(path,'event_uid')
            evt_uid.text = '{}_{}'.format(attribs['date'],attribs['time'])
            if attribs['phase'] == 'ScS':
            # Only ScS needs source side domains
                dom_id = 'SSide_{}_{}_{}'.format(attribs['stat'],
                                                 attribs['date'],attribs['time'])
                sside_op = self.domain2operator(dom_id,'ScS',attribs['evdp'],
                                                attribs['gcarc'],0)    
                # Add SSide Operator to path
                scs_dwn_op = self.configure_scs_operators(row, 'Downgoing', ldom_id)
                scs_up_op = self.configure_scs_operators(row, 'Upgoing', ldom_id)
                stat_op = self.configure_stat_corr_operator(row)
                path.append(sside_op)
                path.append(scs_dwn_op)
                path.append(scs_up_op)
                path.append(stat_op)
            elif attribs['phase'] in ['SKS', 'SKKS']:
                #Use the old way for now, maybe replace later?
                _ , lmm_azi = vincenty_dist(row.LOWMM_LAT, row.LOWMM_LON, row.STLA, row.STLO)
                lowmm_op = self.domain2operator(ldom_id, attribs['phase'], attribs['evdp'],
                                                attribs['gcarc'], lmm_azi)
                # Add D`` Operator
                path.append(lowmm_op)
                stat_op = self.configure_stat_corr_operator(row)
                path.append(stat_op)
            else:
                raise(ValueError('Unknown Phase'))
            p+=1
        #Now write out the Pathset XML
        tree = ElementTree.ElementTree(self.pathset_root)
        ElementTree.indent(tree, space="\t", level=0)
        pathset_out = f'{self.opath}/{self.pathset_xml}.xml'
        print(f'Pathset written to {pathset_out}')
        tree.write(pathset_out, encoding="utf-8")

    def gen_Model_XML(self, stations, mod_name=None, Low_Domains=None, ):
        '''
        Generates the model XML file required by matisse, with lowermost mantle domains and seperate upper mantle domains for every station requested. Source Side corrections for ScS are auto-added
        
        Args:
            stations (list) - list of all Stations to create domains for. 
            mod_name (str) - name of the model. The model XML will be written to a file mod_name.xml
            Low_Domains (array-like) - [optional] array of D`` domains to use in model. If none we assume this is a single domain run 
        Returns:
            model XML written to the file mod_name.xml
        '''
        if mod_name is None:
            mod_name = input('No model name provided. Enter one now :')
        modpath = '/Users/ja17375/Projects/Epac_fast_anom/Models'
        # Write the initial xml
        root = ElementTree.Element('MatisseML')
        tree = ElementTree.ElementTree(root)
        root.set("xmlns",xmlns['mtsML'])
        root.append(ElementTree.Comment('Generated by gen_Model from pathset.py'))
        m = ElementTree.SubElement(root,'model')
        m_uid = ElementTree.SubElement(m,'model_uid')
        m_uid.text = 'Epac_fast_anomaly Inversion'
        
                        
        if Low_Domains:
            for ldom in Low_Domains:
                # Loop over requested lower domains
                print(ldom)
                dom = bin2domain(f'Lower_{ldom}')
                m.append(dom)
        
        if stations == 'all':
            stations =  self.df.STAT.unique()
        
        for stat in stations:
            dom = add_station_correction(stat)
            m.append(dom)
            
        for i,row in self.df.iterrows():
            if (row.STAT in stations) & (row.PHASE == 'ScS'):
                date, time, stat = row.DATE, row.TIME, row.STAT
                print(f'{date}, {time}, {stat}')
                sdom = add_sside_correction(date, time, stat)
                m.append(sdom)

        else:
            print('No Lowermost Mantle Domains')
            
    
        print('Model Generated')
        self.model = m
        _write_pretty_xml(root,file='{}/{}.xml'.format(modpath,mod_name))        
            

 
def gen_MTS_Info(opath, config_uid):
    '''
    Function to create the simple XML file "MTS_Info"
    '''
    root = ElementTree.Element('MatisseML')
    tree = ElementTree.ElementTree(root)
    root.set("xmlns",xmlns['mtsML'])
    root.append(ElementTree.Comment('Generated by gen_MTS_Info from pathset.py'))
    samplerinfo = ElementTree.SubElement(root,'SamplerInfo')
    c_uid = ElementTree.SubElement(samplerinfo,'config_uid')
    c_uid.text = config_uid
    # Now write the XML to a file
    write_pretty_xml(root,file='{}/MTS_Info.xml'.format(opath))

def gen_MTS_Config(config_uid, pathset_xml, opath, options=None):
    '''
    Function to create the MTS config file. N.B this function will just generate a default config file. It is easier to edit/copy an existing file as only a few changes are needed to tweak a run.
    
    Args:
        options (dict) - options to be set in the Config file. N.B all options must be specified

    '''
    model = modelxml.split('/')[0]

    root = ElementTree.Element('MatisseML')
    tree = ElementTree.ElementTree(root)
    root.set("xmlns",xmlns['mtsML'])
    root.append(ElementTree.Comment('Generated by gen_MTS_Config from pathset.py'))
    config = ElementTree.SubElement(root,'config')
    c_uid = ElementTree.SubElement(config,'config_uid')
    c_uid.text = config_uid
    config.append(ElementTree.Comment(' input files'))
    mf = ElementTree.SubElement(config,'model_file') # Creat model file tag
    mf.text = model
    ps = ElementTree.SubElement(config,'pathset_file')
    ps.text = '{}.xml'.format(pathset_xml)
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
    _write_pretty_xml(root,file='{}/MTSConfig.xml'.format(opath))


def _write_pretty_xml(root,file):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(root)
    reparsed = minidom.parseString(rough_string)
    print(reparsed)
    with open('{}'.format(file),'w') as writer :
        reparsed.writexml(writer,indent="    ",addindent="    ",newl="\n",encoding="utf-8")

        print('XML written to {}'.format(file))

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
    Set = PathSetter(phasefile='/Users/ja17375/Projects/Epac_fast_anom/HQ_data/HQ_phases_on_fast_anom.sdb',
                      odir='/Users/ja17375/Projects/Epac_fast_anom/HQ_data',
                      model='/Users/ja17375/Projects/Epac_fast_anom/Models/EP_fast_anom_Model.xml')
    # Make Pathsets for Epac Data
    Set.gen_PathSet_XML(stations='all',lower_dom_no='fast_anom',
                        phases=['ScS','SKS','SKKS'],fname='EP_fast_anom_Paths_ScS_fix')
# try expanding the dataset
    # Set = PathSetter(phasefile='/Users/ja17375/Projects/Epac_fast_anom/Expanded_Dataset/Full_Send.sdb',
    #                   odir='/Users/ja17375/Projects/Epac_fast_anom/Expanded_Dataset/',
    #                   model='/Users/ja17375/Projects/Epac_fast_anom/Models/EP_full_send.xml')
    # Set.gen_Model_XML(stations='all', mod_name='EP_full_send', Low_Domains=['fast_anom'])
    #Set.gen_PathSet_XML(stations='all',lower_dom_no='fast_anom',
                 #       phases=['ScS','SKS','SKKS'],fname='EP_full_send_pathset')
    