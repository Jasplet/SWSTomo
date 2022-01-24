from pathlib import Path
from xml.etree import ElementTree 
from shutil import copy
from obspy.signal.rotate import rotate2zne

def get_mts(mts_file,stat=None,phase=None, syn=False):
    '''Function to get the .mts file for a phase and read in the xml.

        Args:
            mts_file- filename (minus extension) of the data XML we want to read
            stat - the station code (data is stored by station code)
            phase (str) - phase code (ScS,SKS,SKKS, etc.) to fetch data for
        
        Returns:
            data (ElementTree) - ElementTree XML structure containing the required data XML

        Examples:
            >>> Set.get_mts('SKS')

    '''
    if syn:
        mts_file = '{}.mts'.format(mts_file)
        
    xml = ElementTree.parse(mts_file) # Parse the xml file (output from sheba as .mts)
    data = xml.getroot() # Gets the root element of the XML. In this case (the .mts) this is the tag <data> which we want to inject into the
                                 # the bigger Pathset XML file
    for file in ['file1','file2','file3']:# Loop over file tags to add in pointer to data dir
        f = data.find(file).text
        # print(f)
        f_new = 'data/{}'.format(f)
        data.find(file).text = f_new

    return data

def get_sac(fileID,dst_path,src_path,stat,phase, syn=False):
    '''Function to copy data to the local data directory "data" if it is not already there (mainly useful for test/first runs).

       Args:
            fileID - filename (minus extension) of the sac data to copy
            stat - the station code (data is stored by station code)
            phase (str) - phase code (ScS,SKS,SKKS, etc.) 
       Returns:

       Examples:
           >>> Set.get_sac('SKS')
    '''
    for comp in ['E','N','Z']:
        f = Path('/Users/ja17375/Projects/Epac_fast_anom/data/{}.BH{}'.format(fileID,comp))
        if f.is_file():
            #print(f'{f} exists, not copying'.format(fileID,comp))
            pass
        else:
            print(f'File not found, copying from {src_path} if possible')     
            if syn:
                file = f'{fileID}.BH{comp}'
                f = fileID.split('/')[-1]
                dst = '/Users/ja17375/Projects/Epac_fast_anom/data/{}.BH{}'.format(f,comp)
                print(f'Copy synthetic data {file}')
            else:
                dst = f'{dst_path}/data/{fileID}.BH{comp}'
                file = f'{src_path}/{stat}/{phase}/{fileID}.BH{comp}'
                    
            _ = copy(file, dst)
        
def rotate_traces(st):
    '''
    Function to rotate an obspy stream (assuming a SAC file read in) to ZNE using rotate2zne
    '''
    bh1 = st.select(channel = 'BH1')
    bh2 = st.select(channel = 'BH2')
    bhz = st.select(channel = 'BHZ')
    
    bh1_data = bh1[0].data
    bh1_inc = bh1[0].stats.sac['cmpinc']
    bh1_az = bh1[0].stats.sac['cmpaz']
    
    bh2_data = bh2[0].data
    bh2_inc = bh2[0].stats.sac['cmpinc']
    bh2_az = bh2[0].stats.sac['cmpaz']
    
    bhz_data = bhz[0].data
    bhz_inc = bhz[0].stats.sac['cmpinc']
    bhz_az = bhz[0].stats.sac['cmpaz']
    
    (Z,N,E) = rotate2zne(bhz_data, bhz_az, bhz_inc,
                         bh1_data, bh1_az, bh1_inc,
                         bh2_data, bh2_az, bh2_inc,
                         )
    
    st.select(channel = 'BH1')[0].data = E
    st.select(channel = 'BH2')[0].data = N
    st.select(channel = 'BHZ')[0].data = Z
    
    return st 
    