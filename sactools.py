from pathlib import Path
from xml.etree import ElementTree 
from shutil import copy
from obspy.signal.rotate import rotate2zne

def get_mts(fileID,stat,phase):
    '''Function to get the .mts file for a phase and read in the xml.

        Args:
            fileID - filename (minus extension) of the data XML we want to read
            stat - the station code (data is stored by station code)
            phase (str) - phase code (ScS,SKS,SKKS, etc.) to fetch data for
        
        Returns:
            data (ElementTree) - ElementTree XML structure containing the required data XML

        Examples:
            >>> Set.get_mts('SKS')

    '''
    if phase == 'ScS':
        path = '/Users/ja17375/SWSTomo/ScS_data'
    else:
        path = '/Users/ja17375/SWSTomo/SnKS_data'
    mts = '{}/{}/{}/{}.mts'.format(path,stat,phase,fileID)
    xml = ElementTree.parse(mts) # Parse the xml file (output from sheba as .mts)
    data = xml.getroot() # Gets the root element of the XML. In this case (the .mts) this is the tag <data> which we want to inject into the
                                 # the bigger Pathset XML file
    for file in ['file1','file2','file3']:# Loop over file tags to add in pointer to data dir
        f = data.find(file).text
        # print(f)
        f_new = 'data/{}'.format(f)
        data.find(file).text = f_new

    return data

def get_sac(fileID,stat,phase):
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
        f = Path('/Users/ja17375/SWSTomo/data/{}.BH{}'.format(fileID,comp))
        if f.is_file():
            # print('/Users/ja17375/SWSTomo/BluePebble/E_pacific/data/{}.BH{} exists, not copying'.format(fileID,comp))
            pass
        else:
            # print('File not found, copying from Sheba Run Dir E_pacific if possible')
            if phase == 'ScS':
                path = '/Users/ja17375/SWSTomo/ScS_data'
            else:
                path = '/Users/ja17375/SWSTomo/SnKS_data'

            file = '{}/{}/{}/{}.BH{}'.format(path,stat,phase,fileID,comp)
            dst = '/Users/ja17375/SWSTomo/data/{}.BH{}'.format(fileID,comp)
            p = copy(file, dst)
            
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
    