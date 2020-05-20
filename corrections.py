import pandas as pd
from xml.etree import ElementTree 
from scipy.stats import circmean
from numpy import mean,deg2rad,rad2deg,pi,sin,cos,arctan2
import numpy as np 
def bin2domain(dom_name,ac=None,bc=None,gc=None,sc=None):
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

    if ac :
        alpha = ElementTree.SubElement(domain, 'alpha', type="fixed",value=str(ac))
    else:
        alpha = ElementTree.SubElement(domain, 'alpha', type="fixed",value="90")

    if bc :
        beta = ElementTree.SubElement(domain, 'beta', type="fixed",value=str(bc))
    else:
        beta = ElementTree.SubElement(domain, 'beta', type="fixed",value="0")

    if gc :
        gamma = ElementTree.SubElement(domain, 'gamma', type="fixed", value=str(gc))
    else:
        gamma = ElementTree.SubElement(domain, 'gamma', type="periodic",min="-90",max="90",init="0")

    if sc :
        s = ElementTree.SubElement(domain, 'strength', type="fixed", value=str(sc))
    else:
        s = ElementTree.SubElement(domain, 'strength', type="linear",min="0.00",max="0.05",init="0.0125")

    return domain

def add_sside_correction(date,time,stat):
    '''
    This function adds a source side correction for ScS phases. It looks up ScS corrections from a pre-prepared corrections database
    For source side (ScS) domains it is assumed that each ScS phase has its own source side correction.
    By fixing a 100km path length (in the Pathset) we can approximate that dt = 100*strength (derived from 1 domain synthetic example)
    
    Args:
    
        date (str) - event date in julian day format (yyyyjjj)
        time (str) - event time (hhmm)
        stat (str) - station code
    Returns:
    
        corr_dom (ElementTree) - the XML needed for create the domain with the ScS correction added
        
    Example:
        >>>add_sside_correction('2002321','0453','BRK',1)
            <Element 'domain' at 0x11c3c3310>
    '''
    date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
    cdf = pd.read_csv('/Users/ja17375/SWSTomo/SourceSideCorrs/Jacks_ScS_FINAL_w_SSIDE_corr.rdb',
                    delim_whitespace=True,converters=date_time_convert)
    corr = cdf[(cdf.DATE == date) & (cdf.TIME == time) & (cdf.STAT == stat)]
    if len(corr) == 0:
        print('No Correction for this ScS Phase {} {}_{}'.format(stat, date, time))
    elif len(corr) > 1:
        print(corr)
        
    uid = 'SSide_{}_{}_{}'.format(stat,date,time)
    gamma_rad = circmean(deg2rad(corr.S_FAST.values), low = -pi*0.5, high = pi*0.5)
    gamma = rad2deg(gamma_rad)
    strength = mean(corr.S_TLAG.values) / 100
    corr_dom = bin2domain(uid, gc = gamma, sc = strength)
    
    return corr_dom
    
def check_if_sside_corr(df):
    '''
    This function checks a dataframe containing ScS splitting data and returns a dataframe for which there are ScS corrections
    
    Args:
        df
    
    Returns:
        df_w_corr
    '''
    date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
    cdf = pd.read_csv('/Users/ja17375/SWSTomo/SourceSideCorrs/Jacks_ScS_FINAL_w_SSIDE_corr.rdb',
                    delim_whitespace=True,converters=date_time_convert)
    
    m = pd.merge(df,cdf,how='inner',on=['DATE', 'TIME', 'STAT'])
    labels = {'STLA_x':'STLA','STLO_x':'STLO','EVDP_x':'EVDP','EVLA_x':'EVLA','EVLO_x':'EVLO','DIST_x':'DIST',
              'AZI_x':'AZI','BAZ_x':'BAZ','FAST_x':'FAST','DFAST_x':'DFAST','TLAG_x':'TLAG',
              'EIGCORR_x':'EIGCORR', 'Q_x':'Q', 'SNR_x':'SNR', 'NDF_x':'NDF', 'TKOFFANG_x':'TKOFFANG',
              'SRCPOL_x':'SRCPOL', 'RCVPOL_x':'RCVPOL', 'ENTLAT_x':'ENTLAT', 'ENTLON_x':'ENTLON', 
              'BNCLAT_x':'BNCLAT', 'BNCLON_x':'BNCLON','EXTLAT_x':'EXTLAT', 'EXTLON_x':'EXTLON'  }
    m.rename(labels,axis='columns',inplace=True)
    drp = ['EVLA_y', 'EVLO_y', 'STLA_y', 'STLO_y','EVDP_y', 'DIST_y', 'AZI_y', 'BAZ_y', 'FAST_y',
           'DFAST_y', 'TLAG_y','DTLAG_y', 'SPOL_y', 'DSPOL_y', 'WBEG_y', 'WEND_y', 'EIGORIG_y',
           'EIGCORR_y', 'Q_y', 'SNR_y', 'NDF_y', 'FOLDER', 'TKOFFANG_y', 'SRCPOL_y', 'RCVPOL_y',
           'ENTLAT_y', 'ENTLON_y', 'BNCLAT_y', 'BNCLON_y', 'EXTLAT_y', 'EXTLON_y', 'NI', 
           'R_FAST', 'R_DFAST', 'R_TLAG', 'R_DTLAG','R_QC','WQ','GROUP']
    df_w_corr = m.drop(drp,axis='columns')
    
    return df_w_corr

def weighted_circmean(angles,weights=[]):
    '''
    A function to take a weighted arithmetic mean of a set of angles that falls within a specified range. Default is -90 to 90 degrees.
    
    Args:
        angles (array): input array. angles should be in degrees (they will be converted to radians inside the function)
    
        weights (array): the weighting to be applied to the angles. default is to set weights to 1 (i.e take an unweighted mean)
     
        high (float or int): high bound of range that angles should fall in
        
        low (float or int): low bound of range that angles should fall in 
        
    Returns:
        weighted mean 
    
    Examples:
        >>>weighted_circmean(deg2rad([80,-80]),[1,1])
            -1.5707963267948966
        
        >>>weighted_circmean([0.1,0.2,0.3],[2,1,2])
            0.20000000000000018
            
        >>>weighted_circmean([0.5,0,-np.pi],[2,1,2])
            0.19558559642933782
        
        >>>weighted_circmean([0.5,0,0],[2,1,2])
            0.19558559642933782
    '''
    high=pi/2
    low=-pi/2
    angles = np.asarray(angles)
    weights = np.asarray(weights)
#   Check there are enough weights for angles 
    if weights.size == 0:
        weights = np.ones(angles.size)
    elif weights.size != angles.size:
        raise ValueError('The number of elements in angles ({}) and weights ({}) do not match'.format(angles.size,weights.size))
    # Translate samples to a range of 0 - 2pi
    samps = (angles - low)*2.*pi / (high - low)
    print(rad2deg(angles))
    print(rad2deg(samps))
    s_angles = weights * sin(samps) / weights.sum()
    c_angles = weights * cos(samps) / weights.sum()
    s_sum = s_angles.sum()
    c_sum = c_angles.sum()
    tmp = arctan2(s_sum,c_sum)
    # transform result back to a range of -pi/2 - pi/2
    res = tmp*(high - low)/2/pi + low
    if res < -(pi/2):
        print(rad2deg(res))
        res = pi/2 - abs(pi/2 - abs(res))
    elif res > (pi/2):
        print(rad2deg(res))
        res = -pi/2 + abs(res - pi/2)
    
    return res

def depth_stack_point(mod_point,depth_max=250):
    '''
    This function takes a single gird (lat,lon) point in a surface wave model and returns a weighted average of the modelled fast direction across all the depth slices
    
    Args:
        mod_point (obj) - numpy array containing [depth,lat,lon,phi]
        depth_max (int) - the maximum depth to average up to. E.g. the default is to only average over the upper 250 km of the mantle. 

    '''
    n = (mod_point[:,0] <= depth_max).sum()# number of depth slices to average
    print(n)
    if mod_point[n,0] != depth_max:
        weights = np.zeros(n+1)
        phis = np.zeros(n+1)
        strengths = np.zeros(n+1)
        # if the last point is not the same as the specified depth max we will fix it to be so (this is a little iffy)
        print(weights)
        weights[n] = (depth_max - mod_point[n,0]) # normalized weighting
    else:
        weights = np.zeros(n)
        phis = np.zeros(n)
        strengths = np.zeros(n)
    
    top = 0 
    for ds in range(0,n):
        phis[ds] = mod_point[ds,3]
        weights[ds] = (mod_point[ds,0] - top) # normalized weighting
        strengths[ds] = mod_point[ds,6]
        top = mod_point[ds,0]
    
    cm = weighted_circmean(np.deg2rad(phis),np.deg2rad(weights))
    str_m = np.sum(weights * strengths) / np.sum(weights)
    
    return rad2deg(cm),str_m

def depth_stack_model(model_file,depth_max=250):
    '''
    This function takes one of the Schaffer models and creates a depth averaged version up to a set depth
    '''
    model = np.loadtxt('Schaeffer_Models/SL2016svA_n-k.mod')
    dm = model[:,0],min()
    n = (model[model[:,0]] == dm).shape[0]
    stacked_model = np.zeros((n,4))
    i = 0
    for lon in model[:,1]:
        for lat in model[:,2]:
            mod_point = model[(model[:,1] == lon) & (model[:,2] == lat)]
            phi, strength = depth_stack_point(mod_point,depth_max)
            stacked_model[i,0] = lon
            stacked_model[i,1] = lat
            stacked_model[i,2] = phi
            stacked_model[i,3] = strength
    
    return stacked_model

def resample_model():
    '''
    This function takes the depth averaged model and interpolates it onto our T3 domain grid (triangular bin midpoints)
    To interpolate the phi values shift the range to 0 - 180 and use the %180 
    '''

def add_rside_correction(dom,type):
    '''
    This function looks up a domain correction for the input upper mantle domain (domain IDs assigned by geogeom)
    
    For reciever side domains these corrections come from Schaffer's surface wave models. Some transformations may need to be done to the model averages to
    get them to fit into our scheme
    '''