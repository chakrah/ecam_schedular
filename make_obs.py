import os
import sys
import numpy as np
import numpy as pd
import pandas as pd
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u

#########################################################################################################

default_pth= '/Users/hritam/EulerCam_Catalogues/automatic_obs/' # MODIFY THIS FOR YOUR PC

#########################################################################################################




if os.path.exists(default_pth):
    print('\n '"\x1B[3m" + 'Created OB will be saved in:', default_pth, "\x1B[3m")
else:
    default_pth = os.getcwd() + '/'
    print('\n'"\x1B[3m"+ 'Created OB will be saved in:', default_pth, "\x1B[3m")


target_name= input('\n Target name: ')

filter_name= input('\n Filter name (default - RG): ')

program_number= input('\n Program number (default - 449): ')


# setting defaults
if target_name == "":
    sys.exit('\n FAIL: Please provide a target name')
if filter_name == "":
    filter_name = 'RG'
if program_number == "":
    program_number = "449"


def make_OB(target_name, filter='RG', program=449,
            simbad_identifier=None, file_name=None, library_folder=None, 
            verbose=True):
    
    template = [
    {
        'code': '----',
        'type': '----',
        'alphacat': '--------',
        'deltacat': '--------',
        'alphatar': '--------',
        'deltatar': '--------',
        'alpharef': '--------',
        'deltaref': '--------',
        'equicat': '-------',
        'mualph': '------',
        'mudelt': '------',
        'sequence': '--------',
        'ampname': '-------',
        'mv': '--',
        'ut1': '---',
        'uti': '---',
        'ut2': '---',
        'ute': '---',
        'ut3': '---',
        'piscomod': '--------',
        'noprog': '------'
    },
    {
        'code': 'WASP-127',
        'type': 'CAM_ABTR',
        'alphacat': '10h42m14.0836s',
        'deltacat': '-03:50:06.260',
        'alphatar': '10h42m14.0836s',
        'deltatar': '-03:50:06.260',
        'alpharef': np.nan,
        'deltaref': np.nan,
        'equicat': 2000.0,
        'mualph': 0.000,
        'mudelt': 0.000,
        'sequence': 'RG/10/0.0/0.0/0.0/0.00',
        'ampname': 'ALL',
        'mv': 10.17,
        'ut1': -1,
        'uti': -1,
        'ut2': -1,
        'ute': -1,
        'ut3': -1,
        'piscomod': 'off',
        'noprog': 449
    }
    ]
    
    raw_ob= pd.DataFrame(template)
    

    raw_ob['code'][1] = target_name # setting target name
    
    if simbad_identifier != None:
        simbad_target_name = simbad_identifier
    else:
        simbad_target_name = target_name
    
    customSimbad = Simbad()
    customSimbad.add_votable_fields('flux(V)')    
    customSimbad.add_votable_fields('flux(R)')
    
    simbad_results = customSimbad.query_object(simbad_target_name)
    
    
    ra_str=simbad_results[0][1] # in h:m:s
    dec_str=simbad_results[0][2]# in d:m:s. Here d is in degree. 
    
    simbad_cat_ra= ra_str.split(' ')[0]+ 'h' + ra_str.split(' ')[1]+ 'm' + ra_str.split(' ')[2]+ 's'
    simbad_cat_dec= dec_str.split(' ')[0] + ':' + dec_str.split(' ')[1] + ':' + dec_str.split(' ')[2]
    
    simbad_target_ra= ra_str.split(' ')[0]+ 'h' + ra_str.split(' ')[1]+ 'm' + ra_str.split(' ')[2]+ 's'
    simbad_target_dec= dec_str.split(' ')[0] + ':' + dec_str.split(' ')[1] + ':' + dec_str.split(' ')[2]

    raw_ob['alphacat'][1] = simbad_cat_ra 
    raw_ob['deltacat'][1] = simbad_cat_dec
    raw_ob['alphatar'][1] = simbad_target_ra
    raw_ob['deltatar'][1] = simbad_target_dec
    
    
    
    # SETTING FILTER VALUE IN OB
    all_filters= ['UG', 'B1', 'BG', 'B2', 'V1', 'VG', 'GG', 'RG', 'IC', 'ZG', 'NG', 'OO']
    count=0
    for i in range(len(all_filters)):
        if filter == all_filters[i]:
            count+=1
        else:
            count+=0
       
    if count == 0:
        sys.exit('Specified filter is not available with EulerCam !')
    else:
        raw_ob['sequence'][1] = filter + raw_ob['sequence'][1][2:]
        
        
    
    # SETTING PROGRAM NUMBER IN OB
    all_programs= ['449', '447', '450', '500', '999']
    count=0
    for i in range(len(all_programs)):
        if str(program) == all_programs[i]:
            count+=1
        else:
            count+=0
            
    if count == 0:
        sys.exit('Specified program number is not available with EulerCam !')
    else:
        raw_ob['noprog'][1] =  str(program)
        
        
    
    # SETTING MAGNITUDE IN OB
    mag_v= simbad_results['FLUX_V'][0]
    mag_r= simbad_results['FLUX_R'][0]
    if np.isfinite(mag_v):
        raw_ob['mv'][1] =  np.round(mag_v,2)
    elif np.isfinite(mag_r):
        print('\n Setting R-mag instead of V-mag')
        raw_ob['mv'][1] = np.round(mag_r,2)
    else:
        print('\n Setting v-mag to default (10 mag)')
        raw_ob['mv'][1] = 10.0

    if verbose:
        print('\n')
        print(raw_ob)
        
    if file_name:
        if library_folder != None:
            raw_ob.to_csv(f'{library_folder}{file_name}.rdb', sep='\t', index=False)
        else:
            raw_ob.to_csv({file_name}.rdb, sep='\t', index=False)
    else:
        if library_folder != None:
            raw_ob.to_csv(f'{library_folder}{target_name}.rdb', sep='\t', index=False)
        else:
            raw_ob.to_csv(f'{target_name}.rdb', sep='\t', index=False)

    print(f'\n Generated OB file is successfully saved: {library_folder}{target_name}.rdb')
    print('#######################################################################')

    return


def check_target_name(target_name, ob):
    """check OB target name and target name"""
    
    ob_target_name= ob['code'][1]
    
    if target_name != ob_target_name:
        sys.exit('target name and OB target name doesn\'t match')
    else:
        print('Target and OB name check passed')
    return



def compare_corrdinates_with_simbad(ob, simbad_identifier=None):
    
    if simbad_identifier != None:
        ob_target_name = simbad_identifier
    else:
        ob_target_name= ob['code'][1]
        
    simbad_results = Simbad.query_object(ob_target_name)
    
    ra_str=simbad_results[0][1] # in h:m:s
    dec_str=simbad_results[0][2]# in d:m:s. Here d is in degree. 
    
    # simbad_target_ra= ra_str[0:2]+'h'+ra_str[3:5]+'m'+ra_str[6:]+'s'
    simbad_target_ra_indeg= float(ra_str.split(' ')[0]) * 15 + float(ra_str.split(' ')[1]) * (1/60) + float(ra_str.split(' ')[2]) * (1/3600)
    simbad_target_dec_indeg= float(dec_str.split(' ')[0]) + float(dec_str.split(' ')[1]) * (1/60) + float(dec_str.split(' ')[2]) * (1/3600)
    
    # expanding OB names
    split_ra_hour= ob['alphacat'][1].split('h')
    split_ra_min= split_ra_hour[1].split('m')
    split_ra_sec= split_ra_min[1].split('s')
    
    split_dec_deg= ob['deltacat'][1].split(':')
    ob_target_ra_indeg=  float(split_ra_hour[0]) * 15 + float(split_ra_min[0])  * (1/60) + float(split_ra_sec[0])  * (1/3600)
    ob_target_dec_indeg= float(split_dec_deg[0]) +      float(split_dec_deg[1]) * (1/60) + float(split_dec_deg[2]) * (1/3600)
    
    simbad_target_corr = SkyCoord(ra=simbad_target_ra_indeg*u.degree, dec=simbad_target_dec_indeg*u.degree, frame='icrs')   
    ob_target_corr= SkyCoord(ra=ob_target_ra_indeg*u.degree, dec=ob_target_dec_indeg*u.degree, frame='icrs')
    seperation_between_simbad_and_OB_indeg= simbad_target_corr.separation(ob_target_corr).degree * u.degree
    
    euler_fov_limit= 15 * (1/60) * u.degree
    
    if seperation_between_simbad_and_OB_indeg < euler_fov_limit/2:
        print(f'Target in EULER FoV and seperation from center is: {seperation_between_simbad_and_OB_indeg}')
    else:
        sys.exit('Target is not in the Field of View! Remake OB')
    return

# if len(target_name)> 1:
#     if len(filter_name) >1 and len(filter_name) != len(target_name):
#         sys.exit('Number of target and filters do not match!')
#     if len(filter_name) == 1:
#         for i in range(len(target_name)):
#             make_OB(target_name[i], filter_name, program_number, library_folder=default_pth)
#     else:
#         for i in range(len(target_name)):
#             make_OB(target_name[i], filter_name[i], program_number, library_folder=default_pth)

make_OB(target_name, filter_name, program_number, library_folder=default_pth)
