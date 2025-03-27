from PyAstronomy import pyasl
import datetime as dt
import sys
import pandas as pd
import numpy as np
import matplotlib
import math
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import warnings
from astropy.utils.exceptions import AstropyUserWarning
warnings.simplefilter('ignore', category=AstropyUserWarning)
warnings.simplefilter('ignore', category=FutureWarning)


def load_targets_from_sheets(lib_csv):
    print(lib_csv)
    sheet= pd.read_csv(lib_csv, sep=',', dtype='str', skiprows=[0,1,2])
    
    all_targets= sheet.iloc[:,0]
    all_obs_nights= sheet.iloc[:,1]
    
    all_targets_pyasl=[]
    all_obs_year=[]
    all_obs_month=[]
    all_obs_day=[]
    
    for i, target in enumerate(all_targets):
        
        target_name_pyasl= target[:-1] + ' ' + target[-1] # target name suitable for pysal format
        all_targets_pyasl.append(target_name_pyasl)
        
        year= '20'+ all_obs_nights[i][-2:]
        all_obs_year.append(year)
        month= all_obs_nights[i][-5:-3]
        all_obs_month.append(month)
        day= all_obs_nights[i][0:2]
        all_obs_day.append(day)
        
    d= {'Target': all_targets_pyasl, 
        'Year': all_obs_year, 
        'Month': all_obs_month, 
        'Day': all_obs_day}
    df= pd.DataFrame(d)
        
    return df

def predict_transit_from_sheets(observation_df, obsOffset=1.5, obsOffset_plotting=0.5):
    
    # observation_df= observation_df.drop() # drop certain target from sheets. 
    
    n_targets= observation_df.shape[0]
    
    output_csv= np.empty(shape=(n_targets, 7), dtype=object)
    
    for i in range(n_targets):

        d= dt.datetime(year=int(observation_df['Year'][i]), 
                    month=int(observation_df['Month'][i]),
                    day=int(observation_df['Day'][i]),
                    hour=00)
        date_string= observation_df['Year'][i] + '_' + observation_df['Month'][i] + '_' + observation_df['Day'][i]
        
        jd= pyasl.jdcnv(d)
        
        planet= NasaExoplanetArchive.query_object(observation_df['Target'][i])
        
        data= {'ra': planet['ra'].value,
            'dec': planet['dec'].value, 
            'T0': planet['pl_tranmid'].value, 
            'orbPer': planet['pl_orbper'].value, 
            'orbInc': planet['pl_orbincl'].value, 
            'SMA': planet['pl_orbsmax'].value,
            'RpJ': planet['pl_radj'].value, 
            'RsSun': planet['st_rad'].value,
            'Tdur': planet['pl_trandur'].value/24, 
            }
        
        orbital_params= pd.DataFrame(data)
        
        ind= np.where(np.all(np.isfinite(orbital_params), axis=1))[0]
        # print('Number of reference\'s available:', len(ind))
        planet_name_row= pd.Series({ 'plName': planet['pl_name'][0]})
        # orbital_params= orbital_params.iloc[ind[0],].append(planet_name_row)
        orbital_params = pd.concat([orbital_params.iloc[ind[0]], planet_name_row])
        
        T0_HJD= pyasl.asl.astroTimeLegacy.helio_jd(date=orbital_params['T0']-2.4e6, ra= orbital_params['ra'], dec= orbital_params['dec']) 
        orbital_params['T0'] = float(T0_HJD) + 2.4e6
        
        orbital_params= orbital_params.to_dict()
        dat = pyasl.transitTimes(jd, jd+2., orbital_params,
                                observatory="esolasilla", obsOffset=obsOffset_plotting/24.,
                                minAltitude=10.0, fileOutput='Output/'+ planet['pl_name'][0] + '_' + date_string +'.txt')
        
        all_transits_info=[]
        for j in range(len(dat)):
            j=j+1
            obs_epoch= dat[j]
            day_or_night= obs_epoch['Twilight']
        
            transit_start= pyasl.daycnv(obs_epoch['Obs jd'][1] - orbital_params['Tdur']/2)
            transit_end= pyasl.daycnv(obs_epoch['Obs jd'][1] + orbital_params['Tdur']/2)
            transit_start= time_conversion_pretty(transit_start)
            transit_end= time_conversion_pretty(transit_end)
            
            obs_start= pyasl.daycnv(obs_epoch['Obs jd'][1] - orbital_params['Tdur']/2 - obsOffset/24.)
            obs_end= pyasl.daycnv(obs_epoch['Obs jd'][1] + orbital_params['Tdur']/2 + obsOffset/24.)
        
            # obs_start, transit_mid, obs_end= pyasl.daycnv(obs_epoch['Obs jd'])
            transit_mid= pyasl.daycnv(obs_epoch['Obs jd'][1])
            obs_start= time_conversion_pretty(obs_start)
            transit_mid= time_conversion_pretty(transit_mid)
            obs_end= time_conversion_pretty(obs_end)
        
            dat_info_string= f'({day_or_night}) Start_obs: {obs_start} :: Transit_start: {transit_start} :: Transit_mid {transit_mid} :: Transit_end {transit_end} :: Obs_end: {obs_end}'
            all_transits_info.append(dat_info_string)

        np.savetxt('Output/' + planet['pl_name'][0] + '_' + date_string +'.txt', all_transits_info, fmt='%s')
        
        
        if len(dat) == 0.0:
            print('\n No Transit !')
            output_csv[i]= observation_df['Target'][i], 'No Transit found'
        else:
            string_name= 'airmass_plots_P114/' + planet['pl_name'][0] + '_' + date_string
            pyasl.transitVisibilityPlot(dat, showMoonDist=True, markTransit=True, print2file=f'{string_name}.png') 
            output_csv[i] = observation_df['Target'][i], 'PASSED', obs_start, transit_start, transit_mid, transit_end, obs_end
        
    output_csv_df= pd.DataFrame(output_csv)
    header_string= ['Target', 'Status', 'Obs. Start', 'Transit Start', 'Transit Mid', 'Transit End', 'Obs. End']
    output_csv_df.to_csv('output.csv', sep='\t', index=False, header=header_string)
    
    return

def predict_transit(target_name='HIP67522 b', date='2023-04-30', n_days=2., obsOffset=1., output_pth=None):
    
    planet= NasaExoplanetArchive.query_object(target_name)
    year, month, day= date.split('-')
    
    d= dt.datetime(year=int(year), 
                    month=int(month),
                    day=int(day),
                    hour=00)
    
    jd= pyasl.jdcnv(d)
    jd_end= jd+n_days
    
    greg_start= pyasl.daycnv(jd)
    greg_end= pyasl.daycnv(jd_end)
    print("\x1B[3m" + f'\n Looking for transits from {greg_start[0]}-{greg_start[1]}-{greg_start[2]} to {greg_end[0]}-{greg_end[1]}-{greg_end[2]}'+ "\x1B[3m")
    

    if target_name == 'WASP-39 b':
        print('JWST emphemeris !')
        data= {'ra': 217.3266477,
            'dec': -3.4444994, 
            'T0': 2458441.5714062033, 
            'orbPer': 4.3242697595, 
            'orbInc': 87.7369, 
            'SMA': 0.0474071,
            'RpJ': 1.270000, 
            'RsSun': 0.895000,
            'Tdur': 2.87220432/24, 
            }
        
    elif target_name == 'TOI-942 b':
        
        print('TOI-942 b emphemeris !')
        data= {'ra': 76.6500,
            'dec': -20.2456, 
            'T0': 2458441.5714062033, 
            'orbPer': 4.3242697595, 
            'orbInc': 87.11, 
            'SMA': 0.0398,
            'RpJ': 0.3784, 
            'RsSun': 0.893,
            'Tdur': 3.66/24, 
            }
    elif target_name == 'TOI-942 c':
        print('TOI-942 c emphemeris !')
        data= {'ra': 76.6500,
            'dec': -20.2456, 
            'T0': 2458447.0633021118, 
            'orbPer': 10.1560847252, 
            'orbInc': 87.73418531, 
            'SMA': 0.07026518,
            'RpJ': 0.4276, 
            'RsSun': 0.893,
            'Tdur': 4.46/24, 
            }
    else:
        data= {'ra': planet['ra'].value,
            'dec': planet['dec'].value, 
            'T0': planet['pl_tranmid'].value, 
            'orbPer': planet['pl_orbper'].value, 
            'orbInc': planet['pl_orbincl'].value, 
            'SMA': planet['pl_orbsmax'].value,
            'RpJ': planet['pl_radj'].value, 
            'RsSun': planet['st_rad'].value,
            'Tdur': planet['pl_trandur'].value/24, 
            }
    
    orbital_params= pd.DataFrame([data])
    
    ind= np.where(np.all(np.isfinite(orbital_params), axis=1))[0]
    print('Number of (complete) reference\'s available:', len(ind))
    planet_name_row= pd.Series({ 'plName': planet['pl_name'][0]})
    # orbital_params= orbital_params.iloc[ind[0],].append(planet_name_row)
    orbital_params = pd.concat([orbital_params.iloc[ind[0]], planet_name_row])

        
    T0_HJD= pyasl.asl.astroTimeLegacy.helio_jd(date=orbital_params['T0']-2.4e6, ra= orbital_params['ra'], dec= orbital_params['dec']) 
    orbital_params['T0'] = float(T0_HJD) + 2.4e6
        
    orbital_params= orbital_params.to_dict()
    
    
    if output_pth:
        txt_file= output_pth + planet['pl_name'][0] + '_' + target_name +'.txt'
    else:
        txt_file= planet['pl_name'][0] + '_' + target_name +'.txt'
        
        
    dat = pyasl.transitTimes(jd, jd+n_days, orbital_params,
                            observatory="esolasilla", obsOffset=obsOffset/24.,
                            minAltitude=10.0, fileOutput= txt_file
                            )
    all_transits_info=[]
    for i in range(len(dat)):
        i=i+1
        obs_epoch= dat[i]
        day_or_night= obs_epoch['Twilight']
        
        transit_start= pyasl.daycnv(obs_epoch['Obs jd'][1] - orbital_params['Tdur']/2)
        transit_end= pyasl.daycnv(obs_epoch['Obs jd'][1] + orbital_params['Tdur']/2)
        transit_start= time_conversion(transit_start)
        transit_end= time_conversion(transit_end)
        
        obs_start, transit_mid, obs_end= pyasl.daycnv(obs_epoch['Obs jd'])
        
        obs_start= time_conversion(obs_start)
        transit_mid= time_conversion(transit_mid)
        obs_end= time_conversion(obs_end)
        
        dat_info_string= f'({day_or_night}) Start_obs: {obs_start} :: Transit_start: {transit_start} :: Transit_mid {transit_mid} :: Transit_end {transit_end} :: Obs_end: {obs_end}'
        all_transits_info.append(dat_info_string)

    np.savetxt('output.txt', all_transits_info, fmt='%s', delimiter='\t')
    if len(dat) == 0.0:
        print('\n No Transit !')
    else:
        if output_pth:
            png_file= output_pth + f'{target_name}.png'
        else:
            png_file= f'{target_name}.png'
        pyasl.transitVisibilityPlot(dat, showMoonDist=True, markTransit=True, print2file=png_file) 
    
    return

def time_conversion(time_array):
    year= time_array[0]
    month= time_array[1]
    day= time_array[2]

    hour_minute, hour= math.modf(time_array[3]) 
    minute_second, minute= math.modf(hour_minute*60)
    second= minute_second * 60
    hour= int(hour)
    minute= int(minute)#np.round(hour_minute*60., 2)
    second= int(second)
    
    date_string= str(year) +'-'+ str(month) +'-'+ str(day) +'_'+ str(hour)+':'+str(minute)+':'+str(second)
    return date_string

def time_conversion_pretty(time_array):
    
    hour_minute, hour= math.modf(time_array[3]) 
    minute_second, minute= math.modf(hour_minute*60)
    
    if hour < 10:
        hour= '0'+str(int(hour))
    else:
        hour= int(hour)
        
    if minute < 10:
        minute= '0' + str(int(minute))
    else:
        minute= int(minute)#np.round(hour_minute*60., 2)
    
    date_string= str(hour)+':'+str(minute)
    return date_string

# sheet_name= '/Users/hritam/Documents/PhD/1. EulerCam/scheduling/emphemeris-predictor/observation_sheets/ATMOS_449_P114.csv'
# observation_df= load_targets_from_sheets(sheet_name)

# predict_transit_from_sheets(observation_df, obsOffset=1.5)

# predict_transit(target_name='TOI-1130 c', date='2024-06-06', n_days=5.,)# output_pth='individual_obs/')
predict_transit(target_name='TOI-942 c', date='2024-10-01', obsOffset=1., n_days=30*4)





