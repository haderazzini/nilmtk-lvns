from __future__ import print_function, division
import pandas as pd
import numpy as np
from copy import deepcopy
from os.path import join, isdir, isfile
from os import listdir
import re
import time, datetime, yaml
from sys import stdout
from nilmtk.utils import get_datastore
from nilmtk.datastore import Key
from nilmtk.timeframe import TimeFrame
from nilmtk.measurement import LEVEL_NAMES
from nilmtk.utils import get_module_directory, check_directory_exists
from nilm_metadata import convert_yaml_to_hdf5, save_yaml_to_datastore
import csv
from pandas.core.frame import DataFrame
import os

import matplotlib.pyplot as plt

"""
TODO:
* The bottleneck appears to be CPU.  So could be sped up by using 
  multiprocessing module to use multiple CPU cores to load lvns channels in 
  parallel.
"""


def convert_lvns(results_LVNS,output_filename, format='HDF',timezone='Canada/Mountain',file_groundtruth=''):
    """
    Parameters
    ----------
    lvns_path : str
        The root path of the alva low_freq dataset.
    output_filename : str
        The destination filename (including path and suffix).
    format : str
        format of output. Either 'HDF' or 'CSV'. Defaults to 'HDF'
    """

    # Open DataStore
    store = get_datastore(output_filename, format, mode='w')
    
    if file_groundtruth=='':
        trans_states=False
    else:
        trans_states=True
    
    houses,dict_met_buildings,states=organize_meters(results_LVNS['S_Hs'],results_LVNS['IPhase_Hs'],results_LVNS['app'],trans_states)
    save_metadata_lvns(results_LVNS['app'],dict_met_buildings)
    
    if trans_states:
        save_ground_truth(states,file_groundtruth,timezone)
        
    
    # Convert raw data to DataStore
    _convert(houses, store, timezone)

    # Add metadata
    try:
        save_yaml_to_datastore(join(get_module_directory(), 
                                  'dataset_converters', 
                                  'lvns', 
                                  'metadata'),
                             store)
    except:
        store.close()
        
    store.close()

    print("Done converting alva to HDF5!")

def save_ground_truth(states,folder,timezone):

    
    for house in states.keys():
        
        df_phase_1=pd.DataFrame()
        df_phase_2=pd.DataFrame()
        
        directory=folder+'house_'+str(house)+'\\'
        
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        for appliance in states[house].keys():
            states[house][appliance]['events'].to_csv(directory+appliance+'.csv')
            
            if states[house][appliance]['app_conn']==1:
                df_phase_1=pd.concat([df_phase_1,states[house][appliance]['events']])
            
            elif states[house][appliance]['app_conn']==2:
                df_phase_2=pd.concat([df_phase_2,states[house][appliance]['events']])
                
            else:#==3
                df_phase_1=pd.concat([df_phase_1,states[house][appliance]['events']])
                df_phase_2=pd.concat([df_phase_2,states[house][appliance]['events']])
    
        df_phase_1=df_phase_1.sort_index()
        df_phase_2=df_phase_2.sort_index()
        
        #df_phase_1.plot(style='*')
        #plt.show()

        df_phase_1=df_phase_1.tz_localize(timezone)
        df_phase_2=df_phase_2.tz_localize(timezone)

        
        df_phase_1.to_csv(directory+'phase_1.csv')
        df_phase_2.to_csv(directory+'phase_2.csv')
            
def _convert(houses,store,timezone):
    
    for house in houses.keys():
        print("Loading house", house, end="... ")
        for meter in houses[house].keys():
            key = Key(building=house, meter=meter)
            store.put(str(key),houses[house][meter].tz_localize(timezone))
            print(meter, end=" ")
        print(' ')

        
def build_power_curve(Appliance,sample_period=1):
    '''
    build a dataframe of active power for Appliance
    
    Input: Appliance,
    
    Output: dataframe 
    
    example: build_curve(Appliance)
    
    physical_quantity	power
    type	active
    2017-08-03 00:00:00	5.0
    2017-08-03 00:00:01	5.0
    2017-08-03 00:00:02	5.0
    
    '''
    
    indexs=create_index_today(Appliance['powratio'],sample_period)
    
    if Appliance['app_conn']==3:
        power_curve=[(Appliance['P']/2)*a for a in Appliance['powratio']]
    else:
        power_curve=[Appliance['P']*a for a in Appliance['powratio']]
    
    multiindex_columns = pd.MultiIndex.from_tuples([('power','active')], names=LEVEL_NAMES)
    return pd.DataFrame(data=power_curve, index=indexs, columns=multiindex_columns, dtype=np.float32)

def create_index_today(V,sample_period):
    '''Index are based in today'''
    
    today=datetime.datetime.now()
    dt_base = datetime.datetime(today.year,today.month,today.day,0,0)
    indexs=[dt_base+datetime.timedelta(seconds=i*sample_period) for i,a in enumerate(V)]
    return indexs

def find_states_trasitions(Vstates,step=0.1):
    
    ls=[Vstates[0]]+Vstates
    rs=Vstates+[Vstates[-1]]
    
    diff=[abs(a_i - b_i) for a_i, b_i in zip(ls, rs)]
    
    return [i for i,a in enumerate(diff) if a >= step]

def remove_low_power_delta(pos_trans,df_power,step=0.1,n_samples=3,power_limit=5):
    
    new_pos_trans=[]
    for transition in pos_trans:
        mean_pre=np.mean(df_power.ix[max(transition-n_samples,0):transition])
        mean_pos=np.mean(df_power.ix[transition+1:min(transition+n_samples,len(df_power)-1)])
        
        if abs(mean_pos.values-mean_pre.values) >= power_limit:
            new_pos_trans.append(transition)
        #else:
        #     print('low_diff')
        #    print(abs(mean_pos.values-mean_pre.values))
            
    return new_pos_trans

def phase_main(home_phase,sample_period=1):
    
    
    multiindex_columns_active = pd.MultiIndex.from_tuples([('power','active')], names=LEVEL_NAMES)
    multiindex_columns_reactive = pd.MultiIndex.from_tuples([('power','reactive')], names=LEVEL_NAMES)
    
    P=[S.real for S in home_phase]
    Q=[S.imag for S in home_phase]

    indexs=create_index_today(P,sample_period)
    
    dfp=pd.DataFrame(data=P, index=indexs, columns=multiindex_columns_active, dtype=np.float32)
    dfq=pd.DataFrame(data=Q, index=indexs, columns=multiindex_columns_reactive, dtype=np.float32)
        
    df_phase=pd.concat([dfp,dfq],axis=1)
        
    return df_phase

def concat_harmonics(phase_mains,phase,IPhase_Hs,app):
    df_phase=phase_mains.copy()
    
    harm_order=app.harm_min
    
    while harm_order <= app.harm_max:
        h = int((harm_order + 1.0) / 2.0) - 1

        multiindex_columns = pd.MultiIndex.from_tuples([( str(harm_order)+' harmonic current','')], names=LEVEL_NAMES)
        
        df_h=pd.DataFrame(data=IPhase_Hs[h][phase], index=df_phase.index, columns=multiindex_columns, dtype=np.float32)

        df_phase=pd.concat([df_phase,df_h],axis=1)
        
        harm_order=harm_order+app.harm_step
    return df_phase

def find_repeated_app(app):
    rept={}
    for i,House in enumerate(app.houses):
        apps_names=[]
        gab_names=[]
        rept[i+1]={}
        for Appliance in House['appliances']:
            gab_names.append(Appliance['app_name'])

            if not Appliance['app_name'] in apps_names:
                apps_names.append(Appliance['app_name'])

            elif Appliance['app_name'] in rept[i+1].keys():
                rept[i+1][Appliance['app_name']]=rept[i+1][Appliance['app_name']]+1
            else:
                rept[i+1][Appliance['app_name']]=2
    return rept

def check_repeated_app(name,house,rept):
    
    if rept[house]=={}:
        return 1,rept,name
    elif name not in rept[house].keys():
        return 1,rept,name
    else:
        instance=rept[house][name]
        new_name=name+' '+str(instance)
        if instance==1:
            del rept[house][name]
        else:
            
            rept[house][name]=rept[house][name]-1
        return instance,rept,new_name
    

def organize_meters(S_Hs,IPhase_Hs,app,tran_states=False):

    abrev={'CFL':'compact fluorescent lamp',#'Compact Fluorescent Lamp'
            'EBL':'fluorescent lamp',#'Electric-Ballast Lamp',
            'MBL':'fluorescent lamp',#'Magnetic-Ballast Lamp',
            'INC':'incandescent lamp',
            'PC':'desktop computer',#'Desktop PC',
            'LCD':'computer monitor', #'LCD Computer Monitor',
            'LAP':'laptop computer',#'Laptop',
            'LCDTV':'television',#'LCD Television',
            'CRTTV':'television',#'CRT Television',
            'RFR':'fridge',#'Regular Fridge',
            'FUR':'electric furnace',
            'ASDFR':'fridge',#'ASD-based Fridge',
            'FRE':'freezer',
            'WSH':'washing machine',#'ASD-based Washer',
            'DRY':'tumble dryer',#'Regular Dryer',
            'OVE':'electric oven',
            'MW':'microwave',
            'TOA':'toaster',
            'COF':'coffee maker',
            'FOO':'food processor',
            'VAC':'vacuum cleaner',
            'DSW':'dish washer',
            'STO':'stove'
          }
    
    rept=find_repeated_app(app)
    
    houses={}
    
    states={}
    
    dict_met_buildings={}

    for i in range(1,int(len(S_Hs)/2)+1):
        
        df_ph1=phase_main(S_Hs[2*i-2])
        df_ph1=concat_harmonics(df_ph1,2*i-2,IPhase_Hs,app)
        
        df_ph2=phase_main(S_Hs[2*i-1])
        df_ph2=concat_harmonics(df_ph2,2*i-1,IPhase_Hs,app)
        
        houses[i]={1:df_ph1,2:df_ph2}
        
        dict_met_buildings[i]={'instance':i,'original_name':'lvns_house_'+str(i),'elec_meters':{1:{'device_model': 'synthetic_main','site_meter':True},
                                                                                                2:{'device_model': 'synthetic_main','site_meter':True}}}

        states[i]={}
        
    for i,House in enumerate(app.houses):

        counter=max(houses[i+1].keys())+1
        dict_met_buildings[i+1]['appliances']=[]

        for Appliance in House['appliances']:
            
            instance,rept,name=check_repeated_app(Appliance['app_name'],i+1,rept)
            
            houses[i+1][counter]=build_power_curve(Appliance,sample_period=1)
            
            if tran_states: #to do a groundtruth of detection
                
                ptrans=find_states_trasitions(Appliance['states'],step=0.1)
                #print(Appliance['app_name'])
                ptrans=remove_low_power_delta(ptrans,houses[i+1][counter],step=0.1,n_samples=3,power_limit=5)
                
                indexs_trans=houses[i+1][counter].ix[ptrans].index
                df_ground=pd.DataFrame(data=counter, index=indexs_trans, columns=['Appliance'], dtype=np.float32)
                df_ground.index.name='Events'
                states[i+1][name]={'app_conn':Appliance['app_conn'],
                                        'events':df_ground}
                
                #plt.plot(Appliance['states'])
                #plt.plot(Appliance['powratio'])

                #plt.plot(ptrans,[0]*len(ptrans),'*')
                #plt.show()
                
            else:
                states={}
            
            if Appliance['app_conn']==1:
                dict_met_buildings[i+1]['elec_meters'][counter]={'submeter_of':1, 'device_model': 'synthetic_meter'}
                dict_met_buildings[i+1]['appliances'].append({'instance':instance,'meters':[counter],'original_name':name,'type':abrev[Appliance['app_name']]})

            elif Appliance['app_conn']==2:
                dict_met_buildings[i+1]['elec_meters'][counter]={'submeter_of':2, 'device_model': 'synthetic_meter'}
                dict_met_buildings[i+1]['appliances'].append({'instance':instance,'meters':[counter],'original_name':name,'type':abrev[Appliance['app_name']]})

            else:#Appliance['app_conn']==3
                dict_met_buildings[i+1]['elec_meters'][counter]={'submeter_of':1, 'device_model': 'synthetic_meter'}
                counter+=1
                dict_met_buildings[i+1]['elec_meters'][counter]={'submeter_of':2, 'device_model': 'synthetic_meter'}
                houses[i+1][counter]=build_power_curve(Appliance,sample_period=1)
                dict_met_buildings[i+1]['appliances'].append({'instance':instance,'meters':[counter-1,counter],'original_name':name,'type':abrev[Appliance['app_name']]})

            counter+=1
            
    return houses,dict_met_buildings,states

def save_metadata_lvns(app,dict_met_buildings):
    
    directory=join(get_module_directory(), 
                              'dataset_converters', 
                              'lvns', 
                              'metadata')
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    for house in dict_met_buildings.keys():
    
        with open(directory+'\\building'+str(house)+'.yaml', 'w') as outfile:
            yaml.dump(dict_met_buildings[house], outfile, default_flow_style=False)
    
    
    
    #--------------------------Meter Devices---------------------------------
    meter_devices={'synthetic_main':{'description': 'Measures main meter',
     'max_sample_period': 60,
     'measurements': [{'physical_quantity': 'power','type': 'active','upper_limit': 3200000,'lower_limit': 0},
                     {'physical_quantity': 'power','type': 'reactive','upper_limit': 3200000,'lower_limit': -320000}
                     ],
     'model': 'Synthetic',
     'sample_period': 1},

    'synthetic_meter':{'description': 'Measures loads',
     'max_sample_period': 60,
     'measurements': [{'physical_quantity': 'power','type': 'active','upper_limit': 3200000,'lower_limit': 0}],
     'model': 'Synthetic',
     'sample_period': 1}}
    
    harm_order=app.harm_min
    
    while harm_order <= app.harm_max:
        
        #h = int((harm_order + 1.0) / 2.0) - 1
    
        new_harm={'physical_quantity': str(harm_order)+' harmonic current','upper_limit': 3200000,'lower_limit': -320000}
        
        #print(new_harm)
        meter_devices['synthetic_main']['measurements'].append(new_harm)
        harm_order=harm_order+app.harm_step
    
    with open(directory+'\\meter_devices.yaml', 'w') as outfile:
        yaml.dump(meter_devices, outfile, default_flow_style=False)
        
    
    #--------------------------Dataset---------------------------------
    today=datetime.datetime.now()

    dataset={'contact': 'torquato@ieee.org',
     'creators': [['Torquato, Ricador'],['Shi','Qingxin'],['Xu','Wilsun'],['Freitas','Walmir']],
     'description': 'A Monte Carlo Simulation Platform for Studying Low Voltage Residential Networks',
     'geo_location': {'country': 'Canada',
      'latitude': None,
      'locality': None,
      'longitude': None},
     'institution': 'UNIVERSITY OF ALBERTA AND UNIVERSTITY OF CAMPINAS (UNICAMP)',
     'long_name': 'Low Voltage Network Simulator',
     'name': 'LVNS',
     'number_of_buildings': len(app.houses),
     'publication_date': today.year,
     'related_documents': ['A Monte Carlo Simulation Platform for Studying Low Voltage Residential Networks','http://ieeexplore.ieee.org/abstract/document/6853399/?reload=true'],
     'schema': 'https://github.com/nilmtk/nilm_metadata/tree/v0.2',
     'subject': 'Disaggregated power demand from domestic buildings.',
     'timezone': 'Canada/Mountain'}

    with open(directory+'\\dataset.yaml', 'w') as outfile:
        yaml.dump(dataset, outfile, default_flow_style=False)