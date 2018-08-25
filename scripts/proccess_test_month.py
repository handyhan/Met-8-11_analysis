import numpy as np
import os
import pandas as pd
from datetime import datetime, timedelta
from MSG_read_general import extract_data_test_month_msg
from general_area_analysis import  get_gridded, get_grid_bounds
from Omission_comission_V2 import get_oms_coms_at_time
from general_cluster_find import get_clusters_at_time
from get_col_lin import get_col_lin,get_lat_lon



## parameters --------------------




def assign_seviri_col_row(modislist):
    lat_lons = np.array([modislist.LONGITUDE.values[:],modislist.LATITUDE.values[:]])
    for i in range(0,len(lat_lons[0]),1):
        args = lat_lons[:,i]
        col,line = get_col_lin(args[0],args[1],'Full')
        if i == 0:
            col_lin_arr = np.array([[col,line]])
        else:
            col_lin_arr = np.append(col_lin_arr,[[col,line]],axis=0)
            
    return col_lin_arr

def extract_data_test_month_modis():
    BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox"
    modis = BaseDir + "/TEST_DATA/MODIS_201408_atm.csv"
    if  os.path.exists(modis):
        modis_df = pd.read_csv(modis)
        columns = ['sample','line','lontitude','latitude','FRP','vzen','year','month','day','time','MODIS']
        modis_df = modis_df[columns]
        modis_df['time'] = modis_df['time'].astype(int).astype(str)
        modis_df['time'] = modis_df['time'].apply(lambda x: x.zfill(4))      
        modis_df['hour'] = modis_df['time'].apply(lambda x: x[0:2])      
        modis_df['minute'] = modis_df['time'].apply(lambda x: x[2:4])     
        modis_df['TIME'] = pd.to_datetime(modis_df[['year', 'month', 'day','hour','minute']])
        modis_df = modis_df.drop(columns=['year', 'month', 'day','hour','time','minute'])        
        modis_df = modis_df.rename(columns={'line':'LINE','sample':'COL','latitude':'LATITUDE','lontitude':'LONGITUDE','vzen':'VZA','MODIS':'SAT'})
        modis_df = modis_df.sort_values('TIME')
        modis_df = modis_df.set_index(pd.DatetimeIndex(modis_df['TIME'])).drop(columns=['TIME'])
        col_lin = assign_seviri_col_row(modis_df)    
        modis_df['MSG_COL'] = col_lin[:,0]
        modis_df['MSG_LINE'] = col_lin[:,1]
            
        return modis_df
    elif not os.path.exists(modis):
        print "ERROR when trying to open file, this file does not exist:" + modis    
        modis_df = [pd.DataFrame()]
        return modis_df

BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox" 
#lon_max = 20
#lon_min = -15
#lat_max = -5
#lat_min = -10
#resolution = 1.0


lon_max = 50
lon_min = 10
lat_max = 0
lat_min = -20
resolution = 1.0

msg_df = extract_data_test_month_msg()
#modis_df = extract_data_test_month_modis()
#modis_df.to_csv(BaseDir + "/TEST_DATA/MODIS_with_MSG_columns",header=True,index=True) 

modis_df = pd.read_csv(BaseDir + "/TEST_DATA/MODIS_with_MSG_columns", parse_dates=True,index_col='TIME') 

modis_df = modis_df[(modis_df['LATITUDE'] > lat_min) & (modis_df['LATITUDE'] < lat_max) ]                  
modis_df = modis_df[(modis_df['LONGITUDE'] > lon_min) & (modis_df['LONGITUDE'] < lon_max) ]
msg_df = msg_df[(msg_df['LATITUDE'] > lat_min) & (msg_df['LATITUDE'] < lat_max) ]                  
msg_df = msg_df[(msg_df['LONGITUDE'] > lon_min) & (msg_df['LONGITUDE'] < lon_max) ]

#t = np.array([datetime(2014,8,1,12,15),datetime(2014,8,1,12,15)])                                 # supply one off times

header_1=True
header_2=True
header_3=True
header_4=True
header_5=True
header_6=True
header_7=True
t = np.arange(datetime(2014,8,1), datetime(2014,9,1), timedelta(minutes=15)).astype(datetime)
for k in range(0,len(t),1):
    
    print "running analysis for time ",t[k]
    #per-area
    """
    lat_bounds,lon_bounds = get_grid_bounds(lat_min,lat_max,lon_min,lon_max,resolution)  
    modis_grid_df = get_gridded(t[k],modis_df,lat_bounds,lon_bounds, 'MODIS', resolution)
    msg_grid_df = get_gridded(t[k],msg_df,lat_bounds,lon_bounds, '11', resolution)
    
    #Omission & commision
    sub_msg_processed,mod_omissions, mod_match = get_oms_coms_at_time(t[k], modis_df, msg_df)
    """
    #per-fire    
    mod_clusters_df = get_clusters_at_time(t[k], modis_df, 'LINE', 'COL')
    msg_clusters_df =get_clusters_at_time(t[k], msg_df, 'PIXEL_LINE', 'PIXEL_COL')
    
    """
    if modis_grid_df.empty:
        pass
    else:
        mode = 'w' if header_1 else 'a'
        modis_grid_df.to_csv(BaseDir + "/TEST_DATA/MODIS_gridded", mode=mode,header=header_1,index=True) 
        header_1=False
        
    if msg_grid_df.empty:
        pass
    else:
        mode = 'w' if header_2 else 'a'
        msg_grid_df.to_csv(BaseDir + "/TEST_DATA/MSG_gridded", mode=mode,header=header_2,index=True) 
        header_2=False
        
    if sub_msg_processed.empty:
        pass
    else:
        mode = 'w' if header_3 else 'a'
        sub_msg_processed.to_csv(BaseDir + "/TEST_DATA/correct_n_commission", mode=mode,header=header_3,index=True) 
        header_3=False
        
    if mod_omissions.empty:
        pass
    else:
        mode = 'w' if header_4 else 'a'
        mod_omissions.to_csv(BaseDir + "/TEST_DATA/omission", mode=mode,header=header_4,index=True) 
        header_4=False

    if mod_match.empty:
        pass
    else:
        mode = 'w' if header_5 else 'a'
        mod_match.to_csv(BaseDir + "/TEST_DATA/correct_MODIS", mode=mode,header=header_5,index=True) 
        header_5=False
    
    """    
    if mod_clusters_df.empty:
        pass
    else:
        mode = 'w' if header_6 else 'a'
        mod_clusters_df.to_csv(BaseDir + "/TEST_DATA/MODIS_fires", mode=mode,header=header_6,index=True) 
        header_6=False

    if msg_clusters_df.empty:
        pass
    else:
        mode = 'w' if header_7 else 'a'
        msg_clusters_df.to_csv(BaseDir + "/TEST_DATA/MSG_fires", mode=mode,header=header_7,index=True) 
        header_7=False
   
    
    
    
    
    
















