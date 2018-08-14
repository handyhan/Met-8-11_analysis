# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 15:02:55 2018

@author: Hannah.N
"""

import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import sys
import math
from MSG_read import extract_data,plot_scatter_fit,plot_single_map
from regression_functions import regression_c0,regression_free,odr_regression
from math import sin, cos, sqrt, atan2, radians
from scipy import stats
from datetime import datetime
from netCDF4 import Dataset
from time import gmtime, strftime



    
BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox"


def read_in_from_csv(dates):
    if os.path.isfile(BaseDir +"/txt_data/"+dates[0]+"_m8_all_pixels") & os.path.isfile(BaseDir +"/txt_data/"+dates[0]+"_m11_all_pixels") :
        m8_pixels = pd.read_csv(BaseDir +"/txt_data/"+dates[0]+"_m8_all_pixels")
        m11_pixels = pd.read_csv(BaseDir +"/txt_data/"+dates[0]+"_m11_all_pixels")
        """
        for date in dates[1:]:
            m8_pixels_day=pd.read_csv("./txt_data/"+date+"_m8_all_pixels")
            m11_pixels_day = pd.read_csv("./txt_data/"+date+"_m11_all_pixels")
            m8_pixels = pd.concat([m8_pixels,m8_pixels_day],axis=0, join='outer')
            m11_pixels = pd.concat([m11_pixels,m11_pixels_day],axis=0, join='outer')
        """
        
        print 
        return m8_pixels, m11_pixels
    else:
        print "Requested dates do not exist"




def proccess_grid(df,lon,lat, sat,resolution):
    
    
    lats_binned = np.digitize(df['LATITUDE'].values[:],lat)
    lon_binned = np.digitize(df['LONGITUDE'].values[:],lon)
    
    df_lat = pd.DataFrame({'LAT':lats_binned,'LONG':lon_binned})
    df['GRID_NO'] = df_lat.apply(lambda row: (str(row.LAT).zfill(3) + " " + str(row.LONG).zfill(3)) , axis=1)   
    
    df = df.sort_values(['GRID_NO'])
    fire_ids, pixel_count = np.unique(df['GRID_NO'], return_counts=True)
    df_ids = pd.DataFrame({'GRID_NO':fire_ids,'pixel_count':pixel_count})
    df['GRID_NO'] = df['GRID_NO'].apply(str)
    grid_df = pd.DataFrame(df_ids.GRID_NO.str.split(' ',1).tolist(), columns = ['lat_box','long_box'])
    grid_df.lat_box  = grid_df.lat_box.apply(lambda x: x.lstrip('0')).astype(int)
    grid_df.long_box  = grid_df.long_box.apply(lambda x: x.lstrip('0')).astype(int)
    
    frp_sum = df.groupby(['GRID_NO'])['FRP'].agg('sum')
    df['FRP_UNCERTAINTY_2'] = df['FRP_UNCERTAINTY'].apply(lambda x: x**2 )
    frp_err = np.sqrt(df.groupby(['GRID_NO'])['FRP_UNCERTAINTY_2'].agg('sum'))
    mean_vza = df.groupby(['GRID_NO'])['PIXEL_VZA'].mean() 
    
    lat_max = grid_df.lat_box.apply(lambda x:lat[(x)])
    lat_min= grid_df.lat_box.apply(lambda x:lat[(x-1)])   
    lat_center = lat_min + resolution/2
    
    long_max = grid_df.long_box.apply(lambda x:lon[(x)])
    long_min = grid_df.long_box.apply(lambda x:lon[(x-1)])
    long_center = (long_max + long_min)/2
    
    grid = pd.DataFrame({'GRID_NO':fire_ids, 'summed_FRP': frp_sum.values[:],'FRP_uncertainty':frp_err.values[:], 'pixel_count':pixel_count , 'SAT': sat,'LAT_MIN':lat_min,'LAT_MAX':lat_max,'LAT_CENTER':lat_center, 'LONG_MIN':long_min,'LONG_MAX':long_max,'LONG_CENTER':long_center, })

    vza_1deg = pd.read_csv(BaseDir + "/txt_data/MET_" + sat + "_1Deg_VZA")
    
    grid['loc_marker'] = grid['LAT_CENTER']*1000 + grid['LONG_CENTER']
    vza_1deg['loc_marker'] = vza_1deg['LATITUDE']*1000 + vza_1deg['LONGITUDE']
    grid = grid.sort_values(['loc_marker'])
    unique = np.unique(grid['loc_marker'])
    vza_1deg = vza_1deg[vza_1deg['loc_marker'].isin(unique)]  
    vza_1deg = vza_1deg.sort_values(['loc_marker'])
    grid['VZA'] = vza_1deg['VZA'].values[:]
    grid = grid.drop(columns=['loc_marker'])
    
    return grid,fire_ids
    


def write_gridded_file(lonmin,lonmax,latmin,latmax,m8_df,m11_df):

   
    m8_df = m8_df[(m8_df['LATITUDE'] > latmin) & (m8_df['LATITUDE'] < latmax) ]
    m8_df = m8_df[(m8_df['LONGITUDE'] > lonmin) & (m8_df['LONGITUDE'] < lonmax) ]
    m8_df = m8_df.reset_index(drop=True)
    m11_df = m11_df[(m11_df['LATITUDE'] > latmin) & (m11_df['LATITUDE'] < latmax) ]
    m11_df = m11_df[(m11_df['LONGITUDE'] > lonmin) & (m11_df['LONGITUDE'] < lonmax) ]
    m11_df = m11_df.reset_index(drop=True)
    
    m8_grid,grid_boxes_m8 = proccess_grid(m8_df,lon,lat,'8',resolution)
    m11_grid,grid_boxes_m11 = proccess_grid(m11_df,lon,lat,'11',resolution) 
     
    joint_boxes =  np.concatenate([grid_boxes_m8,grid_boxes_m11])
    joint_boxes = pd.Series(data=joint_boxes )
    joint_boxes = joint_boxes[joint_boxes.duplicated()] 
    m8_grid = m8_grid[m8_grid['GRID_NO'].isin(joint_boxes.values[:])]
    m11_grid = m11_grid[m11_grid['GRID_NO'].isin(joint_boxes.values[:])]
    
    grid_df = pd.DataFrame(joint_boxes.str.split(' ',1).tolist(), columns = ['lat_box','long_box'])
    grid_df.lat_box  = (grid_df.lat_box.apply(lambda x: x.lstrip('0')).astype(int))-1
    grid_df.long_box  = (grid_df.long_box.apply(lambda x: x.lstrip('0')).astype(int))-1
    
    #regrid_vza(lonmin,lonmax,latmin,latmax)
    grid = pd.concat([m8_grid,m11_grid])
    
    return grid
    


resolution = 1.0
lonmin = -20.0
lonmax = 50.0
latmin = -35.0
latmax = 35.0
lon = np.arange(lonmin,(lonmax+resolution),resolution)
lat = np.arange(latmin,(latmax+resolution),resolution)
####################


    #m8_grid = pd.DataFrame({'MIN_LAT':lat_min,'MAX_LAT':lat_max,'MIN_LONG':lon_min,'MAX_LONG':lon_max,'GRID_NO': None,'summed_FRP': 0,'FRP_uncertainty':0,'mean_vza': 0, 'pixel_count':0 , 'SAT': 8})


"""
dates = pd.date_range(start='2018-06-09',end='2018-07-15')
dates = dates.format(formatter=lambda x: x.strftime('%Y%m%d'))
for date in dates:
    print "NOTE: Generating fires for "  + date
    times = pd.date_range(start='00:00',end='23:45', freq='15min')
    times=times.format(formatter=lambda x: x.strftime('%H%M'))
    header_1=True
    for time in times:
        print "Proccesing per area FRP " + time
        m8_df,m11_df = extract_data(date,time)
        
        if m8_df.empty or m11_df.empty:
            print "ERROR while opening file"
            
        else:    
            #grid = 
            grid = write_gridded_file(lonmin,lonmax,latmin,latmax,m8_df,m11_df)
            #vza_1deg = pd.read_csv(BaseDir + "/txt_data/MET_" + "8" + "_1Deg_VZA")
            
            #plot_single_map(vza_1deg['LATITUDE'].values[:],vza_1deg['LONGITUDE'].values[:],vza_1deg['VZA'].values[:],'met8_#_VZA_map')
            #print vza_1deg.sort_values(['VZA'])
            grid['DATE'] = date
            grid['TIME'] = time
            mode = 'w' if header_1 else 'a'
            grid.to_csv(BaseDir +"/txt_data/AREA_grid_"+date, mode=mode,header=header_1,index=False) 
            header_1=False
            
"""           
#print grid 
#fit_and_plot_OLS_bins(grid)
#fit_and_plot_ODR_bins(grid)

#ax1 = grid.plot.scatter(x='M11_summed_FRP',y='M8_summed_FRP')
#plt.savefig('TEST')
 
"""
    dataset = Dataset('data/test.nc','w', format='NETCDF4_CLASSIC')
    
    lat = dataset.createDimension('lat', len(m8_df['LATITUDE']))
    lon = dataset.createDimension('lon', len(m8_df['LONGITUDE']))
    time = dataset.createDimension('time', None)
    
    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float64,('lat',))
    longitudes = dataset.createVariable('longitude', np.float64,('lon',)) 
    frp = dataset.createVariable('Summed_FRP', np.float32,('time','lat','lon'))
    frp_err = dataset.createVariable('FRP_uncertainty', np.float32,('time','lat','lon'))
    
"""
def regrid_vza(lonmin,lonmax,latmin,latmax):
    
    lon = np.arange(lonmin,(lonmax+resolution),resolution)
    lat = np.arange(latmin,(latmax+resolution),resolution)
    deg_1_df = pd.DataFrame({'long':lon,'lat':lat})
    deg_1_df['central_lat'] = deg_1_df['lat'] + 0.5 
    deg_1_df['central_long'] = deg_1_df['long'] + 0.5 
    
    m11_vza = BaseDir +"/SEVIRI_Static/HDF5_LSASAF_MSG_VZA_IODC-Disk_201808131045"
    m11_lat = BaseDir +"/SEVIRI_Static/HDF5_LSASAF_MSG_LAT_IODC-Disk_201711300000"
    m11_long = BaseDir +"/SEVIRI_Static/HDF5_LSASAF_MSG_LON_IODC-Disk_201711300000"
    
    m11_vza = h5py.File(m11_vza,'r+') 
    
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())    
    msg_rows = [0,500,1000,1500,2000,2500,3000]
    m11_lat_lon = pd.read_csv(BaseDir + "/txt_data/TEST_MSG8_lat_lon_")#+str(msg_rows[0]))
    #for file in msg_rows[1:]:
    #    m11_lat_lon_n=pd.read_csv(BaseDir + "/txt_data/MSG_lat_lon_"+str(file))
    #    m11_lat_lon =  m11_lat_lon.append(m11_lat_lon_n,ignore_index=True)

   
    #plot_single_map(m11_lat_lon['LONGITUDE'].values[:],m11_lat_lon['LATITUDE'].values[:],"PLOT_M8_PIX")
    m11_lat_lon = m11_lat_lon[(m11_lat_lon['LONGITUDE'] >= (lonmin*100)) & (m11_lat_lon['LONGITUDE'] <= (lonmax*100))]
    m11_lat_lon = m11_lat_lon[(m11_lat_lon['LATITUDE'] >= (latmin*100)) & (m11_lat_lon['LATITUDE'] <= (latmax*100))]
    m11_lat_lon = m11_lat_lon.sample(len(m11_lat_lon['LONGITUDE'])/50)
    m11_lat_lon.reset_index(drop=True,inplace=True)
    m11_lat_lon['LONGITUDE'] =m11_lat_lon['LONGITUDE']/100
    m11_lat_lon['LATITUDE'] =m11_lat_lon['LATITUDE']/100
    
    m11_lat_lon['lon_index'] = pd.cut(m11_lat_lon['LONGITUDE'],lon,labels=False,include_lowest=True).values[:]
    m11_lat_lon['lat_index'] = pd.cut(m11_lat_lon['LATITUDE'],lat,labels=False,include_lowest=True).values[:]
    m11_lat_lon.reset_index(drop=True,inplace=True)
    
    m11_lat_lon['VZA'] = m11_lat_lon.index    
    m11_lat_lon['VZA'] = (m11_lat_lon.VZA.apply(lambda x: m11_vza['VZA'][m11_lat_lon['row_ix'][x]][m11_lat_lon['col_ix'][x]]))/100
    m11_lat_lon = m11_lat_lon.sort_values(by=['lon_index','lat_index'],ascending=True)    
    m11_lat_lon['bin_id'] =  ((m11_lat_lon['lat_index']*1000)+m11_lat_lon['lon_index']).astype(int)

    vza_mean = m11_lat_lon.groupby(['bin_id'])['VZA'].mean()
    vza_grid = pd.DataFrame({'bin_id':np.unique(m11_lat_lon['bin_id'])})
    vza_grid['VZA'] = vza_mean.values[:]
    vza_grid['lat_index'] = (vza_grid['bin_id'].astype(np.float)/1000).astype(int)
    vza_grid['lon_index'] = vza_grid['bin_id'] - (vza_grid['lat_index']*1000)
    
    vza_grid['LATITUDE'] =  vza_grid.lat_index.apply(lambda x: deg_1_df['central_lat'].iloc[x])
    vza_grid['LONGITUDE'] =  vza_grid.lon_index.apply(lambda x: deg_1_df['central_long'].iloc[x])
    vza_grid = vza_grid[vza_grid['VZA'] > 0 ] 
    #vza_corrupt.to_csv("C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/txt_data/CORRUPT",header=True,index=True)
    vza_grid.to_csv("C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/txt_data/MET_8_1Deg_VZA_TEST2",header=True,index=True)
    
    print vza_grid[(vza_grid['VZA'] < 2) & (vza_grid['VZA'] > -2)]
    vza_grid.plot(x='LONGITUDE',y='LATITUDE',c = 'VZA', kind='scatter')    
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    
    
    
    
regrid_vza(lonmin,lonmax,latmin,latmax)
    
#m11_vza = pd.read_csv("C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/txt_data/MET_8_1Deg_VZA_TEST")
#print m11_vza[(m11_vza['VZA'] < 2) & (m11_vza['VZA'] > -2)]