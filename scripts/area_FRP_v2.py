# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 15:02:55 2018

@author: Hannah.N
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import sys
import math
from MSG_read import extract_data,plot_compare,plot_scatter_fit
from math import sin, cos, sqrt, atan2, radians
from scipy import stats
from datetime import datetime
from netCDF4 import Dataset


def read_in_from_csv(dates):
    if os.path.isfile("./txt_data/"+dates[0]+"_m8_all_pixels") & os.path.isfile("./txt_data/"+dates[0]+"_m11_all_pixels") :
        m8_pixels = pd.read_csv("./txt_data/"+dates[0]+"_m8_all_pixels")
        m11_pixels = pd.read_csv("./txt_data/"+dates[0]+"_m11_all_pixels")
        """
        for date in dates[1:]:
            m8_pixels_day=pd.read_csv("./txt_data/"+date+"_m8_all_pixels")
            m11_pixels_day = pd.read_csv("./txt_data/"+date+"_m11_all_pixels")
            m8_pixels = pd.concat([m8_pixels,m8_pixels_day],axis=0, join='outer')
            m11_pixels = pd.concat([m11_pixels,m11_pixels_day],axis=0, join='outer')
        """
        return m8_pixels, m11_pixels
    else:
        print "Requested dates do not exist"




def proccess_grid(df):
    
    
    fire_pixels, pixel_count = np.unique(df['GRID_NO'], return_counts=True)
    print  fire_pixels, pixel_count
    
    df['GRID_NO']
    print df.apply(lambda row: int(row.GRID_NO[2:4]) , axis=1) #########################
    
    
    
    for j in fire_grid_ix:
        frp_sum = df.groupby(['FRP_GRID']).get_group(j)['FRP'].agg('sum')
        sum_sq_err = df.groupby(['FRP_GRID']).get_group(j)['FRP_UNCERTAINTY'].apply(lambda x: x**2 ).agg('sum')
        mean_vza = df.groupby(['PIXEL_VZA']).get_group(j).mean()
        sum_sq_err = np.sqrt(sum_sq_err)
        fire_grid.at[j,'SUMMED_FRP'] = frp_sum
        fire_grid.at[j,'FRP_UNCERTAINTY'] = sum_sq_err
        fire_grid.at[j,'PIXEL_VZA'] = mean_vza
   
    grid.update(fire_grid)
    return grid
        

def eliminate_zero_grid(m8_grid,m11_grid):
    
    m8_grid_list = np.unique((m8_grid[m8_grid['GRID_NO'] >= 0 ])['GRID_NO']) 
    m11_grid_list = np.unique((m11_grid[m11_grid['GRID_NO'] >= 0 ])['GRID_NO']) 
    grid_list = np.unique(np.concatenate((m11_grid_list,m8_grid_list)))
    
    m8_grid = m8_grid.loc[grid_list]
    
    m11_grid = m11_grid.loc[grid_list]
    return m8_grid, m11_grid
    


def write_gridded_file(dates):
    stepsize = 1   #
    lonmin = 10.0
    lonmax = 30.0
    latmin = -10.0
    latmax = 10.0
    nx = (lonmax - lonmin)/stepsize
    ny = (latmax - latmin)/stepsize
    lon = np.linspace(lonmin,lonmax,nx)
    lat = np.linspace(latmin,latmax,ny)
    print len(lon)
    
    
    m8_df, m11_df = read_in_from_csv(dates)
    m8_df = m8_df[(m8_df['LATITUDE'] >= latmin) & (m8_df['LATITUDE'] <= latmax) ]
    m8_df = m8_df[(m8_df['LONGITUDE'] >= lonmin) & (m8_df['LONGITUDE'] <= lonmax) ]
    m8_df = m8_df.reset_index(drop=True)
    
    lats_binned = np.digitize(m8_df['LATITUDE'].values[:],lat)
    lon_binned = np.digitize(m8_df['LONGITUDE'].values[:],lon)
    
    df_lat = pd.DataFrame({'LAT':lats_binned,'LONG':lon_binned})
    m8_df['GRID_NO'] = df_lat.apply(lambda row: (str(row.LAT).zfill(3) + str(row.LONG).zfill(3)) , axis=1)
    print m8_df['GRID_NO']
    
    
    m8_grid = proccess_grid(m8_df) 
    
    """
    m8_grid = pd.DataFrame({'MIN_LAT':lat_min,'MAX_LAT':lat_max,'MIN_LONG':lon_min,'MAX_LONG':lon_max,'GRID_NO': None})
    m11_grid = pd.DataFrame({'MIN_LAT':lat_min,'MAX_LAT':lat_max,'MIN_LONG':lon_min,'MAX_LONG':lon_max,'GRID_NO': None})
    m8_grid['C_LATITUDE'] = m8_grid.apply(lambda row: (row.MAX_LAT + row.MIN_LAT)/2 , axis=1)
    m11_grid['C_LATITUDE'] = m11_grid.apply(lambda row: (row.MAX_LAT + row.MAX_LAT)/2 , axis=1)    
    m8_grid['C_LONGITUDE'] = m8_grid.apply(lambda row: (row.MAX_LONG + row.MIN_LONG)/2 , axis=1)
    m11_grid['C_LONGITUDE'] = m11_grid.apply(lambda row: (row.MAX_LONG + row.MIN_LONG)/2 , axis=1)
    
    
    m8_df, m11_df = read_in_from_csv(dates)
    m8_df = m8_df[(m8_df['LATITUDE'] >= latmin) & (m8_df['LATITUDE'] <= latmax) ]
    m8_df = m8_df[(m8_df['LONGITUDE'] >= lonmin) & (m8_df['LONGITUDE'] <= lonmax) ]
    m8_df = m8_df.reset_index(drop=True)
    
    m11_df = m11_df[(m11_df['LATITUDE'] >= latmin) & (m11_df['LATITUDE'] <= latmax) ]
    m11_df = m11_df[(m11_df['LONGITUDE'] >= lonmin) & (m11_df['LONGITUDE'] <= lonmax) ]
    m11_df = m11_df.reset_index(drop=True)
    
    
    m11_df= assign_pixels(m11_df, m11_grid)
    m8_df = assign_pixels(m8_df, m8_grid)
    
    m8_grid = proccess_grid(m8_df, m8_grid) 
    
    print m8_grid
    m11_grid = proccess_grid(m11_df, m11_grid) 
    
    #m8_grid,m11_grid = eliminate_zero_grid(m8_grid,m11_grid)
    
    #m8_grid.to_csv("./txt_data/"+date+"_m11_all_grid",header=True,index=True) 
    #m11_grid.to_csv("./txt_data/"+date+"_m8_all_grid",header=True,index=True) 
    """

dates = ["20180609"]#,"20180611"]
write_gridded_file(dates)





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