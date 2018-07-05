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
from MSG_read import extract_data,plot_compare,plot_scatter_fit
from regression_functions import regression_c0,regression_free,odr_regression
from math import sin, cos, sqrt, atan2, radians
from scipy import stats
from datetime import datetime
from netCDF4 import Dataset


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
    
    grid = pd.DataFrame({'GRID_NO':fire_ids, 'summed_FRP': frp_sum.values[:],'FRP_uncertainty':frp_err.values[:],'mean_vza': mean_vza.values[:], 'pixel_count':pixel_count , 'SAT': sat,'LAT_MIN':lat_min,'LAT_MAX':lat_max,'LAT_CENTER':lat_center, 'LONG_MIN':long_min,'LONG_MAX':long_max,'LONG_CENTER':long_center, })
    return grid, fire_ids

        

def eliminate_zero_grid(m8_grid,m11_grid):
    
    m8_grid_list = np.unique((m8_grid[m8_grid['GRID_NO'] >= 0 ])['GRID_NO']) 
    m11_grid_list = np.unique((m11_grid[m11_grid['GRID_NO'] >= 0 ])['GRID_NO']) 
    grid_list = np.unique(np.concatenate((m11_grid_list,m8_grid_list)))
    
    m8_grid = m8_grid.loc[grid_list]
    
    m11_grid = m11_grid.loc[grid_list]
    return m8_grid, m11_grid
    


def write_gridded_file(m8_df,m11_df):

    
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
    """
    m8_grid['M8_summed_FRP'] = m8_grid['summed_FRP'] 
    m8_grid['M8_FRP_uncertainty'] = m8_grid['FRP_uncertainty'] 
    m8_grid['M8_PIXEL_COUNT'] = m8_grid['pixel_count'].astype(int)
    
    m8_grid['M11_summed_FRP'] = m11_grid['summed_FRP'] 
    m8_grid['M11_FRP_uncertainty'] = m11_grid['FRP_uncertainty'] 
    m8_grid['M11_PIXEL_COUNT'] = m11_grid['pixel_count'].astype(int)
    m8_grid['FRP_diff'] = m8_grid['M11_summed_FRP'] - m8_grid['M8_summed_FRP']
    
    grid = m8_grid
    grid = grid.drop(columns=['summed_FRP','pixel_count','mean_vza','SAT','FRP_uncertainty'])
    """
    grid = pd.concat([m8_grid,m11_grid])
    return grid
    


def fit_and_plot_OLS_bins(fire_data):  

    results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_summed_FRP','M8_summed_FRP') 
    print results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual
    ####remove outliers and re-fit
    #fire_data = fire_data[~fire_data.iloc[:].index.isin(outlier_ix)]
    #results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_summed_FRP','M8_summed_FRP')
    #outlier_count_all = len(outlier_ix)
    #outlier_fires = fire_data[fire_data.iloc[:].index.isin(outlier_ix)]
    #plot_bin_maps(outlier_fires,fires_m8,'clusters')
    #######
    
    fire_max = max(fire_data['M11_summed_FRP'].max(),fire_data['M8_summed_FRP'].max())
    bias_all = fire_data['FRP_diff'].mean()
    SD_all = fire_data['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(fire_data['M11_summed_FRP'])
    
    fit_points = (np.array([0,max(fire_data['M11_summed_FRP'].values)]))
    perf_fit = (fit_points)*1 + 0
    z = (fit_points)*slope_all + intercept_all
    textstr1 = 'y = %.2fx + %.2f \n Fire count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope_all, intercept_all,fire_count_all, r2_all, std_err_all)
    textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias_all, SD_all, rmsd_all)
    plot_scatter_fit(fire_data['M11_summed_FRP'],fire_data['M8_summed_FRP'],z,perf_fit,fit_points,'all',('./Plots/FRP_ols'+plot_tag), textstr1,textstr2,fire_max)
    #plot_scatter_residuals(fire_data['M11_summed_FRP'],fire_data['M8_summed_FRP'],z, residual,slope_all,intercept_all,vza_max,fire_count_all,'all',('./Plots/FRP_OLS'+plot_tag))
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT':fire_count_all,'SLOPE':slope_all,'INTERCEPT':intercept_all, 'R2':r2_all,'SE':std_err_all})
    

    return stats_f
       
    

def fit_and_plot_ODR_bins(grid):
    print grid['M11_summed_FRP'],grid['M8_summed_FRP'],grid['M11_FRP_uncertainty'],grid['M8_FRP_uncertainty']
    
    params,param_uncertainty, adjusted_err, residual = odr_regression(grid['M11_summed_FRP'],grid['M8_summed_FRP'],grid['M11_FRP_uncertainty'],grid['M8_FRP_uncertainty'])
    print params,param_uncertainty, adjusted_err, residual
    fire_max = max(grid['M11_summed_FRP'].max(),grid['M8_summed_FRP'].max())
    bias_all = grid['FRP_diff'].mean()
    SD_all = grid['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(grid['M8_summed_FRP'])
    
    fit_points = (np.array([0,max(grid['M11_summed_FRP'].values)]))
    m8_fit = fit_points*params[1] + params[0]
    perf_fit = fit_points*1 + 0
    
    plot_ODR_scatter_residuals(grid['M11_summed_FRP'],grid['M8_summed_FRP'],grid['M11_FRP_uncertainty'],grid['M8_FRP_uncertainty'],m8_fit,perf_fit,fit_points, residual,adjusted_err,params, fire_max,fire_count_all,'all',('C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/Plots/AREA_FRP_ODR'+plot_tag))
        
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT': fire_count_all,'SLOPE':params[1] ,'INTERCEPT':params[0],'SE_slope':param_uncertainty[1],'SE_intercept':param_uncertainty[0]})

    return stats_f


#### Parameters ####
    
BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox"

resolution = 1.0
lonmin = 10.0
lonmax = 30.0
latmin = -10.0
latmax = 10.0
lon = np.arange(lonmin,(lonmax+resolution),resolution)
lat = np.arange(latmin,(latmax+resolution),resolution)
steps_lat = len(lat)
steps_long = len(lon)  
    
lon_min = np.tile(lon[:-1],steps_long)
lon_max = np.tile(lon[1:],steps_long)
lat_min = np.repeat(lat[:-1],steps_lat)
lat_max = np.repeat(lat[1:],steps_lat)
####################


    #m8_grid = pd.DataFrame({'MIN_LAT':lat_min,'MAX_LAT':lat_max,'MIN_LONG':lon_min,'MAX_LONG':lon_max,'GRID_NO': None,'summed_FRP': 0,'FRP_uncertainty':0,'mean_vza': 0, 'pixel_count':0 , 'SAT': 8})

dates = pd.date_range(start='2018-06-09',end='2018-07-03')
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
            grid = write_gridded_file(m8_df,m11_df)
            grid['DATE'] = date
            grid['TIME'] = time
            mode = 'w' if header_1 else 'a'
            grid.to_csv(BaseDir +"/txt_data/AREA_grid_"+date, mode=mode,header=header_1,index=False) 
            header_1=False
            
            
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