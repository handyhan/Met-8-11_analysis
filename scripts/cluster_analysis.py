# -*- coding: utf-8 -*-
"""
Created on Mon May 21 11:33:38 2018

@author: Hannah.N
"""
import os.path
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import sys
import math
from MSG_read import extract_data,plot_scatter_fit, plot_single_map,plot_ODR_scatter_residuals,plot_scatter_area,plot_scatter_log
from math import sin, cos, sqrt, atan2, radians
from regression_functions import regression_c0,regression_free,odr_regression,regression_log
pd.options.mode.chained_assignment = None  
pd.set_option('display.max_columns', 500)


BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox" 

def read_in_from_csv(dates,time_split,side):
    if os.path.isfile(BaseDir +"/txt_data/CLUSTERS_FD_fires_"+dates[0]):
        df_fires = pd.read_csv(BaseDir +"/txt_data/CLUSTERS_FD_fires_"+dates[0])
        df_pixels = pd.read_csv(BaseDir +"/txt_data/CLUSTERS_FD_pixels_"+dates[0])
        for date in dates[1:]:
            df_fires_day=pd.read_csv(BaseDir +"/txt_data/CLUSTERS_FD_fires_"+date)
            df_pixels_day = pd.read_csv(BaseDir +"/txt_data/CLUSTERS_FD_pixels_"+date)
            df_fires =  df_fires.append(df_fires_day,ignore_index=True) # = pd.concat([df_fires,df_fires_day],axis=0, join='outer')
            df_pixels = pd.concat([df_pixels,df_pixels_day],axis=0, join='outer')
        
        if time_split == True:
            fire_data = split_timea(df_fires, "daily_average")
            #fire_data = pd.DataFrame({'DATE':date.values[:],'TIME':time.values[:],'M8_summed_FRP':m8_frp.values[:],'M8_FRP_uncertainty':m8_frp_err.values[:],'M11_summed_FRP':m11_frp.values[:],'M11_FRP_uncertainty':m11_frp_err.values[:],'M8_PIXEL_COUNT':m8_pixel_count.values[:],'M11_PIXEL_COUNT':m11_pixel_count.values[:],'FRP_diff':fire_diff,'VZA_diff':vza_diff,'GRID_NO':clusters})
        
        elif time_split == False:
            vza_ratio,frp_ratio,lat,long, m8_frp,m8_frp_err,m11_frp,m11_frp_err,m8_pixel_count,m11_pixel_count,fire_diff,vza_diff,clusters = procces_df_fires(df_fires,side)    
            fire_data = pd.DataFrame({'VZA_RATIO':vza_ratio,'FRP_RATIO':frp_ratio,'LATITUDE':lat.values[:],'LONGITUDE':long.values[:],'M8_summed_FRP':m8_frp.values[:],'M8_FRP_uncertainty':m8_frp_err.values[:],'M11_summed_FRP':m11_frp.values[:],'M11_FRP_uncertainty':m11_frp_err.values[:],'M8_PIXEL_COUNT':m8_pixel_count.values[:],'M11_PIXEL_COUNT':m11_pixel_count.values[:],'FRP_diff':fire_diff,'VZA_diff':vza_diff,'cluster_no':clusters})
        
        
        return fire_data


def procces_df_fires(df_fires,side):        
    df_fires = df_fires.sort_values(by=['cluster_no','SAT'],ascending=True)        
    m8_frp = df_fires.groupby(['SAT']).get_group(8)['summed_FRP']
    m8_vza = df_fires.groupby(['SAT']).get_group(8)['mean_vza']
    m8_frp_err = df_fires.groupby(['SAT']).get_group(8)['FRP_uncertainty']
    m11_frp = df_fires.groupby(['SAT']).get_group(11)['summed_FRP']
    m11_vza = df_fires.groupby(['SAT']).get_group(11)['mean_vza']
    m11_frp_err = df_fires.groupby(['SAT']).get_group(11)['FRP_uncertainty']
    lat = df_fires.groupby(['SAT']).get_group(8)['mean_lat']
    long = df_fires.groupby(['SAT']).get_group(8)['mean_long']
    clusters = df_fires.groupby(['SAT']).get_group(8)['cluster_no']
    m11_pixel_count = df_fires.groupby(['SAT']).get_group(11)['pixel_count']
    m8_pixel_count = df_fires.groupby(['SAT']).get_group(8)['pixel_count']
    fires_m8 = df_fires[df_fires['SAT'] == 8 ].reset_index(drop=True)
    fires_m11 = df_fires[df_fires['SAT'] == 11 ].reset_index(drop=True)
    fire_diff = m8_frp.values[:] - m11_frp.values[:] 
    vza_diff = m8_vza.values[:] - m11_vza.values[:]                 
    if side == '_negative_':
        frp_ratio = (m11_frp.values[:])/(m8_frp.values[:])
        vza_ratio = (m11_vza.values[:])/(m8_vza.values[:])
    if side == '_positive_':
       frp_ratio = (m8_frp.values[:])/(m11_frp.values[:] )
       vza_ratio = (m8_vza.values[:])/(m11_vza.values[:] )
       
    return vza_ratio,frp_ratio,lat,long, m8_frp,m8_frp_err,m11_frp,m11_frp_err,m8_pixel_count,m11_pixel_count,fire_diff,vza_diff,clusters
        



def fit_and_plot_OLS_log(fire_data,BaseDir, x, y, scaling_factor):  
    m8_variable = y
    m11_variable = x
    results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, m11_variable,m8_variable) 
    ####remove outliers and re-fit
    #fire_data = fire_data[~fire_data.iloc[:].index.isin(outlier_ix)]
    #results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_summed_FRP','M8_summed_FRP')
    #outlier_count_all = len(outlier_ix)
    #outlier_fires = fire_data[fire_data.iloc[:].index.isin(outlier_ix)]
    #plot_bin_maps(outlier_fires,fires_m8,'clusters')
    #######
    vza_max = fire_data['VZA_diff'].max()
    fire_max = max(fire_data[m11_variable].max(),fire_data[m8_variable].max())
    fire_max  = int(math.ceil(fire_max / scaling_factor)) * scaling_factor
    bias_all = fire_data['FRP_diff'].mean()
    SD_all = fire_data['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(fire_data['M11_summed_FRP'])
    fit_points = (np.array([0,fire_max]))
    perf_fit = (fit_points)*1 + 0
    z = (fit_points)*slope_all + intercept_all
    textstr1 = 'y = %.2fx + %.2f \n Fire count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope_all, intercept_all,fire_count_all, r2_all, std_err_all)
    textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias_all, SD_all, rmsd_all)
    #plot_scatter_area(True,fire_data[m11_variable],fire_data[m8_variable],fire_data['VZA_diff'],z,perf_fit,fit_points,'all',(BaseDir+'/Plots/'+var_to_plot+'_colour_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, var_to_plot)
    plot_scatter_area(False,fire_data[m11_variable],fire_data[m8_variable],fire_data['VZA_diff'],z,perf_fit,fit_points,'all',(BaseDir+'/Plots/'+ x +'_OLS_LOG'+plot_tag), textstr1,textstr2,vza_max,fire_max, x,y,scaling_factor)
    
    vza_bins =['all_deg','deg_10','deg10_20','deg20_30','deg30_40','deg40_']
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT':fire_count_all,'SLOPE':slope_all,'INTERCEPT':intercept_all, 'R2':r2_all,'SE':std_err_all})
    
    for bin in vza_bins[1:]:
        actual_bins = np.unique(fire_data['vza_bin'])
        if bin in actual_bins:    
            m8_frp = fire_data.groupby(['vza_bin']).get_group(bin)[m8_variable]
            m11_frp = fire_data.groupby(['vza_bin']).get_group(bin)[m11_variable]
            vza_diff = fire_data.groupby(['vza_bin']).get_group(bin)['VZA_diff']  
            group_fire = pd.DataFrame({m11_variable :m11_frp, m8_variable:m8_frp,'VZA_diff':vza_diff})
            if len(group_fire[m8_variable]) > 1 :
                results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, m11_variable,m8_variable)
                ###remove outliers and re-fit
                #group_fire = group_fire[~group_fire.iloc[:].index.isin(outlier_ix)]
                #results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, m11_variable,'M8_summed_FRP')
                #outlier_count = len(outlier_ix)
                ######
                        
                fire_count = len(group_fire[m8_variable])
                bias = fire_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].mean()
                SD = fire_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].std()
                rmsd = sqrt(bias**2 + SD**2)
                stat_title = ['BIAS','SCATTER','RMSD','FIRE_COUNT','SLOPE','INTERCEPT','R2','SE'] ; stat_value = [bias,SD,rmsd,fire_count,slope,intercept,r2,std_err]
                for i in range(0,8,1):
                    n = stats_f.VZA_BIN[stats_f.VZA_BIN == bin].index            
                    stats_f.at[n,stat_title[i]] = stat_value[i]
                    
                perf_fit = (fit_points)*1 + 0
                z = (fit_points)*slope + intercept
                
                textstr1 ='y = %.2fx + %.2f \n Fire count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope, intercept,fire_count, r2,std_err)
                textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias, SD, rmsd)
                #plot_scatter_area(True,group_fire[m11_variable],group_fire[m8_variable],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,(BaseDir+'/Plots/'+var_to_plot+'_colour_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, var_to_plot)
                plot_scatter_area(False,group_fire[m11_variable],group_fire[m8_variable],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,(BaseDir+'/Plots/'+ x +'_OLS_LOG'+plot_tag), textstr1,textstr2,vza_max,fire_max, x,y,scaling_factor)
    
     
def fit_and_plot_ODR_bins(fire_data,BaseDir):  
    params,param_uncertainty, adjusted_err, residual = odr_regression(fire_data['M11_summed_FRP'],fire_data['M8_summed_FRP'],fire_data['M11_FRP_uncertainty'],fire_data['M8_FRP_uncertainty'])
    vza_max = fire_data['VZA_diff'].max()    
    fire_max = max(fire_data['M11_summed_FRP'].max(),fire_data['M8_summed_FRP'].max())
    fire_max  = int(math.ceil(fire_max / 500.0)) * 500
    bias_all = fire_data['FRP_diff'].mean()
    SD_all = fire_data['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(fire_data['M8_summed_FRP'])
    
    fit_points = (np.array([0,fire_max]))
    m8_fit = fit_points*params[1] + params[0]
    perf_fit = fit_points*1 + 0
    textstr1 ='y = %.2fx + %.2f \n Fire count=%d '%(params[1],  params[0],fire_count_all)
    plot_ODR_scatter_residuals(fire_data['M11_summed_FRP'],fire_data['M8_summed_FRP'],fire_data['M11_FRP_uncertainty'],fire_data['M8_FRP_uncertainty'],m8_fit,perf_fit,fit_points, residual,adjusted_err,params,vza_max, fire_max,fire_count_all,'all',(BaseDir+'/Plots/FRP_ODR'+plot_tag),textstr1)
        
    vza_bins =['all_deg','deg_10','deg10_20','deg20_30','deg30_40','deg40_']
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT': fire_count_all,'SLOPE':params[1] ,'INTERCEPT':params[0],'SE_slope':param_uncertainty[1],'SE_intercept':param_uncertainty[0]})
    for bin in vza_bins[1:]:
        actual_bins = np.unique(fire_data['vza_bin'])
        if bin in actual_bins:    
            m8_frp = fire_data.groupby(['vza_bin']).get_group(bin)['M8_summed_FRP']
            m11_frp = fire_data.groupby(['vza_bin']).get_group(bin)['M11_summed_FRP']
            m8_err = fire_data.groupby(['vza_bin']).get_group(bin)['M8_FRP_uncertainty']
            m11_err = fire_data.groupby(['vza_bin']).get_group(bin)['M11_FRP_uncertainty']
            vza_diff = fire_data.groupby(['vza_bin']).get_group(bin)['VZA_diff']  
            vza_max = vza_diff.max()
            
            params,param_uncertainty, adjusted_err, residual = odr_regression(m11_frp,m8_frp,m11_err,m8_err)
                    
            fire_count = len(m11_err)
            bias = fire_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].mean()
            SD = fire_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].std()
            rmsd = sqrt(bias**2 + SD**2)
            stat_title = ['BIAS','SCATTER','RMSD','FIRE_COUNT','SLOPE','INTERCEPT','SE_slope','SE_intercept'] ; stat_value = [bias,SD,rmsd, fire_count ,params[1],params[0],param_uncertainty[1],param_uncertainty[0]]
            for i in range(0,8,1):
                n = stats_f.VZA_BIN[stats_f.VZA_BIN == bin].index            
                stats_f.at[n,stat_title[i]] = stat_value[i]
            
            m8_fit = fit_points*params[1] + params[0]
            perf_fit = fit_points*1 + 0
            textstr1 ='y = %.2fx + %.2f \n Fire count=%d '%(params[1],  params[0],fire_count)
    #textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias, SD, rmsd)
            plot_ODR_scatter_residuals(m11_frp,m8_frp,m11_err,m8_err,m8_fit,perf_fit,fit_points, residual,adjusted_err,params,vza_max,fire_max,fire_count,bin,(BaseDir+'/Plots/FRP_ODR'+plot_tag),textstr1)
         
    return stats_f

def fit_and_plot_OLS_bins(fire_data,BaseDir, x, y, scaling_factor):  
    m8_variable = y
    m11_variable = x
    results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, m11_variable,m8_variable) 
    ####remove outliers and re-fit
    #fire_data = fire_data[~fire_data.iloc[:].index.isin(outlier_ix)]
    #results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_summed_FRP','M8_summed_FRP')
    #outlier_count_all = len(outlier_ix)
    #outlier_fires = fire_data[fire_data.iloc[:].index.isin(outlier_ix)]
    #plot_bin_maps(outlier_fires,fires_m8,'clusters')
    #######
    vza_max = fire_data['VZA_diff'].max()
    fire_max = max(fire_data[m11_variable].max(),fire_data[m8_variable].max())
    fire_max  = (math.ceil(fire_max / scaling_factor) * scaling_factor) + 10
    
    bias_all = fire_data['FRP_diff'].mean()
    SD_all = fire_data['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(fire_data['M11_summed_FRP'])
    fit_points = (np.array([0,fire_max]))
    perf_fit = (fit_points)*1 + 0
    z = (fit_points)*slope_all + intercept_all
    
    textstr1 = 'y = %.2fx + %.2f \n Fire count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope_all, intercept_all,fire_count_all, r2_all, std_err_all)
    textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias_all, SD_all, rmsd_all)
    #plot_scatter_area(True,fire_data[m11_variable],fire_data[m8_variable],fire_data['VZA_diff'],z,perf_fit,fit_points,'all',(BaseDir+'/Plots/'+var_to_plot+'_colour_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, var_to_plot)
    plot_scatter_area(False,fire_data[m11_variable],fire_data[m8_variable],fire_data['VZA_diff'],z,perf_fit,fit_points,'all',(BaseDir+'/Plots/'+ x +'_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, x,y,scaling_factor)
    
    vza_bins =['all_deg','deg_10','deg10_20','deg20_30','deg30_40','deg40_']
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT':fire_count_all,'SLOPE':slope_all,'INTERCEPT':intercept_all, 'R2':r2_all,'SE':std_err_all})
    
    for bin in vza_bins[1:]:
        actual_bins = np.unique(fire_data['vza_bin'])
        if bin in actual_bins:    
            m8_frp = fire_data.groupby(['vza_bin']).get_group(bin)[m8_variable]
            m11_frp = fire_data.groupby(['vza_bin']).get_group(bin)[m11_variable]
            vza_diff = fire_data.groupby(['vza_bin']).get_group(bin)['VZA_diff']  
            group_fire = pd.DataFrame({m11_variable :m11_frp, m8_variable:m8_frp,'VZA_diff':vza_diff})
            if len(group_fire[m8_variable]) > 1 :
                results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, m11_variable,m8_variable)
                ###remove outliers and re-fit
                #group_fire = group_fire[~group_fire.iloc[:].index.isin(outlier_ix)]
                #results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, m11_variable,'M8_summed_FRP')
                #outlier_count = len(outlier_ix)
                ######
                        
                fire_count = len(group_fire[m8_variable])
                bias = fire_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].mean()
                SD = fire_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].std()
                rmsd = sqrt(bias**2 + SD**2)
                stat_title = ['BIAS','SCATTER','RMSD','FIRE_COUNT','SLOPE','INTERCEPT','R2','SE'] ; stat_value = [bias,SD,rmsd,fire_count,slope,intercept,r2,std_err]
                for i in range(0,8,1):
                    n = stats_f.VZA_BIN[stats_f.VZA_BIN == bin].index            
                    stats_f.at[n,stat_title[i]] = stat_value[i]
                    
                perf_fit = (fit_points)*1 + 0
                z = (fit_points)*slope + intercept
                
                textstr1 ='y = %.2fx + %.2f \n Fire count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope, intercept,fire_count, r2,std_err)
                textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias, SD, rmsd)
                #plot_scatter_area(True,group_fire[m11_variable],group_fire[m8_variable],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,(BaseDir+'/Plots/'+var_to_plot+'_colour_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, var_to_plot)
                plot_scatter_area(False,group_fire[m11_variable],group_fire[m8_variable],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,(BaseDir+'/Plots/'+ x +'_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, x,y,scaling_factor)
    
     

def plot_bin_maps(fire_data,fires_m8,title):
    
    #fires_m8= fires_m8[fires_m8['mean_long'] > 40]    
    #fires_m8= fires_m8[fires_m8['mean_long'] < 60]
    #fires_m8= fires_m8[fires_m8['mean_lat'] > -20]    
    #fires_m8= fires_m8[fires_m8['mean_lat'] < -10]
    #list_fires = np.unique(fires_m8['cluster_no'])
    #select_fires = pixels_m8[pixels_m8['cluster_no'].isin(list_fires)]
    #print pixels_m8 
    #vza_bins =['deg_10','deg10_20','deg20_30','deg30_40','deg40_']    
    #for bin in vza_bins:
    #    if bin in fire_data['vza_bin'].unique():
    #cluster_list = fire_data.groupby(['vza_bin']).get_group(bin)['cluster_no']   
    #cluster_bin = fire_data[fire_data['cluster_no'].isin(cluster_list)]
    latitude = fires_m8['LATITUDE'].values[:]
    longitude = fires_m8['LONGITUDE'].values[:]
    #latitude = fires_m8['mean_lat'].values[:]# fires_m8[fires_m8['cluster_no'].isin(cluster_list)]['mean_lat'].values[:]
    #longitude = fires_m8['mean_long'].values[:] #fires_m8[fires_m8['cluster_no'].isin(cluster_list)]['mean_long'].values[:]
    plot_single_map(latitude,longitude,(title +plot_tag))
  
    

#dates = pd.date_range(start='2018-06-09',end='2018-06-20')

#dates = dates.format(formatter=lambda x: x.strftime('%Y%m%d'))
dates = ['20180610','20180611','20180612','20180613','20180616','20180617','20180618','20180619','20180621','20180624','20180625','20180626','20180630','20180701','20180702','20180704','20180705','20180707']

side = '_positive_'
plot_tag =''

fire_data = read_in_from_csv(dates,False,side)

#fire_data = fire_data[(fire_data['LONGITUDE'] > 5) & (fire_data['LONGITUDE'] < 25)]

#fire_data = fire_data[(fire_data['LATITUDE'] < 5) & (fire_data['LATITUDE'] > -5)]
#print fire_data[fire_data['TIME'] == 1315 ]

fire_data = fire_data[(fire_data['M11_summed_FRP'] > 50) & (fire_data['M8_summed_FRP'] > 50)]


if side == '_negative_':
    fire_data = fire_data[fire_data['VZA_diff'] < 0 ]
    plot_tag = side + plot_tag
if side == '_positive_':
    fire_data = fire_data[fire_data['VZA_diff'] > 0 ]
    plot_tag = side + plot_tag

fire_data['VZA_diff'] = fire_data['VZA_diff'].abs()

fire_data['vza_bin'] = np.where((fire_data.VZA_diff <= 10),
              'deg_10',
              np.where(np.logical_and(fire_data.VZA_diff > 10,fire_data.VZA_diff <= 20),            
              'deg10_20',
              np.where(np.logical_and(fire_data.VZA_diff > 20,fire_data.VZA_diff <= 30),
              'deg20_30',
              np.where(np.logical_and(fire_data.VZA_diff > 30,fire_data.VZA_diff <= 40),
              'deg30_40',
              'deg40_'))))

#list_fires = np.unique(fire_data['cluster_no'])
#select_fires = pixels_m8[pixels_m8['cluster_no'].isin(list_fires)]
#plot_bin_maps(fire_data,select_fires,'TEST_fire_locations')
  
fit_and_plot_ODR_bins(fire_data,BaseDir)
#plot_single_map(fire_data['LATITUDE'].values[:], fire_data['LONGITUDE'].values[:],fire_data['VZA_diff'].values[:], ('FIRES'+plot_tag))
#fire_data['LOG_VZA_RATIO'] = np.log(fire_data['VZA_RATIO'])
#fire_data['LOG_FRP_RATIO'] = np.log(fire_data['FRP_RATIO'])

fit_and_plot_OLS_bins(fire_data,BaseDir,'M11_PIXEL_COUNT', 'M8_PIXEL_COUNT',5.0 )
#labels = ['Bias','Scatter','RMSD']
#plot_stats(stats_f['VZA_BIN'][:],stats_f['BIAS'][:],stats_f['SCATTER'][:],stats_f['RMSD'][:],'./Plots/Bias_scatter_rmsd',plot_tag,labels)
#labels = ['Slope','R$^2$','Standard Error']
#plot_stats(stats_f['VZA_BIN'][:],stats_f['SLOPE'][:],stats_f['R2'][:],stats_f['SE'][:],'./Plots/LinReg_params',plot_tag,labels)
#stats_f.to_csv('./fire_stats_table' +plot_tag +'.csv')



################################ runnint cose



  
"""
def fit_and_plot_OLS_bins_pixel(fire_data):  

    results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_PIXEL_COUNT','M8_PIXEL_COUNT') 
    #remove outliers and re-fit
    #fire_data = fire_data[~fire_data.iloc[:].index.isin(outlier_ix)]
    #results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_PIXEL_COUNT','M8_PIXEL_COUNT')
    #outlier_count_all = len(outlier_ix)
    #outlier_fires = fire_data[fire_data.iloc[:].index.isin(outlier_ix)]
    #plot_bin_maps(outlier_fires,fires_m8,'clusters')
    #######
    
    vza_max = fire_data['VZA_diff'].max()
    fire_max = max(fire_data['M11_PIXEL_COUNT'].max(),fire_data['M8_PIXEL_COUNT'].max())
    print fire_max
    bias_all = fire_data['FRP_diff'].mean()
    SD_all = fire_data['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(fire_data['M11_PIXEL_COUNT'])
    
    fit_points = (np.array([0,max(fire_data['M11_PIXEL_COUNT'].values)]))
    perf_fit = (fit_points)*1 + 0
    z = (fit_points)*slope_all + intercept_all
    textstr1 = 'Fire count=%d \n slope=%.3f \n intercept=%.3f\n R$^{2}$=%.3f\n SE=%.3f'%(fire_count_all,slope_all, intercept_all, r2_all, std_err_all)
    textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias_all, SD_all, rmsd_all)
    plot_scatter_fit(fire_data['M11_PIXEL_COUNT'],fire_data['M8_PIXEL_COUNT'],fire_data['VZA_diff'],z,perf_fit,fit_points,'all',('./Plots/FRP_ols'+plot_tag), textstr1,textstr2,vza_max,fire_max)
    #plot_scatter_residuals(fire_data['M11_summed_FRP'],fire_data['M8_summed_FRP'],z, residual,slope_all,intercept_all,vza_max,fire_count_all,'all',('./Plots/FRP_OLS'+plot_tag))

    
    vza_bins =['all_deg','deg_10','deg10_20','deg20_30','deg30_40','deg40_']
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT':fire_count_all,'SLOPE':slope_all,'INTERCEPT':intercept_all, 'R2':r2_all,'SE':std_err_all})
    
    for bin in vza_bins[1:]:
            
        m8_frp = fire_data.groupby(['vza_bin']).get_group(bin)['M8_PIXEL_COUNT']
        m11_frp = fire_data.groupby(['vza_bin']).get_group(bin)['M11_PIXEL_COUNT']
        vza_diff = fire_data.groupby(['vza_bin']).get_group(bin)['VZA_diff']  
        group_fire = pd.DataFrame({'M11_PIXEL_COUNT':m11_frp,'M8_PIXEL_COUNT':m8_frp,'VZA_diff':vza_diff})
        results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, 'M11_PIXEL_COUNT','M8_PIXEL_COUNT')
        #remove outliers and re-fit
        #group_fire = group_fire[~group_fire.iloc[:].index.isin(outlier_ix)]
        #results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, 'M11_PIXEL_COUNT','M8_PIXEL_COUNT')
        #outlier_count = len(outlier_ix)
        ######
                
        fire_count = len(group_fire['M8_PIXEL_COUNT'])
        bias = fire_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].mean()
        SD = fire_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].std()
        rmsd = sqrt(bias**2 + SD**2)
        stat_title = ['BIAS','SCATTER','RMSD','FIRE_COUNT','SLOPE','INTERCEPT','R2','SE'] ; stat_value = [bias,SD,rmsd,fire_count,slope,intercept,r2,std_err]
        for i in range(0,8,1):
            n = stats_f.VZA_BIN[stats_f.VZA_BIN == bin].index            
            stats_f.at[n,stat_title[i]] = stat_value[i]
            
        perf_fit = (fit_points)*1 + 0
        z = (fit_points)*slope + intercept
        textstr1 ='Fire count=%d \n slope=%.3f \n intercept=%.3f\n R$^{2}$=%.3f\n SE=%.3f'%(fire_count,slope, intercept, r2,std_err)
        textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias, SD, rmsd)
        plot_scatter_fit(group_fire['M11_PIXEL_COUNT'],group_fire['M8_PIXEL_COUNT'],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,('./Plots/PIXEL_ols'+plot_tag), textstr1,textstr2,vza_max,fire_max)
        #plot_scatter_residuals(group_fire['M11_summed_FRP'],group_fire['M8_summed_FRP'],z,residual,slope,intercept,vza_max,fire_count,bin,('./Plots/FRP_OLS'+plot_tag))

    return stats_f
"""
    
    
    
    
    
    
    
    