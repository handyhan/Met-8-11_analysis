import os.path
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import sys
import math
from MSG_read import extract_data,plot_scatter_area, plot_single_map, plot_ODR_scatter_residuals
from regression_functions import regression_c0,regression_free,odr_regression
from math import sin, cos, sqrt, atan2, radians
from pandas.plotting import table


BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox" 

def split_timea(df_grid, analysis_type):  # analysis type can have values daily_average or chronological
    
    if  analysis_type == 'daily_average':
        
        df_grid = df_grid.sort_values(['TIME','GRID_NO','SAT'])
        times = df_grid['TIME'].unique()
        df_time = df_grid[df_grid['TIME'] == times[0] ]
        vza_ratio,frp_ratio,lat,long, m8_frp,m8_frp_err,m11_frp,m11_frp_err,m8_pixel_count,m11_pixel_count,fire_diff,vza_diff,clusters = procces_df_grid(df_time,side)    
        grid_data = pd.DataFrame({'TIME':times[0],'VZA_RATIO':vza_ratio,'FRP_RATIO':frp_ratio,'LATITUDE':lat.values[:],'LONGITUDE':long.values[:],'M8_summed_FRP':m8_frp.values[:],'M8_FRP_uncertainty':m8_frp_err.values[:],'M11_summed_FRP':m11_frp.values[:],'M11_FRP_uncertainty':m11_frp_err.values[:],'M8_PIXEL_COUNT':m8_pixel_count.values[:],'M11_PIXEL_COUNT':m11_pixel_count.values[:],'FRP_diff':fire_diff,'VZA_diff':vza_diff,'GRID_NO':clusters})
        
        for time in times[1:]:
            df_time = df_grid[df_grid['TIME'] == time ]
            vza_ratio,frp_ratio,lat,long, m8_frp,m8_frp_err,m11_frp,m11_frp_err,m8_pixel_count,m11_pixel_count,fire_diff,vza_diff,clusters    = procces_df_grid(df_time,side)    
            grid_data_1 = pd.DataFrame({'TIME':time,'VZA_RATIO':vza_ratio,'FRP_RATIO':frp_ratio,'LATITUDE':lat.values[:],'LONGITUDE':long.values[:],'M8_summed_FRP':m8_frp.values[:],'M8_FRP_uncertainty':m8_frp_err.values[:],'M11_summed_FRP':m11_frp.values[:],'M11_FRP_uncertainty':m11_frp_err.values[:],'M8_PIXEL_COUNT':m8_pixel_count.values[:],'M11_PIXEL_COUNT':m11_pixel_count.values[:],'FRP_diff':fire_diff,'VZA_diff':vza_diff,'GRID_NO':clusters})
            grid_data = pd.concat([grid_data,grid_data_1])
            
        #print grid_data
        grid_data = grid_data.sort_values(['TIME'])
        
        time_agg_m11 = grid_data.groupby(['TIME'])['M11_summed_FRP'].agg('sum')
        time_agg_m8 = grid_data.groupby(['TIME'])['M8_summed_FRP'].agg('sum')
        
        
        return grid_data
        
    elif analysis_type == 'time_series':
        pass
    else:
        print "ERROR: invalid analysis type"
        pass
    
def procces_df_grid(df_grid,side):
    
        df_grid = df_grid.sort_values(by=['GRID_NO','SAT'],ascending=True)        
        m8_frp = df_grid.groupby(['SAT']).get_group(8)['summed_FRP']
        m8_vza = df_grid.groupby(['SAT']).get_group(8)['VZA']
        m8_frp_err = df_grid.groupby(['SAT']).get_group(8)['FRP_uncertainty']
        m11_frp = df_grid.groupby(['SAT']).get_group(11)['summed_FRP']
        m11_vza = df_grid.groupby(['SAT']).get_group(11)['VZA']
        m11_frp_err = df_grid.groupby(['SAT']).get_group(11)['FRP_uncertainty']
        lat = df_grid.groupby(['SAT']).get_group(8)['LAT_CENTER']
        long = df_grid.groupby(['SAT']).get_group(8)['LONG_CENTER']
        clusters = df_grid.groupby(['SAT']).get_group(8)['GRID_NO']
        m11_pixel_count = df_grid.groupby(['SAT']).get_group(11)['pixel_count']
        m8_pixel_count = df_grid.groupby(['SAT']).get_group(8)['pixel_count']
        fires_m8 = df_grid[df_grid['SAT'] == 8 ].reset_index(drop=True)
        fires_m11 = df_grid[df_grid['SAT'] == 11 ].reset_index(drop=True)
        fire_diff = m8_frp.values[:] - m11_frp.values[:] 
        vza_diff = m8_vza.values[:] - m11_vza.values[:]                 
        if side == '_negative_':
            frp_ratio = (m11_frp.values[:])/(m8_frp.values[:])
            vza_ratio = (m11_vza.values[:])/(m8_vza.values[:])
        if side == '_positive_':
            frp_ratio = (m8_frp.values[:])/(m11_frp.values[:] )
            vza_ratio = (m8_vza.values[:])/(m11_vza.values[:] )
        
        return vza_ratio,frp_ratio,lat,long, m8_frp,m8_frp_err,m11_frp,m11_frp_err,m8_pixel_count,m11_pixel_count,fire_diff,vza_diff,clusters
        
        
def read_in_from_csv(dates,time_split,side):
    if os.path.isfile("C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/txt_data/AREA_grid_"+dates[0]):
        df_grid = pd.read_csv("C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/txt_data/AREA_grid_"+dates[0])
        for date in dates[1:]:
            df_grid_day=pd.read_csv("C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/txt_data/AREA_grid_"+date)
            df_grid =  df_grid.append(df_grid_day,ignore_index=True)
        
        
        
        
        if time_split == True:
            grid_data = split_timea(df_grid, "daily_average")
            #grid_data = pd.DataFrame({'DATE':date.values[:],'TIME':time.values[:],'M8_summed_FRP':m8_frp.values[:],'M8_FRP_uncertainty':m8_frp_err.values[:],'M11_summed_FRP':m11_frp.values[:],'M11_FRP_uncertainty':m11_frp_err.values[:],'M8_PIXEL_COUNT':m8_pixel_count.values[:],'M11_PIXEL_COUNT':m11_pixel_count.values[:],'FRP_diff':fire_diff,'VZA_diff':vza_diff,'GRID_NO':clusters})
        
        elif time_split == False:
            vza_ratio,frp_ratio,lat,long, m8_frp,m8_frp_err,m11_frp,m11_frp_err,m8_pixel_count,m11_pixel_count,fire_diff,vza_diff,clusters = procces_df_grid(df_grid,side)    
            grid_data = pd.DataFrame({'VZA_RATIO':vza_ratio,'FRP_RATIO':frp_ratio,'LATITUDE':lat.values[:],'LONGITUDE':long.values[:],'M8_summed_FRP':m8_frp.values[:],'M8_FRP_uncertainty':m8_frp_err.values[:],'M11_summed_FRP':m11_frp.values[:],'M11_FRP_uncertainty':m11_frp_err.values[:],'M8_PIXEL_COUNT':m8_pixel_count.values[:],'M11_PIXEL_COUNT':m11_pixel_count.values[:],'FRP_diff':fire_diff,'VZA_diff':vza_diff,'GRID_NO':clusters})
        
        
        return grid_data





def fit_and_plot_OLS_bins(grid_data,BaseDir, x, y, scaling_factor):  
    m8_variable = y
    m11_variable = x
    results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(grid_data, m11_variable,m8_variable) 
    ####remove outliers and re-fit
    #fire_data = fire_data[~fire_data.iloc[:].index.isin(outlier_ix)]
    #results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_summed_FRP','M8_summed_FRP')
    #outlier_count_all = len(outlier_ix)
    #outlier_fires = fire_data[fire_data.iloc[:].index.isin(outlier_ix)]
    #plot_bin_maps(outlier_fires,fires_m8,'clusters')
    #######
    vza_max = grid_data['VZA_diff'].max()
    fire_max = max(grid_data[m11_variable].max(),grid_data[m8_variable].max())
    fire_max  = int(math.ceil(fire_max / scaling_factor)) * scaling_factor
    bias_all = grid_data['FRP_diff'].mean()
    SD_all = grid_data['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(grid_data['M11_summed_FRP'])
    fit_points = (np.array([0,fire_max]))
    perf_fit = (fit_points)*1 + 0
    z = (fit_points)*slope_all + intercept_all
    textstr1 = 'y = %.2fx + %.2f \n Grid-cell count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope_all, intercept_all,fire_count_all, r2_all, std_err_all)
    textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias_all, SD_all, rmsd_all)
    #plot_scatter_area(True,grid_data[m11_variable],grid_data[m8_variable],grid_data['VZA_diff'],z,perf_fit,fit_points,'all',(BaseDir+'/Plots/'+var_to_plot+'_colour_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, var_to_plot)
    plot_scatter_area(False,grid_data[m11_variable],grid_data[m8_variable],grid_data['VZA_diff'],z,perf_fit,fit_points,'all',(BaseDir+'/Plots/'+ x +'_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, x,y,scaling_factor)
    
    vza_bins =['all_deg','deg_10','deg10_20','deg20_30','deg30_40','deg40_']
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT':fire_count_all,'SLOPE':slope_all,'INTERCEPT':intercept_all, 'R2':r2_all,'SE':std_err_all})
    
    for bin in vza_bins[1:]:
        actual_bins = np.unique(grid_data['vza_bin'])
        if bin in actual_bins:    
            m8_frp = grid_data.groupby(['vza_bin']).get_group(bin)[m8_variable]
            m11_frp = grid_data.groupby(['vza_bin']).get_group(bin)[m11_variable]
            vza_diff = grid_data.groupby(['vza_bin']).get_group(bin)['VZA_diff']  
            group_fire = pd.DataFrame({m11_variable :m11_frp, m8_variable:m8_frp,'VZA_diff':vza_diff})
            if len(group_fire[m8_variable]) > 1 :
                results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, m11_variable,m8_variable)
                ###remove outliers and re-fit
                #group_fire = group_fire[~group_fire.iloc[:].index.isin(outlier_ix)]
                #results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, m11_variable,'M8_summed_FRP')
                #outlier_count = len(outlier_ix)
                ######
                        
                fire_count = len(group_fire[m8_variable])
                bias = grid_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].mean()
                SD = grid_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].std()
                rmsd = sqrt(bias**2 + SD**2)
                stat_title = ['BIAS','SCATTER','RMSD','FIRE_COUNT','SLOPE','INTERCEPT','R2','SE'] ; stat_value = [bias,SD,rmsd,fire_count,slope,intercept,r2,std_err]
                for i in range(0,8,1):
                    n = stats_f.VZA_BIN[stats_f.VZA_BIN == bin].index            
                    stats_f.at[n,stat_title[i]] = stat_value[i]
                    
                perf_fit = (fit_points)*1 + 0
                z = (fit_points)*slope + intercept
                
                textstr1 ='y = %.2fx + %.2f \n Grid-cell count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope, intercept,fire_count, r2,std_err)
                textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias, SD, rmsd)
                #plot_scatter_area(True,group_fire[m11_variable],group_fire[m8_variable],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,(BaseDir+'/Plots/'+var_to_plot+'_colour_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, var_to_plot)
                plot_scatter_area(False,group_fire[m11_variable],group_fire[m8_variable],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,(BaseDir+'/Plots/'+ x +'_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, x,y,scaling_factor)
    
    #plot_scatter_residuals(group_fire['M11_summed_FRP'],group_fire['M8_summed_FRP'],z,residual,slope,intercept,vza_max,fire_count,bin,('./Plots/FRP_OLS'+plot_tag))
     
    
    #return stats_f
    
def fit_and_plot_OLS_log(grid_data,BaseDir, x, y, scaling_factor): 
    print grid_data
    m8_variable = y
    m11_variable = x
    results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(grid_data, m11_variable,m8_variable) 
    ####remove outliers and re-fit
    #fire_data = fire_data[~fire_data.iloc[:].index.isin(outlier_ix)]
    #results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_summed_FRP','M8_summed_FRP')
    #outlier_count_all = len(outlier_ix)
    #outlier_fires = fire_data[fire_data.iloc[:].index.isin(outlier_ix)]
    #plot_bin_maps(outlier_fires,fires_m8,'clusters')
    #######
    print grid_data, "hello1"
    vza_max = grid_data['VZA_diff'].max()
    fire_max = max(grid_data[m11_variable].max(),grid_data[m8_variable].max())
    fire_max  = int(math.ceil(fire_max / scaling_factor)) * scaling_factor
    bias_all = grid_data['FRP_diff'].mean()
    SD_all = grid_data['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(grid_data['M11_summed_FRP'])
    fit_points = (np.array([0,fire_max]))
    perf_fit = (fit_points)*1 + 0
    z = (fit_points)*slope_all + intercept_all
    print grid_data, "hello2"
    textstr1 = 'y = %.2fx + %.2f \n Grid-cell count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope_all, intercept_all,fire_count_all, r2_all, std_err_all)
    textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias_all, SD_all, rmsd_all)
    #plot_scatter_area(True,grid_data[m11_variable],grid_data[m8_variable],grid_data['VZA_diff'],z,perf_fit,fit_points,'all',(BaseDir+'/Plots/'+var_to_plot+'_colour_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, var_to_plot)
    plot_scatter_area(False,grid_data[m11_variable],grid_data[m8_variable],grid_data['VZA_diff'],z,perf_fit,fit_points,'all',(BaseDir+'/Plots/'+ x +'_OLS_LOG'+plot_tag), textstr1,textstr2,vza_max,fire_max, x,y,scaling_factor)
    
    vza_bins =['all_deg','deg_10','deg10_20','deg20_30','deg30_40','deg40_']
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT':fire_count_all,'SLOPE':slope_all,'INTERCEPT':intercept_all, 'R2':r2_all,'SE':std_err_all})
    
    for bin in vza_bins[1:]:
        actual_bins = np.unique(grid_data['vza_bin'])
        if bin in actual_bins:    
            m8_frp = grid_data.groupby(['vza_bin']).get_group(bin)[m8_variable]
            m11_frp = grid_data.groupby(['vza_bin']).get_group(bin)[m11_variable]
            vza_diff = grid_data.groupby(['vza_bin']).get_group(bin)['VZA_diff']  
            group_fire = pd.DataFrame({m11_variable :m11_frp, m8_variable:m8_frp,'VZA_diff':vza_diff})
            if len(group_fire[m8_variable]) > 1 :
                results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, m11_variable,m8_variable)
                ###remove outliers and re-fit
                #group_fire = group_fire[~group_fire.iloc[:].index.isin(outlier_ix)]
                #results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, m11_variable,'M8_summed_FRP')
                #outlier_count = len(outlier_ix)
                ######
                        
                fire_count = len(group_fire[m8_variable])
                bias = grid_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].mean()
                SD = grid_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].std()
                rmsd = sqrt(bias**2 + SD**2)
                stat_title = ['BIAS','SCATTER','RMSD','FIRE_COUNT','SLOPE','INTERCEPT','R2','SE'] ; stat_value = [bias,SD,rmsd,fire_count,slope,intercept,r2,std_err]
                for i in range(0,8,1):
                    n = stats_f.VZA_BIN[stats_f.VZA_BIN == bin].index            
                    stats_f.at[n,stat_title[i]] = stat_value[i]
                    
                perf_fit = (fit_points)*1 + 0
                z = (fit_points)*slope + intercept
                
                textstr1 ='y = %.2fx + %.2f \n Grid-cell count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope, intercept,fire_count, r2,std_err)
                textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias, SD, rmsd)
                #plot_scatter_area(True,group_fire[m11_variable],group_fire[m8_variable],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,(BaseDir+'/Plots/'+var_to_plot+'_colour_OLS'+plot_tag), textstr1,textstr2,vza_max,fire_max, var_to_plot)
                plot_scatter_area(False,group_fire[m11_variable],group_fire[m8_variable],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,(BaseDir+'/Plots/'+ x +'_OLS_LOG'+plot_tag), textstr1,textstr2,vza_max,fire_max, x,y,scaling_factor)
    
    
def fit_and_plot_ODR_bins(grid_data,BaseDir):  
    params,param_uncertainty, adjusted_err, residual = odr_regression(grid_data['M11_summed_FRP'],grid_data['M8_summed_FRP'],grid_data['M11_FRP_uncertainty'],grid_data['M8_FRP_uncertainty'])
    vza_max = grid_data['VZA_diff'].max()    
    fire_max = max(grid_data['M11_summed_FRP'].max(),grid_data['M8_summed_FRP'].max())
    fire_max  = int(math.ceil(fire_max / 2000.0)) * 2000
    bias_all = grid_data['FRP_diff'].mean()
    SD_all = grid_data['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(grid_data['M8_summed_FRP'])
    
    fit_points = (np.array([0,fire_max]))
    m8_fit = fit_points*params[1] + params[0]
    perf_fit = fit_points*1 + 0
    plot_ODR_scatter_residuals(grid_data['M11_summed_FRP'],grid_data['M8_summed_FRP'],grid_data['M11_FRP_uncertainty'],grid_data['M8_FRP_uncertainty'],m8_fit,perf_fit,fit_points, residual,adjusted_err,params,vza_max, fire_max,fire_count_all,'all',(BaseDir+'/Plots/FRP_ODR'+plot_tag))
        
    vza_bins =['all_deg','deg_10','deg10_20','deg20_30','deg30_40','deg40_']
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT': fire_count_all,'SLOPE':params[1] ,'INTERCEPT':params[0],'SE_slope':param_uncertainty[1],'SE_intercept':param_uncertainty[0]})
    for bin in vza_bins[1:]:
        actual_bins = np.unique(grid_data['vza_bin'])
        if bin in actual_bins:    
            m8_frp = grid_data.groupby(['vza_bin']).get_group(bin)['M8_summed_FRP']
            m11_frp = grid_data.groupby(['vza_bin']).get_group(bin)['M11_summed_FRP']
            m8_err = grid_data.groupby(['vza_bin']).get_group(bin)['M8_FRP_uncertainty']
            m11_err = grid_data.groupby(['vza_bin']).get_group(bin)['M11_FRP_uncertainty']
            vza_diff = grid_data.groupby(['vza_bin']).get_group(bin)['VZA_diff']  
            vza_max = vza_diff.max()
            
            params,param_uncertainty, adjusted_err, residual = odr_regression(m11_frp,m8_frp,m11_err,m8_err)
                    
            fire_count = len(m11_err)
            bias = grid_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].mean()
            SD = grid_data.groupby(['vza_bin']).get_group(bin)['FRP_diff'].std()
            rmsd = sqrt(bias**2 + SD**2)
            stat_title = ['BIAS','SCATTER','RMSD','FIRE_COUNT','SLOPE','INTERCEPT','SE_slope','SE_intercept'] ; stat_value = [bias,SD,rmsd, fire_count ,params[1],params[0],param_uncertainty[1],param_uncertainty[0]]
            for i in range(0,8,1):
                n = stats_f.VZA_BIN[stats_f.VZA_BIN == bin].index            
                stats_f.at[n,stat_title[i]] = stat_value[i]
            
            m8_fit = fit_points*params[1] + params[0]
            perf_fit = fit_points*1 + 0
            #textstr1 ='y = %.2fx + %.2f \n Grid-cell count=%d'%(slope, intercept,fire_count)
            textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias, SD, rmsd)
            plot_ODR_scatter_residuals(m11_frp,m8_frp,m11_err,m8_err,m8_fit,perf_fit,fit_points, residual,adjusted_err,params,vza_max,fire_max,fire_count,bin,(BaseDir+'/Plots/FRP_ODR'+plot_tag))
         
    return stats_f



dates = pd.date_range(start='2018-06-10',end='2018-07-15')
dates = dates.format(formatter=lambda x: x.strftime('%Y%m%d'))

dates = ['20180610','20180611','20180612' ,'20180613','20180616','20180617','20180618','20180619','20180621','20180624','20180625','20180626','20180630','20180701','20180702','20180704','20180705','20180707']

side = '_positive_'
plot_tag =''

grid_data = read_in_from_csv(dates,False,side)

#grid_data = grid_data[(grid_data['LONGITUDE'] > 5) & (grid_data['LONGITUDE'] < 25)]

#grid_data = grid_data[(grid_data['LATITUDE'] < 5) & (grid_data['LATITUDE'] > -5)]
#print grid_data[grid_data['TIME'] == 1315 ]

grid_data = grid_data[(grid_data['M11_summed_FRP'] > 50) & (grid_data['M8_summed_FRP'] > 50)]

if side == '_negative_':
    grid_data = grid_data[grid_data['VZA_diff'] < 0 ]
    plot_tag = side + plot_tag
if side == '_positive_':
    grid_data = grid_data[grid_data['VZA_diff'] > 0 ]
    plot_tag = side + plot_tag


grid_data['VZA_diff'] = grid_data['VZA_diff'].abs()

grid_data['vza_bin'] = np.where((grid_data.VZA_diff <= 10),
              'deg_10',
              np.where(np.logical_and(grid_data.VZA_diff > 10,grid_data.VZA_diff <= 20),            
              'deg10_20',
              np.where(np.logical_and(grid_data.VZA_diff > 20,grid_data.VZA_diff <= 30),
              'deg20_30',
              np.where(np.logical_and(grid_data.VZA_diff > 30,grid_data.VZA_diff <= 40),
              'deg30_40',
              'deg40_'))))

grid_data['LOG_VZA_RATIO'] = np.log(grid_data['VZA_RATIO'])
grid_data['LOG_FRP_RATIO'] = np.log(grid_data['FRP_RATIO'])


fit_and_plot_OLS_log(grid_data,BaseDir,'LOG_VZA_RATIO', 'LOG_FRP_RATIO',1 )

#plot_single_map(grid_data['LATITUDE'].values[:], grid_data['LONGITUDE'].values[:],grid_data['VZA_diff'].values[:], ('AREA'+plot_tag))
#fit_and_plot_OLS_bins(grid_data,BaseDir,'VZA_RATIO', 'FRP_RATIO',5.0 )
#fit_and_plot_ODR_bins(grid_data,BaseDir)
       
 
