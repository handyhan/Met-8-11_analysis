# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 14:26:46 2018

@author: Hannah.N
"""



def regression_free(df,x_name,y_name):
    mod = ols(formula = y_name +' ~ '+ x_name , data = df)
    res= mod.fit()
    results = (res.summary())    
    test = res.outlier_test()
    outlier_ix = test[test['bonf(p)'] < 0.05].index
    slope = res.params[x_name]
    intercept = res.params.Intercept
    r_squared = res.rsquared
    se = res.bse[x_name]
    residual = res.resid
    
    return results, outlier_ix,slope, intercept,r_squared, se, residual

def regression_c0(df,x_name,y_name):
    mod = ols(formula = y_name +' ~ '+ x_name + '+0' , data = df)
    res= mod.fit()
    results = (res.summary())    
    test = res.outlier_test()
    outlier_ix = test[test['bonf(p)'] < 0.05].index
    slope = res.params[x_name]
    r_squared = res.rsquared
    se = se = res.bse[x_name]
    intercept = 0.00
    residual = res.resid
    
    return results, outlier_ix, slope , intercept, r_squared, se, residual



def odr_regression(x,y,x_err,y_err):
        
    def linear(p,x) :        # A linear function with:      #   Constant Background          : p[0]       #   Slope                        : p[1]
        return p[0]+p[1]*x

    func=linear                           # set the function to an object
    p_guess = (0,1)    # can choose first guess intercept to be zero by changing linear (replate p[0] with 0, then only fitting slope)      
    data = scodr.RealData(x=x, y=y, sx=x_err, sy=y_err)     
    model = scipy.odr.Model(func)                    
    odr = scodr.ODR(data, model, p_guess, maxit=5000,job=10)      
    output = odr.run()
    params = output.beta 	# 'beta' is an array of the parameter estimates
    param_uncertainty = output.sd_beta # parameter standard uncertainties    
    delta = output.delta  # estimated x-component of the residuals
    eps   = output.eps    # estimated y-component of the residuals
    xstar = x_err*np.sqrt( ((y_err*delta)**2) / ( (y_err*delta)**2 + (x_err*eps)**2 ) )
    ystar = y_err*np.sqrt( ((x_err*eps)**2) / ( (y_err*delta)**2 + (x_err*eps)**2 ) )
    adjusted_err = np.sqrt(xstar**2 + ystar**2)
    residual = np.sign(y-func(params,x))*np.sqrt(delta**2 + eps**2)
    return params,param_uncertainty, adjusted_err, residual



def fit_and_plot_OLS_bins(fire_data):  

    results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_summed_FRP','M8_summed_FRP') 
    ####remove outliers and re-fit
    #fire_data = fire_data[~fire_data.iloc[:].index.isin(outlier_ix)]
    #results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_summed_FRP','M8_summed_FRP')
    #outlier_count_all = len(outlier_ix)
    #outlier_fires = fire_data[fire_data.iloc[:].index.isin(outlier_ix)]
    #plot_bin_maps(outlier_fires,fires_m8,'clusters')
    #######
    
    vza_max = fire_data['VZA_diff'].max()
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
    plot_scatter_fit(fire_data['M11_summed_FRP'],fire_data['M8_summed_FRP'],fire_data['VZA_diff'],z,perf_fit,fit_points,'all',('./Plots/FRP_ols'+plot_tag), textstr1,textstr2,vza_max,fire_max)
    #plot_scatter_residuals(fire_data['M11_summed_FRP'],fire_data['M8_summed_FRP'],z, residual,slope_all,intercept_all,vza_max,fire_count_all,'all',('./Plots/FRP_OLS'+plot_tag))

    
    vza_bins =['all_deg','deg_10','deg10_20','deg20_30','deg30_40','deg40_']
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT':fire_count_all,'SLOPE':slope_all,'INTERCEPT':intercept_all, 'R2':r2_all,'SE':std_err_all})
    
    for bin in vza_bins[1:]:
        actual_bins = np.unique(fire_data['vza_bin'])
        if bin in actual_bins:    
            m8_frp = fire_data.groupby(['vza_bin']).get_group(bin)['M8_summed_FRP']
            m11_frp = fire_data.groupby(['vza_bin']).get_group(bin)['M11_summed_FRP']
            vza_diff = fire_data.groupby(['vza_bin']).get_group(bin)['VZA_diff']  
            group_fire = pd.DataFrame({'M11_summed_FRP':m11_frp,'M8_summed_FRP':m8_frp,'VZA_diff':vza_diff})
            if len(group_fire['M8_summed_FRP']) > 1 :
                results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, 'M11_summed_FRP','M8_summed_FRP')
                ###remove outliers and re-fit
                #group_fire = group_fire[~group_fire.iloc[:].index.isin(outlier_ix)]
                #results, outlier_ix ,slope, intercept, r2, std_err, residual = regression_free(group_fire, 'M11_summed_FRP','M8_summed_FRP')
                #outlier_count = len(outlier_ix)
                ######
                        
                fire_count = len(group_fire['M8_summed_FRP'])
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
                plot_scatter_fit(group_fire['M11_summed_FRP'],group_fire['M8_summed_FRP'],group_fire['VZA_diff'],z,perf_fit,fit_points,bin,('./Plots/FRP_ols'+plot_tag), textstr1,textstr2,vza_max,fire_max)
                #plot_scatter_residuals(group_fire['M11_summed_FRP'],group_fire['M8_summed_FRP'],z,residual,slope,intercept,vza_max,fire_count,bin,('./Plots/FRP_OLS'+plot_tag))
       

    return stats_f
       
    
def fit_and_plot_ODR_bins(fire_data):
    
    
    params,param_uncertainty, adjusted_err, residual = odr_regression(fire_data['M11_summed_FRP'],fire_data['M8_summed_FRP'],fire_data['M11_FRP_uncertainty'],fire_data['M8_FRP_uncertainty'])

    vza_max = fire_data['VZA_diff'].max()    
    fire_max = max(fire_data['M11_summed_FRP'].max(),fire_data['M8_summed_FRP'].max())
    bias_all = fire_data['FRP_diff'].mean()
    SD_all = fire_data['FRP_diff'].std()
    rmsd_all = sqrt(bias_all**2 + SD_all**2)
    fire_count_all = len(fire_data['M8_summed_FRP'])
    
    fit_points = (np.array([0,max(fire_data['M11_summed_FRP'].values)]))
    m8_fit = fit_points*params[1] + params[0]
    perf_fit = fit_points*1 + 0
    plot_ODR_scatter_residuals(fire_data['M11_summed_FRP'],fire_data['M8_summed_FRP'],fire_data['M11_FRP_uncertainty'],fire_data['M8_FRP_uncertainty'],m8_fit,perf_fit,fit_points, residual,adjusted_err,params,vza_max, fire_max,fire_count_all,'all',('./Plots/FRP_ODR'+plot_tag))
        
    vza_bins =['all_deg','deg_10','deg10_20','deg20_30','deg30_40','deg40_']
    stats_f = pd.DataFrame({'VZA_BIN':vza_bins,'BIAS':bias_all,'SCATTER':SD_all,'RMSD':rmsd_all,'FIRE_COUNT': fire_count_all,'SLOPE':params[1] ,'INTERCEPT':params[0],'SE_slope':param_uncertainty[1],'SE_intercept':param_uncertainty[0]})
    print "hello"
    for bin in vza_bins[1:]:
        actual_bins = np.unique(fire_data['vza_bin'])
        if bin in actual_bins:    
            print "hello2"
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
            textstr1 ='y = %.2fx + %.2f \n Fire count=%d'%(slope, intercept,fire_count)
            textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias, SD, rmsd)
            plot_ODR_scatter_residuals(m11_frp,m8_frp,m11_err,m8_err,m8_fit,perf_fit,fit_points, residual,adjusted_err,params,vza_max,fire_max,bin,('./Plots/FRP_ODR'+plot_tag))
         
    return stats_f



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
  
    

dates = pd.date_range(start='2018-06-09',end='2018-06-20')
dates = dates.format(formatter=lambda x: x.strftime('%Y%m%d'))
plot_tag = "_post_QP_negative_" 

fire_data, fires_m8, fires_m11, pixels_m8, pixels_m11 = read_in_from_csv(dates)
fire_data = fire_data[fire_data['VZA_diff'] < 0 ]
fire_data = fire_data[fire_data['M11_summed_FRP'] > 50 ]
fire_data = fire_data[fire_data['M8_summed_FRP'] > 50 ]



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


list_fires = np.unique(fire_data['cluster_no'])
select_fires = pixels_m8[pixels_m8['cluster_no'].isin(list_fires)]
plot_bin_maps(fire_data,select_fires,'fire_locations')
  
#stats_f = fit_and_plot_ODR_bins(fire_data)
stats_f = fit_and_plot_OLS_bins(fire_data)
#labels = ['Bias','Scatter','RMSD']
#plot_stats(stats_f['VZA_BIN'][:],stats_f['BIAS'][:],stats_f['SCATTER'][:],stats_f['RMSD'][:],'./Plots/Bias_scatter_rmsd',plot_tag,labels)
#labels = ['Slope','R$^2$','Standard Error']
#plot_stats(stats_f['VZA_BIN'][:],stats_f['SLOPE'][:],stats_f['R2'][:],stats_f['SE'][:],'./Plots/LinReg_params',plot_tag,labels)
#stats_f.to_csv('./fire_stats_table' +plot_tag +'.csv')

