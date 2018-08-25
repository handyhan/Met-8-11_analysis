import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import math
from general_plots import plot_scatter_simple
from regression_functions import regression_c0,regression_free,odr_regression


def get_grid_bounds(latmin,latmax,lonmin,lonmax,resolution):
    lon = np.arange(lonmin,(lonmax+resolution),resolution)
    lat = np.arange(latmin,(latmax+resolution),resolution)
    return lat,lon

def get_gridded(time,df,lat_bounds,lon_bounds, sat ,resolution):  #dataframe of pixels, location bounds, sat name as str and grid resolution
    
    start_time = time - pd.Timedelta(minutes=8)
    end_time = time + pd.Timedelta(minutes=8)    
    df = df.loc[start_time:end_time]
    df = df.reset_index()
    
    if df.empty :
        print "No pixels at time " + str(time)
        df_empty = pd.DataFrame()
        return df_empty                                                                       # pass if no concurrent pixels
    else:
       
        lats_binned = np.digitize(df['LATITUDE'].values[:],lat_bounds)       # generate an array with indicise indicating in which grid cell the lats/lons belong
        lon_binned = np.digitize(df['LONGITUDE'].values[:],lon_bounds)
        
        binned_df = pd.DataFrame({'LAT':lats_binned,'LONG':lon_binned})
        df['GRID_NO'] = binned_df.apply(lambda row: (str(row.LAT).zfill(3) + " " + str(row.LONG).zfill(3)) , axis=1)  # used the binned indecies to generate a unique marker for each cell 
        #make a gridcell  data frame from the cell markers
        df = df.sort_values(['GRID_NO'])
        fire_ids, pixel_count = np.unique(df['GRID_NO'], return_counts=True)                                #count the pixels in each grid cell
        df_ids = pd.DataFrame({'GRID_NO':fire_ids,'pixel_count':pixel_count})                               
        df['GRID_NO'] = df['GRID_NO'].apply(str)
        grid_df = pd.DataFrame(df_ids.GRID_NO.str.split(' ',1).tolist(), columns = ['lat_box','long_box'])
        grid_df.lat_box  = grid_df.lat_box.apply(lambda x: x.lstrip('0')).astype(int)
        grid_df.long_box  = grid_df.long_box.apply(lambda x: x.lstrip('0')).astype(int)
        # sum the FRP or the pixels in each gridcell and calculate the uncertainty
        frp_sum = df.groupby(['GRID_NO'])['FRP'].agg('sum')
        if sat == 'MODIS':
            pass
        else:
            df['FRP_UNCERTAINTY_2'] = df['FRP_UNCERTAINTY'].apply(lambda x: x**2 )                  #calculate uncertanty as root of sum of squared individual uncertainty
            frp_err = np.sqrt(df.groupby(['GRID_NO'])['FRP_UNCERTAINTY_2'].agg('sum'))
        #make cell bounds and centra
        lat_max = grid_df.lat_box.apply(lambda x:lat_bounds[(x)])                                           
        lat_min= grid_df.lat_box.apply(lambda x:lat_bounds[(x-1)])   
        lat_center = lat_min + resolution/2
        long_max = grid_df.long_box.apply(lambda x:lon_bounds[(x)])
        long_min = grid_df.long_box.apply(lambda x:lon_bounds[(x-1)])
        long_center = (long_max + long_min)/2
        # make df grid with all info
        grid = pd.DataFrame({'GRID_NO':fire_ids, 'summed_FRP': frp_sum.values[:], 'pixel_count':pixel_count , 
                             'SAT': sat,'LAT_MIN':lat_min,'LAT_MAX':lat_max,'LAT_CENTER':lat_center, 'LONG_MIN':long_min,'LONG_MAX':long_max,'LONG_CENTER':long_center, 'TIME_WINDOW':time })
        #get grid_cell VZA from precomputed gridded VZA map
        if sat == 'MODIS':
            pass
        else:
            grid['FRP_uncertainty'] = frp_err.values[:]
            vza_1deg = pd.read_csv(BaseDir + "/txt_data/MET_" + sat + "_1Deg_VZA")
            grid['loc_marker'] = grid['LAT_CENTER']*1000 + grid['LONG_CENTER']
            vza_1deg['loc_marker'] = vza_1deg['LATITUDE']*1000 + vza_1deg['LONGITUDE']  
            grid = grid.sort_values(['loc_marker'])
            unique = np.unique(grid['loc_marker'])
            vza_1deg = vza_1deg[vza_1deg['loc_marker'].isin(unique)]                # get only the VZA of grid cells with pixels in them
            vza_1deg = vza_1deg.sort_values(['loc_marker'])
            grid['VZA'] = vza_1deg['VZA'].values[:]
            grid = grid.drop(columns=['loc_marker'])
            
        return grid
    
    


def match_grids(date_range, msg_df, modis_df):    # takes input of MODIS_gridded and MSG_gridded files generated my proccesing script in df format with time as axis    
    first=True
    for k in range(0,len(date_range),1):        
        time = date_range[k]    
        start_time = time - pd.Timedelta(minutes=8)
        end_time = time + pd.Timedelta(minutes=8)    
        sub_mod = modis_df.loc[start_time:end_time]
        sub_msg = msg_df.loc[start_time:end_time]
        sub_msg = msg_df.loc[start_time:end_time]        
        if sub_mod.empty or sub_msg.empty:
            print "No Concurent gridcells at time " + str(time)                      ## pass if no concurrent pixels
            sub_msg_processed = pd.DataFrame()
        else:            
            unique_mod_cells = np.unique(sub_mod['GRID_NO'])                                # find grid cells at this time unique to both MODIS and MSG
            matching_msg_cells = sub_msg[sub_msg['GRID_NO'].isin(unique_mod_cells)]
            unique_msg_match_cells = np.unique(matching_msg_cells['GRID_NO'])
            mod_full_match = sub_mod[sub_mod['GRID_NO'].isin(unique_msg_match_cells)]
            unique_mod_cells = np.unique(mod_full_match['GRID_NO'])
            msg_full_match = matching_msg_cells[matching_msg_cells['GRID_NO'].isin(unique_mod_cells)]            
            msg_full_match.sort_values('GRID_NO')
            mod_full_match.sort_values('GRID_NO')
            if first == True:
                msg_match = msg_full_match
                mod_match = mod_full_match
            else:
                msg_match = msg_match.append(msg_full_match)                                # append all matches from all times together
                mod_match = mod_match.append(mod_full_match)                           
            first = False
    # combine matched grid cells into a single df
    matching_cells = mod_match
    matching_cells = matching_cells.rename(columns={'summed_FRP':'MODIS_FRP','pixel_count':'MODIS_PIXEL_COUTN'})
    matching_cells['MSG_FRP'] = msg_match['summed_FRP'].values[:]    
    matching_cells['MSG_PIXEL_COUNT'] = msg_match['pixel_count'].values[:]
    matching_cells = matching_cells.drop(columns=['SAT'])        
    matching_cells['FRP_DIFF'] = matching_cells['MODIS_FRP'] - matching_cells['MSG_FRP']       # calculate the difference in FRP values
    return matching_cells


def fit_and_plot_OLS_bins(matching_cells,x_string,y_string,scaling_factor):     #in put is df with matching cells, x and y are df labels for x and y variables, scaling factor, axis to the nearest scaling factor 
    x = x_string
    y = y_string
    results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(matching_cells, x,y) # run the regression ( regression_c0 can replace and has intecept at 0)
    ####remove outliers and re-fit
    #fire_data = fire_data[~fire_data.iloc[:].index.isin(outlier_ix)]
    #results, outlier_ix ,slope_all, intercept_all, r2_all, std_err_all, residual = regression_free(fire_data, 'M11_summed_FRP','M8_summed_FRP')
    #outlier_count_all = len(outlier_ix)
    #outlier_fires = fire_data[fire_data.iloc[:].index.isin(outlier_ix)]
    #plot_bin_maps(outlier_fires,fires_m8,'clusters')
    #######    
    fire_max = max(matching_cells[x].max(),matching_cells[y].max())
    fire_max  = int(math.ceil(fire_max / scaling_factor)) * scaling_factor
    bias_all = matching_cells['FRP_DIFF'].mean()                                        # calculate the mean bias
    SD_all = matching_cells['FRP_DIFF'].std()                                           # standard deviation
    rmsd_all = math.sqrt(bias_all**2 + SD_all**2)                                       # and rmsd
    fire_count_all = len(matching_cells[x])
    fit_points = (np.array([0,fire_max]))                                               # make points to plot the fit and the 1:1 line
    perf_fit = (fit_points)*1 + 0                                                       # 1:1
    z = (fit_points)*slope_all + intercept_all                                          # ols fit
    textstr1 = 'y = %.2fx + %.2f \n Grid-cell count=%d \n R$^{2}$=%.3f\n SE=%.3f'%(slope_all, intercept_all,fire_count_all, r2_all, std_err_all)
    textstr2 = 'Bias=%.3f \n  Scatter =%.3f \n RMSD=%.3f'%(bias_all, SD_all, rmsd_all)
    plot_scatter_simple(matching_cells[x],matching_cells[y],z,perf_fit,fit_points,(BaseDir+'/Plots/AREA_OLS'), textstr1,textstr2, fire_max, x,y,scaling_factor)  # plot the graph and save in set loc





"""
BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox" 
lon_max = 20
lon_min = -15
lat_max = -5
lat_min = -10
resolution = 1.0

    
modis_df = pd.read_csv(BaseDir + "/TEST_DATA/MODIS_gridded_full", parse_dates=True,index_col='TIME_WINDOW') 
msg_df = pd.read_csv(BaseDir + "/TEST_DATA/MSG_gridded_full", parse_dates=True,index_col='TIME_WINDOW') 

t = np.arange(datetime(2014,8,1), datetime(2014,9,1), timedelta(minutes=15)).astype(datetime)
matching_cells = match_grids(t, msg_df, modis_df)
matching_cells.to_csv(BaseDir + "/TEST_DATA/matching_cells", header=True,index=True) 

matching_cells = matching_cells[(matching_cells['MSG_FRP'] > 200) & (matching_cells['MODIS_FRP'] > 200)]
fit_and_plot_OLS_bins(matching_cells,'MODIS_FRP','MSG_FRP',5000)

"""





    



























