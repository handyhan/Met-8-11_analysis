import h5py
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import sys
import math
from MSG_read import plot_single_map
from MSG_read import extract_data #,plot_compare,plot_scatter_fit,plot_single_map
from math import sin, cos, sqrt, atan2, radians
from scipy import stats
from datetime import datetime
pd.set_option('display.max_columns', 500)
pd.options.mode.chained_assignment = None  

cluster_threshold = 10
match_threshold = 10


def distance_on_unit_sphere(lat1, long1, lat2, long2):
    earth_R = float(6371) 
    lat1 = radians(lat1) ; lat2 = radians(lat2) ;lon1 = radians(long1); lon2 = radians(long2) 
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = earth_R * c
    return distance

def cluster_find(df,sat):
    df['LATITUDE'] = df['LATITUDE'].astype(float).round(4)
    df['LONGITUDE'] = df['LONGITUDE'].astype(float).round(4)
    df['cluster_no'] = ((df['LATITUDE']*10).round(0)) +((df['LONGITUDE']*100000).round(decimals=-5))  
     
    #df = df[['LATITUDE','LONGITUDE','cluster_no']]
    #print df#['cluster_no']
    #for i in range(0,len(df['cluster_no'])):
    #    print df['cluster_no'][i]
    #df['cluster_no'] = df['cluster_no'].round(decimals=-4)  
    unique, counts = np.unique(df['cluster_no'], return_counts=True)
    
    cluster_top = pd.DataFrame({'cluster_no':unique,'pixel_count': counts})
    #print cluster_top['pixel_count'].sort_values()
    cluster_top = cluster_top[cluster_top['pixel_count'] >= 3]
    df = df[df['cluster_no'].isin(cluster_top['cluster_no'].values[:])]
    df['SAT'] = sat
    #df_pos_pix = df_pixels[df_pixels['cluster_no'] > 0]2800000.0
    #df['cluster_no'] = day+time + df_pos_pix['cluster_no'].astype(str) ;     df_pos_pix['cluster_no'] =df_pos_pix['cluster_no'].astype(float)
    
    #print df
    #plot_single_map(df['LATITUDE'].values[:],df['LONGITUDE'].values[:],("TEST_pix_lats"+sat))
    return df



def proccess_cluster_data(df_pixels, sat , date,time):
    day = date[4:8]
    datetime_str = str(date)+str(time)    
    date_time_obj = datetime.strptime(datetime_str, '%Y%m%d%H%M')
    #print df_pixels[df_pixels['cluster_no'] < 0]
    df_neg_pix = df_pixels[df_pixels['cluster_no'] < 0]
    df_neg_pix['cluster_no'] = df_neg_pix['cluster_no'].abs()
    df_neg_pix['cluster_no'] = '-'+day+time + df_neg_pix['cluster_no'].astype(str);     df_neg_pix['cluster_no'] =df_neg_pix['cluster_no'].astype(float)
    df_pos_pix = df_pixels[df_pixels['cluster_no'] > 0]
    df_pos_pix['cluster_no'] = day+time + df_pos_pix['cluster_no'].astype(str) ;     df_pos_pix['cluster_no'] =df_pos_pix['cluster_no'].astype(float)
    df_pixels = pd.concat([df_neg_pix,df_pos_pix])
    df_pixels = df_pixels.sort_values(['cluster_no'],ascending = True)
    unique, counts = np.unique(df_pixels['cluster_no'], return_counts=True) 
    df_1_frp_sum = df_pixels.groupby(['cluster_no'])['FRP'].agg('sum')
    df_1_vza_mean = df_pixels.groupby(['cluster_no'])['PIXEL_VZA'].mean()
    df_1_lat_mean = df_pixels.groupby(['cluster_no'])['LATITUDE'].mean()
    df_1_lon_mean = df_pixels.groupby(['cluster_no'])['LONGITUDE'].mean()
    df_pixels['FRP_UNCERTAINTY_2'] = df_pixels['FRP_UNCERTAINTY'].apply(lambda x: x**2 )
    df_1_frp_err = np.sqrt(df_pixels.groupby(['cluster_no'])['FRP_UNCERTAINTY_2'].agg('sum'))
    df_pixels.drop('FRP_UNCERTAINTY_2', axis=1,inplace=True)  
    df_1_frp_sum = pd.DataFrame({'cluster_no':unique, 'summed_FRP': df_1_frp_sum.values[:],'FRP_uncertainty':df_1_frp_err.values[:],'mean_vza': df_1_vza_mean, 'pixel_count':counts , 'SAT': sat,'mean_lat':df_1_lat_mean,'mean_long':df_1_lon_mean })
    df_1_frp_sum['DATE_TIME'] = date_time_obj
    df_1_frp_sum = df_1_frp_sum.dropna()
    df_clusters = df_1_frp_sum.reset_index(drop=True)
    #df_clusters = df_clusters[df_clusters['cluster_no'] == 53011301999906.0 ]#5301130234.1
    #df_pixels = df_pixels[df_pixels['cluster_no'] == 53011301999906.0 ]
    #print df_pixels['cluster_no'][838],df_pixels['cluster_no'][411],df_pixels['cluster_no'][835],df_pixels['cluster_no'][422]
    #plot_single_map(df_pixels['LATITUDE'].values[:],df_pixels['LONGITUDE'].values[:],("TEST_pix_lats_pix"+sat))
    #plot_single_map(df_clusters['mean_lat'].values[:],df_clusters['mean_long'].values[:],("TEST_pix_lats"+sat))

    return df_clusters, df_pixels



def eliminate_non_match(df1_pixels,df2_pixels,df1_fires,df2_fires):
    df1_unique = np.unique(df1_fires['cluster_no']) 
    df2_unique = np.unique(df2_fires['cluster_no']) 
    clusters = np.concatenate((df1_unique,df2_unique))
    cluster = pd.DataFrame(data = clusters)
    cluster = np.unique(cluster[cluster.duplicated()])
    #unique_fires = np.unique(unique_fires)
    #print len(cluster)
    
    df1_pixels = df1_pixels[df1_pixels['cluster_no'].isin(cluster)]
    df2_pixels = df2_pixels[df2_pixels['cluster_no'].isin(cluster)]    
    df1_fires = df1_fires[df1_fires['cluster_no'].isin(cluster)]
    df2_fires = df2_fires[df2_fires['cluster_no'].isin(cluster)]
    df_fires = pd.concat([df1_fires,df2_fires])
    df_fires = df_fires.reset_index(drop=True)
    df_fires = df_fires.sort_values(by =['cluster_no','SAT'],ascending=True)
    df_pixels = pd.concat([df1_pixels,df2_pixels])
    df_pixels = df_pixels.reset_index(drop=True)
    df_pixels = df_pixels.sort_values(by =['cluster_no','SAT'],ascending=True)
    
    #plot_single_map(df_pixels['LATITUDE'].values[:],df_pixels['LONGITUDE'].values[:],("TEST__pix"))
    #plot_single_map(df_fires['mean_lat'].values[:],df_fires['mean_long'].values[:],("TEST__fires"))
    
    return  df_fires,df_pixels
    
    #print df1_pixels
    

dates = pd.date_range(start='2018-06-09',end='2018-06-20')
dates = dates.format(formatter=lambda x: x.strftime('%Y%m%d'))
#dates = ["20180609"]
for date in dates:
    print "NOTE: Generating fires for "  + date
    times = pd.date_range(start='00:00',end='23:45', freq='15min')
    times=times.format(formatter=lambda x: x.strftime('%H%M'))
    header_1=True
    header_2=True
    
    for time in times:

        print "Finding and matching clusters " + time
        m8_df,m11_df = extract_data(date,time)
        
        if m8_df.empty or m11_df.empty:
            print "ERROR while opening file"
            
        else:
            df_m8_pixels = cluster_find(m8_df,'8')
            df_m11_pixels = cluster_find(m11_df,'11')
            
            if df_m8_pixels.empty or df_m11_pixels.empty:
                print "NO CLUSTERS FOUND AT TIME " + time
            else:
                
                df_m8_fires,df_m8_pixels = proccess_cluster_data(df_m8_pixels, '8' ,date,time)
                df_m11_fires,df_m11_pixels = proccess_cluster_data(df_m11_pixels, '11' ,date,time)
                
                df_fires,df_pixels = eliminate_non_match(df_m8_pixels,df_m11_pixels,df_m8_fires,df_m11_fires)
                mode = 'w' if header_1 else 'a'
                df_fires.to_csv("./txt_data/"+date+"_fires", mode=mode,header=header_1,index=False) 
                header_1=False
                mode = 'w' if header_2 else 'a'
                df_pixels.to_csv("./txt_data/"+date+"_pixels", mode=mode,header=header_2,index=False) 
                header_2=False
              

    
   
"""
    
    #write_daily_fires(date)
    
    
    
    

def cluster_find(df):
    cluster_list = np.zeros(shape=(2,3))
    n = 1
    for i in range(0,(len(df['LATITUDE'])),1):
        for j in range(0,(len(df['LATITUDE'])),1):
            point_distance = distance_on_unit_sphere(df['LATITUDE'][i],df['LONGITUDE'][i], df['LATITUDE'][j],df['LONGITUDE'][j])
            if 0 <= point_distance < cluster_threshold:
                if cluster_list[-1,0] == i:                    
                    new_row = np.array([i,j,n])
                else:
                    n=n+1 ; new_row = np.array([i,j,n])
                 
                cluster_list = np.vstack((cluster_list, new_row))
                   
    cluster_list = np.delete(cluster_list,(0,1),axis=0) 
    cluster_list = pd.DataFrame(cluster_list, columns = ['pixel_i','pixel_j','cluster_no'])
    unique, counts = np.unique(cluster_list['pixel_i'], return_counts=True)
    count_df = pd.DataFrame({'pixel_#':unique,'pixel_count': counts})
    cluster_top = count_df#[0:keep_top_n]
    cluster_top = cluster_top[cluster_top['pixel_count'] >= 5]
    df_top =  np.array(cluster_top['pixel_#'].values[:])    
    cluster_list = cluster_list[cluster_list['pixel_i'].isin(df_top)]
    cluster_list = cluster_list.dropna()
    cluster_list=cluster_list.drop_duplicates('pixel_j')
    cluster_list = cluster_list.sort_values(by=['pixel_j'],ascending=True)    
    cluster_df = cluster_list['pixel_j'].values[:]
    cluster_df = df.loc[cluster_df]
    cluster_df['cluster_no'] = cluster_list['cluster_no'].values[:]
    cluster_df = cluster_df.reset_index(drop=True)
    unique, counts = np.unique(cluster_df['cluster_no'], return_counts=True)
    cluster_top = pd.DataFrame({'cluster_no':unique,'pixel_count': counts})
    cluster_top = cluster_top[cluster_top['pixel_count'] >= 5]
    df_top =  np.array(cluster_top['cluster_no'].values[:])    
    cluster_df = cluster_df[cluster_df['cluster_no'].isin(df_top)]
    cluster_df = cluster_df.reset_index(drop=True)
    lat_mean = cluster_df.groupby(['cluster_no'])['LATITUDE'].mean()
    lon_mean = cluster_df.groupby(['cluster_no'])['LONGITUDE'].mean()
    cluster_no = np.unique(cluster_df['cluster_no'])
    df_check = pd.DataFrame({'cluster_no':cluster_no,'mean_LAT': lat_mean.values[:],'mean_LON':lon_mean.values[:]})
    
    # check clusters don't belong to same cluster
    for i in range(0,len(df_check),1):
        for j in range(0,len(df_check),1):
            point_distance = distance_on_unit_sphere(df_check['mean_LAT'][i],df_check['mean_LON'][i],df_check['mean_LAT'][j],df_check['mean_LON'][j])
            if (0 <= point_distance < cluster_threshold) and (df_check['cluster_no'][i] != df_check['cluster_no'][j]):                
                pixels =  cluster_df[cluster_df['cluster_no'] == df_check['cluster_no'][j]] 
                pixels['cluster_no'][:] = df_check['cluster_no'][i]
                cluster_df = cluster_df[cluster_df['cluster_no'] != df_check['cluster_no'][j]]
                cluster_df = pd.concat([cluster_df,pixels])
                df_check.at[j,'cluster_no'] = df_check['cluster_no'][i]   
    
    cluster_df = cluster_df.sort_values(by=['cluster_no'],ascending=True)
    cluster_df = cluster_df.reset_index(drop=True)                
    return cluster_df


def cluster_match(df_1_clusters,df_2_full, sat_1, sat_2 , ):

    cluster_list = np.zeros(shape=(2,3))
    n = 1
    for i in range(0,(len(df_1_clusters['LATITUDE'])),1):
         for j in range(0,(len(df_2_full['LATITUDE'])),1):
             point_distance = distance_on_unit_sphere(df_1_clusters['LATITUDE'][i],df_1_clusters['LONGITUDE'][i], df_2_full['LATITUDE'][j],df_2_full['LONGITUDE'][j])
             if 0 <= point_distance < match_threshold:
                 if cluster_list[-1,0] == i:                    
                     new_row = np.array([i,j,df_1_clusters['cluster_no'][i]])
                 else:
                     n=n+1 ; new_row = np.array([i,j,df_1_clusters['cluster_no'][i]])
                 
                 cluster_list = np.vstack((cluster_list, new_row))
         
    cluster_list = np.delete(cluster_list,(0,1),axis=0) 
    cluster_list = pd.DataFrame(cluster_list, columns = ['pixel_'+sat_1 ,'pixel_'+sat_2,'cluster_no'])
    cluster_list = cluster_list.drop_duplicates('pixel_'+sat_2)
    cluster_list = cluster_list.sort_values(by=['pixel_'+sat_2],ascending=True)
    cluster_df = cluster_list['pixel_'+sat_2].values[:]
    cluster_df = df_2_full.loc[cluster_df]
    cluster_df['cluster_no'] = cluster_list['cluster_no'].values[:]
    list_unique_match = np.unique(cluster_df['cluster_no'])               # if unmatch found clusters
    df_1_clusters = df_1_clusters[df_1_clusters['cluster_no'].isin(list_unique_match)]
    cluster_df = cluster_df.sort_values(by=['cluster_no'],ascending=True)
    cluster_df = cluster_df.reset_index(drop=True)
    #lat_mean = cluster_df.groupby(['cluster_no'])['LATITUDE'].mean()
    #lon_mean = cluster_df.groupby(['cluster_no'])['LONGITUDE'].mean()
    
    return df_1_clusters, cluster_df
 


def proccess_cluster_data(cluster_marker,df_1_clusters, df_2_matches,sat_1,sat_2, date,time):
    day = date[4:8]
    datetime_str = str(date)+str(time)    
    date_time_obj = datetime.strptime(datetime_str, '%Y%m%d%H%M')
    df_1_clusters['cluster_no'] = day+time+cluster_marker + df_1_clusters['cluster_no'].astype(str); df_1_clusters['cluster_no'] =df_1_clusters['cluster_no'].astype(float)
    df_1_unique, df_1_counts = np.unique(df_1_clusters['cluster_no'],return_counts=True)    
    #cluster_top = pd.DataFrame({'cluster_no':df_1_unique,'pixel_count': df_1_counts})
    df_1_clusters['SAT'] = sat_1 
    df_1_frp_sum = df_1_clusters.groupby(['cluster_no'])['FRP'].agg('sum')
    df_1_vza_mean = df_1_clusters.groupby(['cluster_no'])['PIXEL_VZA'].mean()
    df_1_lat_mean = df_1_clusters.groupby(['cluster_no'])['LATITUDE'].mean()
    df_1_lon_mean = df_1_clusters.groupby(['cluster_no'])['LONGITUDE'].mean()
    df_1_clusters['FRP_UNCERTAINTY_2'] = df_1_clusters['FRP_UNCERTAINTY'].apply(lambda x: x**2 )
    df_1_frp_err = np.sqrt(df_1_clusters.groupby(['cluster_no'])['FRP_UNCERTAINTY_2'].agg('sum'))
    df_1_clusters.drop('FRP_UNCERTAINTY_2', axis=1,inplace=True)  
    df_1_frp_sum = pd.DataFrame({'cluster_no':df_1_unique, 'summed_FRP': df_1_frp_sum.values[:],'FRP_uncertainty':df_1_frp_err.values[:],'mean_vza': df_1_vza_mean, 'pixel_count':df_1_counts , 'SAT': sat_1,'mean_lat':df_1_lat_mean,'mean_long':df_1_lon_mean })
    df_1_frp_sum['DATE_TIME'] = date_time_obj 
    
    df_1_frp_sum = df_1_frp_sum.reset_index(drop=True)
    
    df_2_matches['cluster_no'] = day+time +cluster_marker+ df_2_matches['cluster_no'].astype(str); df_2_matches['cluster_no'] = df_2_matches['cluster_no'].astype(float)    
    df_2_unique, df_2_counts = np.unique(df_2_matches['cluster_no'],return_counts=True)
    #cluster_top = pd.DataFrame({'cluster_no':df_2_unique,'pixel_count': df_2_counts})
    df_2_matches['SAT'] = sat_2            
    df_2_frp_sum = df_2_matches.groupby(['cluster_no'])['FRP'].agg('sum')
    df_2_vza_mean = df_2_matches.groupby(['cluster_no'])['PIXEL_VZA'].mean()    
    df_2_lat_mean = df_2_matches.groupby(['cluster_no'])['LATITUDE'].mean()
    df_2_lon_mean = df_2_matches.groupby(['cluster_no'])['LONGITUDE'].mean()
    df_2_matches['FRP_UNCERTAINTY_2'] = df_2_matches['FRP_UNCERTAINTY'].apply(lambda x: x**2 )
    df_2_frp_err = np.sqrt(df_2_matches.groupby(['cluster_no'])['FRP_UNCERTAINTY_2'].agg('sum'))
    df_2_matches.drop('FRP_UNCERTAINTY_2', axis=1,inplace=True)
    df_2_frp_sum = pd.DataFrame({'cluster_no':df_2_unique, 'summed_FRP': df_2_frp_sum.values[:], 'FRP_uncertainty':df_2_frp_err.values[:],'mean_vza': df_2_vza_mean, 'pixel_count':df_2_counts, 'SAT': sat_2, 'mean_lat':df_2_lat_mean,'mean_long':df_2_lon_mean })
    df_2_frp_sum['DATE_TIME'] = date_time_obj 
    df_2_frp_sum = df_2_frp_sum.reset_index(drop=True)
    
    return df_1_clusters, df_2_matches, df_1_frp_sum, df_2_frp_sum


def write_daily_fires(date): # valud arguments 'M8' or 'M11'
    
    times = pd.date_range(start='10:30',end='10:30', freq='15min')
    times=times.format(formatter=lambda x: x.strftime('%H%M'))
    header_1=True
    header_2=True
    
    for time in times:
        date_time= date + time
        m8_df,m11_df = extract_data(date,time)
        if m8_df.empty or m11_df.empty:
            print "ERROR while opening file"
            
        else:
            print "Finding and matching clusters " + time
            
            m11_clusters = cluster_find(m11_df) 
            m11_clusters, m11_matchs = cluster_match(m11_clusters,m8_df, '11','8')
            df_m11_pixels,df_m11_pixel_matches, df_m11_fires, df_m11_fire_matches = proccess_cluster_data('1', m11_clusters, m11_matchs,'11','8', date,time)

            m8_clusters = cluster_find(m8_df)
            m8_clusters, m8_matchs = cluster_match(m8_clusters,m11_df, '8','11')
            df_m8_pixels,df_m8_pixel_matches, df_m8_fires, df_m8_fire_matches = proccess_cluster_data('8', m8_clusters, m8_matchs,'8','11', date,time)

            df_m11_fires['same_marker'] = 0 ; df_m11_fire_matches['same_marker'] = 0 ; df_m8_fires['same_marker'] = 0 ; df_m8_fire_matches['same_marker'] = 0
            for i in range(0,len(df_m11_fires.mean_lat),1):
                for j in range(0,len(df_m8_fires.mean_lat),1):
                    lat1 =  df_m11_fires['mean_lat'].loc[i] ; long1 =  df_m11_fires['mean_long'].loc[i]
                    lat2 =  df_m8_fires['mean_lat'].loc[j] ; long2 =  df_m8_fires['mean_long'].loc[j]
                    point_distance = distance_on_unit_sphere(lat1,long1,lat2,long2)
                    if 0 <= point_distance < match_threshold:
                        df_m11_fires.at[i,'same_marker'] = np.float(str(i)+str(j))
                        df_m11_fire_matches.at[i,'same_marker'] = np.float(str(i)+str(j))
                        df_m8_fires.at[j,'same_marker'] = np.float(str(i)+str(j))
                        df_m8_fire_matches.at[j,'same_marker'] = np.float(str(i)+str(j))
                        
                        
            print "hello3"
                        
                        
            df_unique_m11_fires = df_m11_fires.loc[df_m11_fires['same_marker'] == 0 ]  #; df2_m11_fires.drop('same_marker', axis=1, inplace=True)
            df_unique_m11_fire_matches = df_m11_fire_matches.loc[df_m11_fire_matches['same_marker'] == 0 ] #; df2_m8_fires.drop('same_marker', axis=1, inplace=True)
            df_unique_m8_fires = df_m8_fires.loc[df_m8_fires['same_marker'] == 0 ] #; df2_m8_fires.drop('same_marker', axis=1, inplace=True)
            df_unique_m8_fire_matches = df_m8_fire_matches.loc[df_m8_fire_matches['same_marker'] == 0 ]  #; df2_m11_fires.drop('same_marker', axis=1, inplace=True)


            df_co_m11_fires = df_m11_fires.loc[df_m11_fires['same_marker'] != 0 ].dropna()  #; df2_m11_fires.drop('same_marker', axis=1, inplace=True)
            df_co_m11_fire_matches = df_m11_fire_matches.loc[df_m11_fire_matches['same_marker'] != 0 ].dropna() #; df2_m8_fires.drop('same_marker', axis=1, inplace=True)

            df_co_m8_fires = df_m8_fires.loc[df_m8_fires['same_marker'] != 0 ].dropna()#; df2_m8_fires.drop('same_marker', axis=1, inplace=True)

            df_co_m8_fire_matches = df_m8_fire_matches.loc[df_m8_fire_matches['same_marker'] != 0 ].dropna()  #; df2_m11_fires.drop('same_marker', axis=1, inplace=True)

            
            list_unique_m11_pixles = df_co_m11_fire_matches['same_marker'][:]
            missing = df_co_m8_fire_matches[~df_co_m8_fire_matches['same_marker'].isin(list_unique_m11_pixles)]
            print missing
            print df_co_m11_fire_matches
            print df_co_m8_fire_matches
            
            
            
            list_unique_m11_pixles = df_unique_m11_fires['cluster_no'][:]
            df_unique_m11_pixles = df_m11_pixels[df_m11_pixels['cluster_no'].isin(list_unique_m11_pixles)]
            df_unique_m11_pixle_match = df_m11_pixel_matches[df_m11_pixel_matches['cluster_no'].isin(list_unique_m11_pixles)]
            list_unique_m8_pixles = df_unique_m8_fires['cluster_no'][:]
            df_unique_m8_pixles = df_m8_pixels[df_m8_pixels['cluster_no'].isin(list_unique_m8_pixles)]
            df_unique_m8_pixle_match = df_m8_pixel_matches[df_m8_pixel_matches['cluster_no'].isin(list_unique_m8_pixles)]
            
            list_co_m11_pixle_match = df_co_m11_fire_matches['cluster_no'][:]
            df_co_m11_pixle_match = df_m11_pixel_matches[df_m11_pixel_matches['cluster_no'].isin(list_co_m11_pixle_match)]
            list_co_pixle_matches = np.array([df_co_m8_fire_matches['cluster_no'].values[:],df_co_m11_fire_matches['cluster_no'].values[:]])
            
            for i in range(0,len(list_co_pixle_matches[0][:]),1):
                match_pixels_m8 = df_m8_pixel_matches[df_m8_pixel_matches['cluster_no'] == (list_co_pixle_matches[0][i])]
                match_pixels_m8['cluster_no'] = list_co_pixle_matches[1][i]
                df_co_m11_pixle_match = pd.concat([df_co_m11_pixle_match,match_pixels_m8])
                
            
            df_all_top_pixels = pd.concat([df_co_m11_pixle_match,df_unique_m11_pixles,df_unique_m8_pixles,df_unique_m8_pixle_match])
            df_all_top_pixels = df_all_top_pixels.sort_values(by=['cluster_no', 'SAT'],ascending=True).reset_index(drop=True)
            df_all_top_pixels.reset_index(drop=True,inplace = True)
            
            df_co_m8_fire_matches['cluster_no'][:] = df_co_m11_fire_matches['cluster_no'][:]
            df_all_top_fires = pd.concat([df_co_m8_fire_matches, df_co_m11_fire_matches,df_unique_m11_fires,df_unique_m11_fire_matches,df_unique_m8_fires,df_unique_m8_fire_matches])
            df_all_top_fires = df_all_top_fires.sort_values(by=['cluster_no', 'SAT'],ascending=True).reset_index(drop=True)
            df_all_top_pixels.reset_index(drop=True,inplace = True)
            

            mode = 'w' if header_1 else 'a'
            df_all_top_fires.to_csv("./txt_data/"+date+"_fires", mode=mode,header=header_1,index=False) 
            header_1=False
            mode = 'w' if header_2 else 'a'
            df_all_top_pixels.to_csv("./txt_data/"+date+"_pixels", mode=mode,header=header_2,index=False) 
            header_2=False
          

    
    
#dates = ["20180609"]#,'20180602','20180603']#["20180515","20180516","20180517","20180518","20180519","20180520","20180521","20180522","20180523","20180524","20180525","20180526","20180527","20180528","20180529","20180530","20180531"]


#plot_single_map(df_co_m11_fires, df_co_m8_fires,df_co_m11_fire_matches,df_co_m8_fire_matches, "all_clusters")
"""