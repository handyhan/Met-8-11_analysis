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
#from sklearn import datasets, linear_model
#from sklearn.metrics import mean_squared_error, r2_score

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


def cluster_find(df,keep_top_n):
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
    count_df = count_df.sort_values(by=['pixel_count'],ascending=False)
    cluster_top = count_df#[0:keep_top_n]
    cluster_top = cluster_top[cluster_top['pixel_count'] >= 3]
    df_top =  np.array(cluster_top['pixel_#'].values[:])    
    cluster_list = cluster_list[cluster_list['pixel_i'].isin(df_top)]
    cluster_list = cluster_list.dropna()
    cluster_list=cluster_list.drop_duplicates('pixel_j')
    cluster_list = cluster_list.sort_values(by=['pixel_j'],ascending=True)    
    cluster_df = cluster_list['pixel_j'].values[:]
    cluster_df = df.loc[cluster_df]
    cluster_df['cluster_no'] = cluster_list['cluster_no'].values[:]
    cluster_df = cluster_df.reset_index(drop=True)
    return cluster_df

def cluster_match(df_1_clusters,df_2_full, sat_1, sat_2 , ):
    cluster_list = np.zeros(shape=(2,4))
    n = 1
    for i in range(0,(len(df_1_clusters['LATITUDE'])),1):
         for j in range(0,(len(df_2_full['LATITUDE'])),1):
             point_distance = distance_on_unit_sphere(df_1_clusters['LATITUDE'][i],df_1_clusters['LONGITUDE'][i], df_2_full['LATITUDE'][j],df_2_full['LONGITUDE'][j])
             if 0 <= point_distance < match_threshold:
                 if cluster_list[-1,0] == i:                    
                     new_row = np.array([i,j,n,df_1_clusters['cluster_no'][i]])
                 else:
                     n=n+1 ; new_row = np.array([i,j,n,df_1_clusters['cluster_no'][i]])
                 
                 cluster_list = np.vstack((cluster_list, new_row))
         
    cluster_list = np.delete(cluster_list,(0,1),axis=0) 
    cluster_list = pd.DataFrame(cluster_list, columns = ['pixel_'+sat_1 ,'pixel_'+sat_2,'cluster_no',sat_1+'_cluster_no'])
    cluster_list = cluster_list.drop_duplicates('pixel_'+sat_2)
    cluster_list = cluster_list.sort_values(by=['pixel_'+sat_2],ascending=True)
    cluster_df = cluster_list['pixel_'+sat_2].values[:]
    cluster_df = df_2_full.loc[cluster_df]
    cluster_df['cluster_no'] = cluster_list['cluster_no'].values[:]
    cluster_df[sat_1+'_cluster_no'] = cluster_list[sat_1+'_cluster_no'].values[:]
    unique, counts = np.unique(cluster_list['cluster_no'], return_counts=True)
    cluster_df = cluster_df.reset_index(drop=True)
    return cluster_df
        
def proccess_cluster_data(cluster_marker,df_1_clusters, df_2_matches,sat_1,sat_2, date,time):
    day = date[6:8]
    datetime_str = str(date)+str(time)
    date_time_obj = datetime.strptime(datetime_str, '%Y%m%d%H%M')
    seconds_since = (date_time_obj-datetime(1970,1,1)).total_seconds()
    largest_clusters = df_1_clusters['cluster_no'].value_counts().index.tolist()[0:10]
    top_df_2_matches = df_2_matches[df_2_matches[sat_1+'_cluster_no'].isin(largest_clusters)]    
    top_df_1_clusters =df_1_clusters[df_1_clusters['cluster_no'].isin(top_df_2_matches[sat_1+'_cluster_no'].unique())]
    df_1_unique, df_1_counts = np.unique(top_df_1_clusters['cluster_no'],return_counts=True)
    top_df_1_clusters['SAT'] = sat_1 
    top_df_1_clusters['cluster_no'] = day+time+cluster_marker + top_df_1_clusters['cluster_no'].astype(str); top_df_1_clusters['cluster_no'] =top_df_1_clusters['cluster_no'].astype(float)
   
    df_2_unique, df_2_counts = np.unique(top_df_2_matches[sat_1+'_cluster_no'],return_counts=True)
    top_df_2_matches['cluster_no'].values[:] = top_df_2_matches[sat_1+'_cluster_no'] .values[:]
    top_df_2_matches.drop(columns = [sat_1+'_cluster_no'],inplace=True)
    top_df_2_matches['SAT'] = sat_2            
    top_df_2_matches['cluster_no'] = day+time +cluster_marker+ top_df_2_matches['cluster_no'].astype(str); top_df_2_matches['cluster_no'] = top_df_2_matches['cluster_no'].astype(float)
    df_1_frp_sum = top_df_1_clusters.groupby(['cluster_no'])['FRP'].agg('sum')
    df_1_vza_mean = top_df_1_clusters.groupby(['cluster_no'])['PIXEL_VZA'].mean()
    df_1_lat_mean = top_df_1_clusters.groupby(['cluster_no'])['LATITUDE'].mean()
    df_1_lon_mean = top_df_1_clusters.groupby(['cluster_no'])['LONGITUDE'].mean()
    top_df_1_clusters['FRP_UNCERTAINTY_2'] = top_df_1_clusters['FRP_UNCERTAINTY'].apply(lambda x: x**2 )
    df_1_frp_err = np.sqrt(top_df_1_clusters.groupby(['cluster_no'])['FRP_UNCERTAINTY_2'].agg('sum'))
    top_df_1_clusters.drop('FRP_UNCERTAINTY_2', axis=1,inplace=True)
    
    df_1_frp_sum = pd.DataFrame({'cluster_no':df_1_unique, 'summed_FRP': df_1_frp_sum.values[:],'FRP_uncertainty':df_1_frp_err.values[:],'mean_vza': df_1_vza_mean, 'pixel_count':df_1_counts , 'SAT': sat_1,'mean_lat':df_1_lat_mean,'mean_long':df_1_lon_mean })
    df_1_frp_sum['cluster_no'] = day+time +cluster_marker + df_1_frp_sum['cluster_no'].astype(str) ;df_1_frp_sum['cluster_no'] =  df_1_frp_sum['cluster_no'].astype(float)
    df_1_frp_sum['DATE_TIME'] = date_time_obj 
    df_1_frp_sum['SECONDS_SINCE'] = seconds_since
    df_1_frp_sum = df_1_frp_sum.reset_index(drop=True)
    
    df_2_frp_sum = top_df_2_matches.groupby(['cluster_no'])['FRP'].agg('sum')
    df_2_vza_mean = top_df_2_matches.groupby(['cluster_no'])['PIXEL_VZA'].mean()    
    df_2_lat_mean = top_df_2_matches.groupby(['cluster_no'])['LATITUDE'].mean()
    df_2_lon_mean = top_df_2_matches.groupby(['cluster_no'])['LONGITUDE'].mean()
    top_df_2_matches['FRP_UNCERTAINTY_2'] = top_df_2_matches['FRP_UNCERTAINTY'].apply(lambda x: x**2 )
    df_2_frp_err = np.sqrt(top_df_2_matches.groupby(['cluster_no'])['FRP_UNCERTAINTY_2'].agg('sum'))
    top_df_2_matches.drop('FRP_UNCERTAINTY_2', axis=1,inplace=True)
    df_2_frp_sum = pd.DataFrame({'cluster_no':df_2_unique, 'summed_FRP': df_2_frp_sum.values[:], 'FRP_uncertainty':df_2_frp_err.values[:],'mean_vza': df_2_vza_mean, 'pixel_count':df_2_counts, 'SAT': sat_2, 'mean_lat':df_2_lat_mean,'mean_long':df_2_lon_mean })
    df_2_frp_sum['cluster_no'] =   day+time+cluster_marker + df_2_frp_sum['cluster_no'].astype(str) ;df_2_frp_sum['cluster_no'] =  df_2_frp_sum['cluster_no'].astype(float)
    df_2_frp_sum['DATE_TIME'] = date_time_obj 
    df_2_frp_sum['SECONDS_SINCE'] = seconds_since
    df_2_frp_sum = df_2_frp_sum.reset_index(drop=True)

    return top_df_1_clusters, top_df_2_matches, df_1_frp_sum, df_2_frp_sum


def write_daily_fires(date,original_fires,matched_fires): # valud arguments 'M8' or 'M11'
    
    times = pd.date_range(start='10:00',end='10:00', freq='15min')
    times=times.format(formatter=lambda x: x.strftime('%H%M'))
    header_1=True
    header_2=True
    
    
    for time in times:
        cluster_marker = "8"
        date_time= date + time
        formated_dt = date_time[0:4]+'-'+date_time[4:6]+'-'+date_time[6:8]+'_'+ date_time[8:10]+':'+ date_time[10:12] +':00'
        m8_df,m11_df = extract_data(date,time)
        if m8_df.empty or m11_df.empty:
            print "ERROR while opening file"
            
        else:
            print "Finding and matching clusters " + time
            df_1 = m8_df ; df_2 = m11_df
            df_1_clusters = cluster_find(df_1,60)
            df_2_matches = cluster_match(df_1_clusters,df_2, original_fires ,matched_fires)
            df_m8_top_pixels,df_m11_pixel_matches, df1_m8_fires, df1_m11_fires = proccess_cluster_data(cluster_marker, df_1_clusters, df_2_matches,original_fires,matched_fires, date,time)

            cluster_marker = "9"
            df_1 = m11_df ; df_2 = m8_df
            df_1_clusters = cluster_find(df_1,60)
            df_2_matches = cluster_match(df_1_clusters,df_2, matched_fires,original_fires)
            df_m11_top_pixels,df_m8_pixel_matches, df2_m11_fires, df2_m8_fires = proccess_cluster_data(cluster_marker,df_1_clusters, df_2_matches,matched_fires,original_fires, date,time)
            
            df2_m11_fires['same_marker'] = 0 ; df2_m8_fires['same_marker'] = 0 ; df1_m11_fires['same_marker'] = 0 ; df1_m8_fires['same_marker'] = 0
            for i in range(0,len(df1_m8_fires.mean_lat),1):
                for j in range(0,len(df2_m11_fires.mean_lat),1):
                    lat1 =  df1_m8_fires['mean_lat'].loc[i] ; long1 =  df1_m8_fires['mean_long'].loc[i]
                    lat2 =  df2_m11_fires['mean_lat'].loc[j] ; long2 =  df2_m11_fires['mean_long'].loc[j]
                    point_distance = distance_on_unit_sphere(lat1,long1,lat2,long2)
                    if 0 <= point_distance < match_threshold:
                        df2_m11_fires.at[j,'same_marker'] = np.float(str(i)+str(j))
                        df2_m8_fires.at[j,'same_marker'] = np.float(str(i)+str(j))
                        df1_m11_fires.at[i,'same_marker'] = np.float(str(i)+str(j))
                        df1_m8_fires.at[i,'same_marker'] = np.float(str(i)+str(j))
                        
                        
            df2_unique_m11_fires = df2_m11_fires.loc[df2_m11_fires['same_marker'] == 0 ]  #; df2_m11_fires.drop('same_marker', axis=1, inplace=True)
            df2_unique_m8_fires = df2_m8_fires.loc[df2_m8_fires['same_marker'] == 0 ] #; df2_m8_fires.drop('same_marker', axis=1, inplace=True)
            df1_unique_m11_fires = df1_m11_fires.loc[df1_m11_fires['same_marker'] == 0 ]  #; df2_m11_fires.drop('same_marker', axis=1, inplace=True)
            df1_unique_m8_fires = df1_m8_fires.loc[df1_m8_fires['same_marker'] == 0 ] #; df2_m8_fires.drop('same_marker', axis=1, inplace=True)

            df2_match_m11_fires = df2_m11_fires.loc[df2_m11_fires['same_marker'] != 0 ]  #; df2_m11_fires.drop('same_marker', axis=1, inplace=True)
            df2_match_m8_fires = df2_m8_fires.loc[df2_m8_fires['same_marker'] != 0 ] #; df2_m8_fires.drop('same_marker', axis=1, inplace=True)
            df1_match_m11_fires = df1_m11_fires.loc[df1_m11_fires['same_marker'] != 0 ]  #; df2_m11_fires.drop('same_marker', axis=1, inplace=True)
            df1_match_m8_fires = df1_m8_fires.loc[df1_m8_fires['same_marker'] != 0 ] #; df2_m8_fires.drop('same_marker', axis=1, inplace=True)
            
            df2_match_m8_fires.sort_values(by=['same_marker'],ascending=True,inplace = True) ; df1_match_m11_fires.sort_values(by=['same_marker'],ascending=True,inplace = True)

            #print "This is the m8 matches:" ,df2_match_m8_fires, '\n',df1_match_m8_fires
            #print "This is the m11 matches:", df1_match_m11_fires,'\n',df2_match_m11_fires
            
            if len(df2_match_m8_fires) > len(df1_match_m8_fires):
                df2_match_m8_fires = df2_match_m8_fires[df2_match_m8_fires['same_marker'].isin(df1_match_m8_fires['same_marker'][:]) ]
                df2_match_m11_fires = df2_match_m11_fires[df2_match_m11_fires['same_marker'].isin(df1_match_m8_fires['same_marker'][:]) ]
            if len(df2_match_m8_fires) < len(df1_match_m8_fires):
                df1_match_m8_fires = df1_match_m8_fires[df1_match_m8_fires['same_marker'].isin(df2_match_m8_fires['same_marker'][:]) ]
                df1_match_m11_fires = df1_match_m11_fires[df1_match_m11_fires['same_marker'].isin(df2_match_m8_fires['same_marker'][:]) ]

            
            df2_match_m8_fires['cluster_no'][:] = df1_match_m11_fires['cluster_no'][:]
            df_all_top_fires = pd.concat([df2_match_m8_fires, df1_match_m11_fires,df2_unique_m11_fires,df2_unique_m8_fires,df1_match_m11_fires,df1_unique_m8_fires])
            df_all_top_fires = df_all_top_fires.sort_values(by=['cluster_no', 'SAT'],ascending=True).reset_index(drop=True)
            
            df_all_top_pixels = pd.concat([df_m8_top_pixels,df_m11_pixel_matches, df_m11_top_pixels, df_m8_pixel_matches])
            unique_fires = df_all_top_fires['cluster_no'][:]
            df_all_top_pixels = df_all_top_pixels[df_all_top_pixels['cluster_no'].isin(unique_fires) ]
            df_all_top_pixels = df_all_top_pixels.sort_values(by=['cluster_no', 'SAT'],ascending=True)
            df_all_top_pixels.reset_index(drop=True,inplace = True)
        
            mode = 'w' if header_1 else 'a'
            df_all_top_fires.to_csv("./txt_data/"+date+"_fires", mode=mode,header=header_1,index=False) 
            header_1=False
            mode = 'w' if header_2 else 'a'
            df_all_top_pixels.to_csv("./txt_data/"+date+"_pixels", mode=mode,header=header_2,index=False) 
            header_2=False


dates = ["20180601"]#,'20180602','20180603']#["20180515","20180516","20180517","20180518","20180519","20180520","20180521","20180522","20180523","20180524","20180525","20180526","20180527","20180528","20180529","20180530","20180531"]
for date in dates:
    print "NOTE: Generating fires for "  + date
    write_daily_fires(date, "8", "11")

##### need to fix ilimiting matching clusters!!!!!!


#plot_compare([pixels_m8,pixels_m11],date )


#slope, intercept, r_value, p_value, std_err = stats.linregress(m8_frp,m11_frp)

#line_start = slope*0 + intercept
#line_end = slope*1500 + intercept

#plot_scatter(m8_frp,m11_frp,m8_vza,m11_vza,date)

