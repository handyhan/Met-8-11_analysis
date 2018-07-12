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
    

BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox" 
dates = ['20180610','20180611','20180612','20180613','20180616','20180617','20180618','20180619','20180621','20180624','20180625','20180626','20180630','20180701','20180702','20180704','20180705','20180707']

#dates = pd.date_range(start='2018-06-10',end='2018-06-10')
#dates = dates.format(formatter=lambda x: x.strftime('%Y%m%d'))
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
                df_fires.to_csv(BaseDir +"/txt_data/CLUSTERS_FD_fires_" +date, mode=mode,header=header_1,index=False) 
                header_1=False
                mode = 'w' if header_2 else 'a'
                df_pixels.to_csv(BaseDir +"/txt_data/CLUSTERS_FD_pixels_" +date, mode=mode,header=header_2,index=False) 
                header_2=False
              

  