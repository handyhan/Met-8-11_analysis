import numpy as np
import pandas as pd
import math
from itertools import product
from Omission_comission_V2 import neighbours3x3
import matplotlib.pyplot as plt
#pd.set_option('display.max_columns', 500)
#pd.options.mode.chained_assignment = None  

#  takes empty grid to assign identities, banary grid of pixels and the grid dimentions
def find_clusters(cluster_grid,grid,grid_dims):       
    n=1
    equivalance = []
    #First pass, scan rows and assign to any left or top nieghbours the same cluster 
    # keep a list of equivalent clusters when a conflict occurs
    print "Running first pass"
    for  x,y  in product(range(grid_dims), range(grid_dims)):
        if grid[x, y] == 0:   #pass if no AF pix
            pass            
        else:
            if y == 0 and x > 0 and grid[x-1, y] == 1 :
                cluster_grid[x, y] = cluster_grid[x-1, y]
                #print "first top"
            elif x == 0 and y > 0 and grid[x, y-1] == 1 :
                cluster_grid[x, y] = cluster_grid[x, y-1]
                #print "first left"
            elif x > 0 and y > 0 and grid[x-1, y] == 1 :
                top = cluster_grid[x-1, y]
                left = cluster_grid[x, y-1]
                if grid[x, y-1] == 0:
                    #print "top"
                    cluster_grid[x, y] = top
                elif grid[x, y-1] == 1:
                    #print "conflict top"
                    cluster_grid[x, y] = top
                    equivalance.append([left,top])
            elif x > 0 and y > 0 and grid[x, y-1] == 1 :
                #print "left"
                cluster_grid[x, y] = cluster_grid[x, y-1]
            elif x > 0 and y > 0 and grid[x-1, y-1] == 1 :
                #print "left diag"
                cluster_grid[x, y] = cluster_grid[x-1, y-1]   
            else: 
                #print "new gorup"
                cluster_grid[x, y] = n 
                n = n+1   
    #second pass, going from high to low re-assign cluster numbers to clusters in which there was a conflict
    equivalance = sorted(equivalance)
    print "Running second pass"
    for i in reversed(range(0,len(equivalance),1)):
        for  x,y  in product(range(grid_dims), range(grid_dims)):
            if cluster_grid[x, y] == 0:
                pass
            elif cluster_grid[x, y] == equivalance[i][0]:
                #print "updated"                    
                cluster_grid[x, y] = equivalance[i][1]
    return cluster_grid
    

# takes a datetime object, panda data frame with its names for the line and column names (str)
def get_clusters_at_time(time, satellite_df, line_label, column_label):    
    #select pixels only in time range
    start_time = time - pd.Timedelta(minutes=8)
    end_time = time + pd.Timedelta(minutes=8)    
    sub_satellite_df = satellite_df.loc[start_time:end_time]                        
    if sub_satellite_df.empty:
        print "No Concurent pixels at time " + str(time)   
        df_empty = pd.DataFrame()
        return df_empty
    else:
        # get range of line and cols (for ease of gridding) prepare pixels for re-grid to smaller binary mask
        subtract_factor = min(sub_satellite_df[line_label].min(),sub_satellite_df[column_label].min())  
        sub_satellite_df['mask_line'] = sub_satellite_df[line_label] - subtract_factor
        sub_satellite_df['mask_col'] = sub_satellite_df[column_label] - subtract_factor
        grid_dims = int(math.ceil((max(sub_satellite_df['mask_line'].max(),sub_satellite_df['mask_col'].max())) / 100.0) * 100.0) + 1
        # linked component clustering 
        cluster_grid = np.zeros((grid_dims,grid_dims))                          # emplt matrix to fill with the cluster identity
        grid = np.zeros((grid_dims,grid_dims))                                    
        for index, row in sub_satellite_df.iterrows():
            grid[int(row['mask_col']),int(row['mask_line'])] = 1                # make binary mask of where AF pixels are   
        # get identaty grid
        cluster_grid = find_clusters(cluster_grid,grid,grid_dims)               # apply linked component algo to get independant fire clusters
        cluster_grid = cluster_grid.flatten('F')                                # flatten column wise
        grid_col = np.array(range(0,grid_dims,1))                               # make arrays of the column and line index
        grid_col = np.tile(grid_col,grid_dims)
        grid_line = np.array(range(0,grid_dims,1))
        grid_line = np.repeat(grid_line,grid_dims)
        
        clusters_df = pd.DataFrame({'mask_col':grid_col,'mask_line':grid_line,'cluster_no':cluster_grid})   # match index to correct cluster_#
        clusters_df = clusters_df[clusters_df['cluster_no'] > 0 ].sort_values('cluster_no')        
        clusters_df[line_label] = clusters_df['mask_line'] + subtract_factor                                # revert back to original columns and rows
        clusters_df[column_label] = clusters_df['mask_col'] + subtract_factor
        clusters_df = clusters_df.sort_values(by=['mask_line','mask_col'])
        sub_satellite_df = sub_satellite_df.sort_values(by=['mask_line','mask_col'])
        sub_satellite_df['cluster_no']  = clusters_df['cluster_no'].values[:]                               # add identities to the pixels 
        sub_satellite_df = sub_satellite_df.sort_values(by=['cluster_no']).drop(columns = ['mask_line','mask_col'])
        
        #print sub_satellite_df
        return sub_satellite_df




def cluster_matchup(date_range, modis_df, msg_df):
    first=True    
    for k in range(0,len(date_range),1):     
        time = date_range[k]   
        time_str = time.strftime("%Y%m%d%H%M")
        print time_str
        start_time = time - pd.Timedelta(minutes=8)
        end_time = time + pd.Timedelta(minutes=8)    
        sub_mod = modis_df.loc[start_time:end_time]
        sub_msg = msg_df.loc[start_time:end_time]
        
        if (sub_mod.empty) or (sub_msg.empty):
            #print "No Concurent pixels at time " + str(time)                      ## pass if no concurrent pixels
            print "No Concurent fires at time " + str(time)                
            sub_msg_processed = pd.DataFrame()
            mod_omissions = pd.DataFrame()
            mod_match= pd.DataFrame()
            #return sub_msg_processed,mod_omissions,mod_match
        else:                    
            #print sub_mod
            
            sub_mod['MSG_TAG'] = sub_mod['MSG_COL']*10000 + sub_mod['MSG_LINE']                          # combine line and colum to get unique tag for seviri pixels
            sub_msg['MSG_TAG'] = sub_msg['PIXEL_COL']*10000 + sub_msg['PIXEL_LINE']       
            
            
            unique_mod_pix = np.unique(sub_mod['MSG_TAG'])                                           # find grid cells at this time unique to both MODIS and MSG
            matching_msg_pix = sub_msg[sub_msg['MSG_TAG'].isin(unique_mod_pix)]
            unique_msg_match_pix = np.unique(matching_msg_pix['MSG_TAG'])
            mod_pix_full_match = sub_mod[sub_mod['MSG_TAG'].isin(unique_msg_match_pix)]
            unique_full_mod_pix = np.unique(mod_pix_full_match['MSG_TAG'])
            msg_pix_full_match = matching_msg_pix[matching_msg_pix['MSG_TAG'].isin(unique_full_mod_pix)]  
            mod_pix_full_match = mod_pix_full_match.sort_values('MSG_TAG')
            msg_pix_full_match = msg_pix_full_match.sort_values('MSG_TAG')
            mod_pix_full_match['MSG_CLUSTER'] = mod_pix_full_match['MSG_TAG'].apply(lambda x: msg_pix_full_match[msg_pix_full_match['MSG_TAG'] == x]['cluster_no'][0])
            msg_pix_full_match['MODIS_CLUSTER'] = msg_pix_full_match['MSG_TAG'].apply(lambda x: mod_pix_full_match[mod_pix_full_match['MSG_TAG'] == x]['cluster_no'][0])
            mod_pix_full_match = mod_pix_full_match.sort_values('MSG_CLUSTER')                      #grouping by MSG cluster as it is bigge and more likely to contain multiple MODIS clusters            
            msg_pix_full_match = msg_pix_full_match.sort_values('cluster_no')
            
            
            unique, counts = np.unique(msg_pix_full_match['cluster_no'], return_counts=True) 
            
            mod_frp_sum = mod_pix_full_match.groupby(['MSG_CLUSTER'])['FRP'].agg('sum')
            #mod_vza_mean = mod_pix_full_match.groupby(['MSG_CLUSTER'])['PIXEL_VZA'].mean()
            mod_lat_mean = mod_pix_full_match.groupby(['MSG_CLUSTER'])['LATITUDE'].mean()
            mod_lon_mean = mod_pix_full_match.groupby(['MSG_CLUSTER'])['LONGITUDE'].mean()
            
            msg_frp_sum = msg_pix_full_match.groupby(['cluster_no'])['FRP'].agg('sum')
            #msg_vza_mean = msg_pix_full_match.groupby(['cluster_no'])['PIXEL_VZA'].mean()            
            msg_lat_mean = mod_pix_full_match.groupby(['MSG_CLUSTER'])['LATITUDE'].mean()
            msg_lon_mean = msg_pix_full_match.groupby(['cluster_no'])['LONGITUDE'].mean()
            
            fires_match = pd.DataFrame({'cluster_no':unique, 'MODIS_FRP': mod_frp_sum.values[:],'MODIS_LAT':mod_lat_mean.values[:],'MODIS_LONG':mod_lon_mean.values[:],'MSG_FRP': msg_frp_sum.values[:],'MSG_LAT':msg_lat_mean.values[:],'MSG_LONG':msg_lon_mean.values[:] })        
            fires_match['cluster_no'] = time_str + fires_match['cluster_no'].astype(str)   # make time unique identifier

            if first == True:
                all_fire_match = fires_match
            else:                                                                                            # append all matches from all times together
                all_fire_match = all_fire_match.append(fires_match)                           
            first = False
            
    print all_fire_match
                    
            #grouping by MSG cluster as it is bigge and more likely to contain multiple MODIS clusters
            

    """
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
    """
                    
          
BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox" 
lon_max = 20
lon_min = -15
lat_max = -5
lat_min = -10
resolution = 1.0

    
modis_df = pd.read_csv(BaseDir + "/TEST_DATA/MODIS_fires_test", parse_dates=True,index_col='TIME').sort_index()
msg_df = pd.read_csv(BaseDir + "/TEST_DATA/MSG_fires_test", parse_dates=True,index_col='TIME').sort_index()
t = np.arange(datetime(2014,8,1), datetime(2014,8,2), timedelta(minutes=15)).astype(datetime)

cluster_matchup(t, modis_df,msg_df )






"""
cluster_group = sub_mod[sub_mod['cluster_no'] == unique_mod_clusters]
unique_mod_msg_tags = np.unique(cluster_group['MSG_TAG'])

matching_msg = sub_msg[sub_msg['MSG_TAG'].isin(unique_mod_msg_tags)]
print matching_msg

for index, row in sub_msg.iterrows():                    
neighbour_list = np.array(neighbours3x3(row['MSG_COL'],row['MSG_LINE']))
if index == 0:
    neighbour_list_full = neighbour_list
else:
    neighbour_list_full = np.concatenate((neighbour_list_full,neighbour_list),axis=1)           # put all into array

neighbour_list_full = pd.DataFrame(data=neighbour_list_full.T,columns = ['PIXEL_COL','PIXEL_LINE'])
neighbour_list_full['MSG_TAG'] = neighbour_list_full['PIXEL_COL']*10000 + neighbour_list_full['PIXEL_LINE']
neighbour_list_full= neighbour_list_full.drop_duplicates('MSG_TAG')
"""
#unique_mod_clusters = np.unique(neighbour_list_full['MSG_TAG'])


            
                
            ##do neighbour find
    
    






















"""

BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox" 
lon_max = 20
lon_min = -15
lat_max = -5
lat_min = -10

t = np.array([datetime(2018,6,16,12,15),datetime(2018,6,16,12,15)])                                 # supply one off times
#t = np.arange(datetime(2018,6,18), datetime(2018,6,19), timedelta(minutes=15)).astype(datetime)
for k in range(0,len(t),1):
    
    # read in txt data with columns  FRP,LATITUDE,LONGITUDE,VZA,FP_sample,FP_line,SAT,TIME,FRP_UNCERTAINTY
    # NOTE: time must be in format eg. 2018-06-12 09:51:34, lat,long headers may vary (re-name) line/column headers can vary
    modis_df = pd.read_csv(BaseDir + "/txt_data/MODIS/MODIS_1806160000_1806192300", parse_dates=True,index_col='TIME') 
    modis_df = modis_df[(modis_df['LATITUDE'] > lat_min) & (modis_df['LATITUDE'] < lat_max) ]                   #cut region
    modis_df = modis_df[(modis_df['LONGITUDE'] > lon_min) & (modis_df['LONGITUDE'] < lon_max) ]
    
    sub_mod = get_clusters_at_time(t[k], modis_df, 'FP_line','FP_sample')           # run clustering algo, may take some minutes
    plt.scatter(sub_mod['FP_line'].values[:],sub_mod['FP_sample'].values[:],c=sub_mod['cluster_no'].values[:])
    plt.show()
    
    
"""     