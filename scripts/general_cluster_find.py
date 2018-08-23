import numpy as np
import pandas as pd
import math
from itertools import product
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
        pass
    
    else:
        # get range of line and cols (for ease of gridding) prepare pixels for banary mask
        subtract_factor = min(sub_satellite_df[line_label].min(),sub_satellite_df[column_label].min())  
        sub_satellite_df['mask_line'] = sub_satellite_df[line_label] - subtract_factor
        sub_satellite_df['mask_col'] = sub_satellite_df[column_label] - subtract_factor
        grid_dims = int(math.ceil((max(sub_satellite_df['mask_line'].max(),sub_satellite_df['mask_col'].max())) / 100.0) * 100.0)
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
        sub_satellite_df = sub_satellite_df.sort_values(by=['cluster_no'])#.drop(columns = ['mask_line','mask_col'])
        return sub_satellite_df



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
    
    
        