import os
import h5py
import numpy as np
import pandas as pd
from time import gmtime, strftime


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
        pass                                                                       # pass if no concurrent pixels
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
        grid = pd.DataFrame({'GRID_NO':fire_ids, 'summed_FRP': frp_sum.values[:],'FRP_uncertainty':frp_err.values[:], 'pixel_count':pixel_count , 
                             'SAT': sat,'LAT_MIN':lat_min,'LAT_MAX':lat_max,'LAT_CENTER':lat_center, 'LONG_MIN':long_min,'LONG_MAX':long_max,'LONG_CENTER':long_center, 'TIME_WINDOW':time })
        #get grid_cell VZA from precomputed gridded VZA map
        if sat == 'MODIS':
            pass
        else:
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

"""
def write_gridded_file(lonmin,lonmax,latmin,latmax,m8_df,m11_df):   
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
    
    grid_df = pd.DataFrame(joint_boxes.str.split(' ',1).tolist(), columns = ['lat_box','long_box'])
    grid_df.lat_box  = (grid_df.lat_box.apply(lambda x: x.lstrip('0')).astype(int))-1
    grid_df.long_box  = (grid_df.long_box.apply(lambda x: x.lstrip('0')).astype(int))-1
    
    #regrid_vza(lonmin,lonmax,latmin,latmax)
    grid = pd.concat([m8_grid,m11_grid])    
    return grid
"""   


BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox" 
lon_max = 20
lon_min = -15
lat_max = -5
lat_min = -10
resolution = 1.0

print isinstance('bob', basestring)

t = np.array([datetime(2018,6,16,12,15),datetime(2018,6,16,12,15)])                                 # supply one off times
#t = np.arange(datetime(2018,6,18), datetime(2018,6,19), timedelta(minutes=15)).astype(datetime)
for k in range(0,len(t),1):


    modis_df = pd.read_csv(BaseDir + "/txt_data/MODIS/MODIS_1806160000_1806192300", parse_dates=True, index_col='TIME') 
    modis_df = modis_df[(modis_df['LATITUDE'] > lat_min) & (modis_df['LATITUDE'] < lat_max) ]                   #cut region
    modis_df = modis_df[(modis_df['LONGITUDE'] > lon_min) & (modis_df['LONGITUDE'] < lon_max) ]
    
    lat_bounds,lon_bounds = get_grid_bounds(lat_min,lat_max,lon_min,lon_max,resolution)   
    grid_df = get_gridded(t[k],modis_df,lat_bounds,lon_bounds, 'MODIS' ,resolution)
    
    
    



























