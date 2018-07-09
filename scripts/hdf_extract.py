import numpy as np
import pandas as pd
import sys
import math
import os
import h5py

lon = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/SEVIRI_Static/HDF5_LSASAF_MSG_LON_IODC-Disk_201711300000"
lat = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/SEVIRI_Static/HDF5_LSASAF_MSG_LAT_IODC-Disk_201711300000"
lon = h5py.File(lon,'r')
lat = h5py.File(lat,'r')
lon = lon[u'LON']
lat = lat[u'LAT']


row = np.arange(0,3712,1)
col = np.arange(0,3712,1)
steps_lat = len(row)
steps_long = len(col)  
    
row_ix = np.tile(row[:],steps_long)
col_ix = np.repeat(col[:],steps_long)



disk = pd.DataFrame({'row_ix':row_ix,'col_ix':col_ix,'LATITUDE':np.nan,'LONGITUDE':np.nan})

def find_lat_lon(df):
    #df['LATITUDE'] =
    for row in df.itertuples():   
        
        i = getattr(row, "row_ix")
        j = getattr(row, "col_ix")
        latitude = lat[i][j] 
        longitude = lon[i][j]
        if (latitude == 9100) or (longitude == 9100.0):
            pass
        else:
            df.at[row.Index,'LATITUDE'] = latitude
            df.at[row.Index, 'LONGITUDE'] = longitude

    return df
     

disk = find_lat_lon(disk)


print disk

disk.to_csv("C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/txt_data/TEST_MSG8_lat_lon_", header=True,index=False)

"""
latitudes = np.array([[0,0,0,0]])

section = [500,1000,1500,2000,2500,3000]
for k in range(0,6,1):
    
    for i in range(section[k],section[(k+1)],1):
        for j in range(0,3712,1):
            if lat[i][j] & lon[i][j]!= 9100 :
                latitudes = np.append(latitudes, [[lat[i][j],lon[i][j],i,j]], axis = 0 )
                
    df = pd.DataFrame(data = latitudes[1:,:], columns = ['LATITUDE', 'LONGITUDE','I_val', 'J_val'] )
    df.to_csv("C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/txt_data/MSG8_lat_lon_" + str(section[k]), header=True,index=False)
"""