from matplotlib import rcParams
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import matplotlib as mpl
import matplotlib.colors as colors
from MSG_read import extract_data,plot_compare,plot_scatter_fit, plot_stats,plot_single_map,plot_ODR_scatter_residuals,plot_scatter_residuals
import h5py
import numpy as np
import pandas as pd
import sys
import math
import os
from math import sin, cos, sqrt, atan2, radians,atan, degrees


lon = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/SEVIRI_Static/HDF5_LSASAF_MSG_LON_IODC-Disk_201711300000"
lat = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/SEVIRI_Static/HDF5_LSASAF_MSG_LAT_IODC-Disk_201711300000"
lon = h5py.File(lon,'r')
lat = h5py.File(lat,'r')
lon = lon[u'LON']
lat = lat[u'LAT']
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
df = pd.read_csv("./txt_data/MSG_lat_lon_1000")
#df_2 = pd.read_csv("./txt_data/MSG_lat_lon_2000")
#df_3 = pd.read_csv("./txt_data/MSG_lat_lon_3000")
#df_4 = pd.read_csv("./txt_data/MSG_lat_lon_4000")

#df = pd.concat([df_1,df_2,df_3,df_4])
print df

df['LATITUDE'] = df['LATITUDE']/100
df['LONGITUDE'] = df['LONGITUDE']/100

print df
lat = df['LATITUDE'].values[:]
lon = df['LONGITUDE'][:].values[:]

plot_single_map(lat,lon,'test','disk')

"""