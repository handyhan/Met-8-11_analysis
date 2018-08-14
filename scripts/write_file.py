import os
import sys
import time
import configparser
import numpy as np
from netCDF4 import Dataset,num2date, date2num

os.chdir('C:/Users/Tianran/Dropbox/Research/CAMS_GFAS/VZA/gridding')
configFile = 'gridding.config'
conf = configparser.ConfigParser()
conf.read(configFile)
dir_write  = conf.get('environment', 'dir_write')


def create_NCfile(NCFilename,gridsetting):
    os.chdir(dir_write)
    dataset = Dataset('%s.nc'%(NCFilename),'w',formate='NETCDF4_CLASSIC') 
    gridsize,llat,ulat,Nlat,llon,ulon,Nlon = gridsetting
    ###create dimensions for dataset###
    datasetLat = dataset.createDimension('lat',Nlat)
    datasetLon = dataset.createDimension('lon',Nlon)
    datasettim = dataset.createDimension('time',None)
    ###create variables for dataset###
    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float32,('lat',))
    longitudes = dataset.createVariable('longitude', np.float32,('lon',)) 
    FRPARR = dataset.createVariable('FRP', np.float32,('time','lat','lon'))
    FRPnumARR = dataset.createVariable('FRPnum', np.int,('time','lat','lon'))
    VZAARR = dataset.createVariable('VZA', np.float32,('time','lat','lon'))
    STDARR = dataset.createVariable('STD', np.float32,('time','lat','lon'))
    ###create attributes for dataset###
    dataset.description = 'gridding FRP product'
    dataset.history = 'Created ' + time.ctime(time.time())
    dataset.source = 'please contact: tianran.zhang@kcl.ac.uk'
    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    times.units = 'hours since 0001-01-01 00:00:00'
    times.calendar = 'gregorian'
    FRPARR.units = 'MW'
    FRPnumARR.units = 'int number'
    VZAARR.units = 'degree'
    STDARR.units = 'MW'
    latitudes[:] = np.arange(llat,ulat,gridsize)+gridsize/2
    longitudes[:] = np.arange(llon,ulon,gridsize)+gridsize/2
    return (dataset,times,FRPARR,FRPnumARR,VZAARR,STDARR)