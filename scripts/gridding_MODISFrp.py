import os
import sys
import glob
import time
import datetime
import configparser
import numpy as np
from pyhdf.SD import SD
from scipy import interpolate
from netCDF4 import Dataset,num2date, date2num

dir_local = os.getcwd()
os.chdir(dir_local)
from write_file import create_NCfile

configFile = 'gridding.config'
conf = configparser.ConfigParser()
conf.read(configFile)

dir_modFRP = conf.get('source', 'dir_modFRP')
dir_mod14c = conf.get('source', 'dir_mod14c')
dir_mydFRP = conf.get('source', 'dir_mydFRP')
dir_myd14c = conf.get('source', 'dir_myd14c')
dir_tcwv   = conf.get('source', 'dir_tcwv')
#dir_npy    = conf.get('source', 'dir_npy')

dir_write  = conf.get('environment', 'dir_write')

gridsize = float(conf.get('gridding', 'gridsize'))
llat     = float(conf.get('gridding', 'llat'))
ulat     = float(conf.get('gridding', 'ulat'))
llon     = float(conf.get('gridding', 'llon'))
ulon     = float(conf.get('gridding', 'ulon'))
Nlat = int((ulat-llat)/gridsize)
Nlon = int((ulon-llon)/gridsize)
gridsetting = [gridsize,llat,ulat,Nlat,llon,ulon,Nlon]

starttime  = datetime.datetime.strptime(conf.get('run', 'starttime'),'%Y%m%d%H%M')
endtime    = datetime.datetime.strptime(conf.get('run', 'endtime'),'%Y%m%d%H%M')
interval   = datetime.timedelta(minutes=int(conf.get('run', 'interval')))

def get_frp_atm_correction(dir_tcwv,gridtime,modislist):
    def read_tcwv(dir_tcwv,gridtime):
        nc_file=os.path.join(dir_tcwv,'netcdf-global-20120101-20161231.nc')
        fh=Dataset(nc_file,mode='r')
        tcwvtimelist=fh.variables['time'][:]
        idx = (np.abs(tcwvtimelist-(gridtime-datetime.datetime(1900,1,1,0,0)).days*24)).argmin()
        tcwvArr = fh.variables['tcwv'][idx,:,:]
        fh.close()
        return (tcwvArr)
    def get_transmittance(tcwv,vzalist):
        vzalist = abs(vzalist)
        vangle = abs(np.deg2rad(vzalist))
        ### modis tcwv list
        H2O_list=[5,10,20,30,40,50,60]
        TAU_list=[0.1381,0.1399,0.1445,0.1504,0.1576,0.166, 0.1757]
        PO_list =[-0.0651622145,-0.0643401721,-0.0629781888,-0.0615412044,-0.0579979622,-0.0603463025,-0.0671102317]
        P1_list =[1.1533658383,1.1517193011,1.1485233606,1.1450919982,1.1382583085,1.1407634177,1.1498490808]
        P2_list =[-0.0739194221,-0.0732436675,-0.0718538574,-0.0703834389,-0.0675955496,-0.068377691,-0.0717979974]         
        tau = interpolate.UnivariateSpline(H2O_list, TAU_list)
        derivTau = tau.derivative()
        p0  = interpolate.UnivariateSpline(H2O_list, PO_list)
        derivP0 = p0.derivative()
        p1  = interpolate.UnivariateSpline(H2O_list, P1_list)
        derivP1 = p1.derivative()
        p2  = interpolate.UnivariateSpline(H2O_list, P2_list)
        derivP2 = p2.derivative()
        vzaRadians = np.deg2rad(vzalist)
        sumP = p0(tcwv) + p1(tcwv) * vzaRadians + p2(tcwv) * vzaRadians ** 2
        sumDerivP = derivP0(tcwv) + derivP1(tcwv) * vzaRadians + derivP2(tcwv) * vzaRadians ** 2
        cosSumP = np.cos(sumP)
        sinSumP = np.sin(sumP)
        transmittance = np.exp(-tau(tcwv) / cosSumP)
        eB = 10e-5 * transmittance * (710.51117 - 8.37751 * vzalist + 0.92238 * vzalist ** 2 -
                                      0.02525 * vzalist ** 3 + 0.00027 * vzalist ** 4)
        dTi = -(cosSumP * derivTau(tcwv) + tau(tcwv) * sinSumP * sumDerivP) / cosSumP ** 2
        eUtcwv = 0.24287 + 0.11172 * tcwv - 0.00090 * tcwv ** 2
        devATM = eB ** 2 / transmittance ** 2 + dTi ** 2 * eUtcwv ** 2
        return transmittance, devATM
    def get_FRP_dev(modislist,devATM):
        def get_radiance(brightnessTemp):
            wavelength = 3.95
            speedLight = 299792458.e0
            planckConst = 6.6260755e-34
            boltzConst = 1.380658e-23
            wavelength *= 1.e-6
            c1 = 2 * planckConst * speedLight ** 2  # [W.m2]
            c2 = planckConst * speedLight / boltzConst  # [K.m]
            radiance = c1 / (wavelength ** 5 * (np.exp(c2 / (wavelength * brightnessTemp)) - 1)) * 1.e-6
            return radiance
        frplist = modislist[:,0]
        t21list = modislist[:,7]
        madt21list  = modislist[:,5]
        meant21list = modislist[:,4]
        
        backgDiff = (get_radiance(t21list) -get_radiance(meant21list))
        backgStd = madt21list * 1.253314
        backgStd = (get_radiance(meant21list+backgStd) - get_radiance(meant21list-backgStd)) / 2.
        devBack = backgStd / backgDiff
        devNoise = np.zeros(t21list.shape)
        devNoise[t21list <= 331] = 0.0019
        devNoise[t21list > 331] = 0.0154
        devNoise /= backgDiff
        devFRP = np.sqrt((0.1 ** 2 + devATM ** 2 +
                             devNoise ** 2 + devBack ** 2).astype(float))*frplist
        return devFRP
        
    tcwvArr = read_tcwv(dir_tcwv,gridtime)
    frplist = modislist[:,0]
    latlist = modislist[:,1]
    lonlist = modislist[:,2]
    vzalist = modislist[:,3]
    it = (1.0*(90-latlist)/1.5).astype(int)
    jt = (1.0*(lonlist-(-180))/1.5).astype(int)
    tcwv = tcwvArr[it,jt]
    transmittance,devATM = get_transmittance(tcwv,vzalist)
    frplist/=transmittance
    devFRP = get_FRP_dev(modislist,devATM)
    modislist = np.concatenate((modislist,np.array([devFRP]).T),axis=1)
    return modislist
	
def grid_hour_modis(dir_grid,gridtime,gridsetting):
    def read_modis_sd(filename):
        product = SD(filename)
        Nfire = product.attributes()['FirePix']
        if Nfire==0:
            return np.zeros((1,8))
        datasets = ['FP_power', 'FP_latitude','FP_longitude',
                    'FP_ViewZenAng','FP_MeanT21','FP_MAD_T21','FP_sample', 'FP_T21']
        datasetList = []
        for sds in datasets:
            selection = product.select(sds).get()
            datasetList.append(selection)
        return (np.array(datasetList).transpose())
    
    gridsize,llat,ulat,Nlat,llon,ulon,Nlon = gridsetting
    frparr = np.zeros((Nlat,Nlon))
    numarr = np.zeros((Nlat,Nlon),dtype=int)
    vzaarr = np.zeros((Nlat,Nlon))
    stdarr = np.zeros((Nlat,Nlon))
    
    os.chdir(dir_grid)
    gridlist = glob.glob(gridtime.strftime('%Y/%j/*%Y%j.%H*.hdf'))
    if len(gridlist)==0:
        return (frparr,numarr,vzaarr,stdarr)
    for ind,filename in enumerate(gridlist):
        readlist = read_modis_sd(filename)
        if ind == 0:
            modislist = readlist
        else:
            modislist = np.concatenate((modislist,readlist))      
    modislist = modislist[modislist[:,0]>0,:]    
    if modislist.size==0:
        return (frparr,numarr,vzaarr,stdarr)
    modislist = get_frp_atm_correction(dir_tcwv,gridtime,modislist)
    frplist = modislist[:,0]
    vzalist = modislist[:,3]
    stdlist = modislist[:,-1]
    
    for ind in np.arange(len(frplist)):
        it = int((modislist[ind,1]-llat)/gridsize)
        jt = int((modislist[ind,2]-llon)/gridsize)    
        frparr[it,jt]+=frplist[ind]
        numarr[it,jt]+=1
        vzaarr[it,jt]+=vzalist[ind]
        stdarr[it,jt]+=stdlist[ind]**2
    vzaarr[numarr>0] = vzaarr[numarr>0]/numarr[numarr>0]
    stdarr[numarr>0] = (stdarr[numarr>0]/numarr[numarr>0])**0.5
    return (frparr,numarr,vzaarr,stdarr)

# ## for debugging
# timlist = np.arange(starttime,endtime,interval).astype(datetime.datetime)
# frparr,numarr,vzaarr,stdarr = grid_hour_modis(dir_modFRP,dir_mod14c,timlist[1],gridsetting)


######MAYBE don't need this####
### for MODIS TERRA
timlist = np.arange(starttime,endtime,interval).astype(datetime.datetime)
NCFilename = 'MOD_%s_%s'%(timlist[0].strftime('%y%m%d%H%M'),timlist[-1].strftime('%y%m%d%H%M'))

dataset,times,FRPARR,FRPnumARR,VZAARR,STDARR = create_NCfile(NCFilename,gridsetting)
times[:]  = date2num(timlist, units = times.units, calendar = times.calendar)



""""
for i,tim in enumerate(timlist):
    gridtime = tim
    frparr,numarr,vzaarr,stdarr = grid_hour_modis(dir_modFRP,gridtime,gridsetting)
    FRPARR[i,:,:]    = frparr
    FRPnumARR[i,:,:] = numarr
    VZAARR[i,:,:]    = vzaarr
    STDARR[i,:,:]    = stdarr
    print (gridtime,int(np.sum(frparr)),np.sum(numarr),np.sum(stdarr))
dataset.close()

### for MODIS AQUA
timlist = np.arange(starttime,endtime,interval).astype(datetime.datetime)
NCFilename = 'MYD_%s_%s'%(timlist[0].strftime('%y%m%d%H%M'),timlist[-1].strftime('%y%m%d%H%M'))

dataset,times,FRPARR,FRPnumARR,VZAARR,STDARR = create_NCfile(NCFilename,gridsetting)
times[:]  = date2num(timlist, units = times.units, calendar = times.calendar)

for i,tim in enumerate(timlist):
    gridtime = tim
    frparr,numarr,vzaarr,stdarr = grid_hour_modis(dir_modFRP,dir_mod14c,gridtime,gridsetting)
    FRPARR[i,:,:]    = frparr
    FRPnumARR[i,:,:] = numarr
    VZAARR[i,:,:]    = vzaarr
    STDARR[i,:,:]    = stdarr
    print (gridtime,int(np.sum(frparr)),np.sum(numarr),np.sum(stdarr))
dataset.close()

"""

