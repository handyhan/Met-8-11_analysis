from pyhdf.SD import SD,SDC
from datetime import datetime



DataDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/Data/MODIS"
fname= "MYD14.A2018161.1120.006.2018162115317.hdf"

sat = fname[0:3]
start_time = fname[15:19]
start_year = int(fname[7:11])
start_day = int(fname[11:14])



def read_MODIS_file(fname):
    os.chdir(DataDir)
    hdf = SD(fname,SDC.READ)
    metadata=hdf.attributes()
    firepix=metadata['FirePix']
    if firepix == 0:
        print ('no fire pixel in %s'%(fname))
        hdf.end()
        return []
    print ('%s fire pixel in %s'%(firepix,fname))
    data = np.zeros((firepix,5))
    sds = hdf.select('FP_latitude')
    latdata = sds.get()
    sds = hdf.select('FP_longitude')
    londata = sds.get()
    sds = hdf.select('FP_power')
    FRPdata = sds.get()
    sds = hdf.select('FP_confidence')
    FRPcondata = sds.get()
    sds = hdf.select('FP_ViewZenAng')
    vzadata = sds.get()
    data[:,0] = latdata
    data[:,1] = londata
    data[:,2] = FRPdata
    data[:,3] = vzadata
    data[:,4] = FRPcondata
    
    hdf.end()
    return data


data_mod = read_MODIS_file(fname)

start_date = datetime(start_year, 1, 1) + timedelta(start_day - 1)
start_date = start_date.replace(hour= int(start_time[0:2]), minute=int(start_time[2:4]))
end_date = start_date + timedelta(minutes=5)

df_pixel = pd.DataFrame(data=data_mod, columns = ['LATITUDE','LONGITUDE','FRP','FRP_CONFIDENCE','PIXEL_VZA'])
df_pixel['SAT']= sat
df_pixel['START_DATE_TIME'] = start_date
df_pixel['END_DATE_TIME'] = end_date

print df_pixel



