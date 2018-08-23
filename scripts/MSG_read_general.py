from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap, cm
import matplotlib as mpl
import matplotlib.colors as colors
import h5py
import numpy as np
import pandas as pd
import sys
import math
import os

pd.set_option('display.max_columns', 500)



def roundup(x):
    return int(math.ceil(x / 100.0)) * 100


# Get hdf FRP data into pandas #



def extract_data(datetime):
    BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox"
    msg = BaseDir +"/Data/SEVIRI/HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_MSG-Disk_" + datetime
    msg_QP = BaseDir + "/Data/SEVIRI/HDF5_LSASAF_MSG_FRP-PIXEL-QualityProduct_MSG-Disk_" + datetime
    
    if  os.path.exists(msg):
        msg_QP = h5py.File(msg_QP,'r+') 
    if os.path.exists(msg):
        msg_f = h5py.File(msg,'r+') 
        msg_d = {'LATITUDE':msg_f['LATITUDE'][:],'LONGITUDE':msg_f['LONGITUDE'][:],'FRP':msg_f['FRP'][:],'FRP_UNCERTAINTY':msg_f['FRP_UNCERTAINTY'][:],'PIXEL_SIZE':msg_f['PIXEL_SIZE'][:],'PIXEL_VZA':msg_f['PIXEL_VZA'][:],'PIXEL_LINE':msg_f['ABS_LINE'][:],'PIXEL_COL':msg_f['ABS_PIXEL'][:]}
        msg_df = pd.DataFrame(data=msg_d,dtype=np.float32)
        msg_df['TIME']=datetime
        msg_df['TIME'] = pd.to_datetime(msg_df['TIME'],format="%Y%m%d%H%M") 
        cols = ['TIME','LATITUDE','LONGITUDE','FRP','FRP_UNCERTAINTY','PIXEL_SIZE','PIXEL_VZA','PIXEL_LINE','PIXEL_COL']
        msg_df = msg_df[cols]    
        # account for scaling factors ect #
        msg_df['LATITUDE'] = msg_df['LATITUDE']/100
        msg_df['LONGITUDE'] = msg_df['LONGITUDE']/100
        msg_df['FRP']/10
        msg_df['FRP_UNCERTAINTY'] = msg_df['FRP_UNCERTAINTY']/100
        msg_df['PIXEL_SIZE'] = msg_df['PIXEL_SIZE']/100
        msg_df['PIXEL_VZA'] = msg_df['PIXEL_VZA']/100
        msg_df['AQUIRE_SEC'] = np.round((900.0/3712)*(msg_df['PIXEL_LINE']))           # plus the time for the scan to reach that line,col
        msg_df['TIME'] = msg_df['TIME'] + pd.to_timedelta(msg_df['AQUIRE_SEC'], unit='s')
        msg_df.drop(columns=['AQUIRE_SEC'],inplace=True)
        
        # get rid of saturated pixels    
        drop_list=[]
        for i in range(0,len(msg_df['PIXEL_LINE'] -1),1):
            pix_col = int(msg_df['PIXEL_COL'][i]) -1
            pix_LINE = int(msg_df['PIXEL_LINE'][i]) - 1
            Q_flag =  msg_QP['QUALITYFLAG'][pix_LINE][pix_col]
            if Q_flag == 2 :
                drop_list.append(i)
        msg_df = msg_df.drop(msg_df.index[drop_list])        
        msg_df = msg_df.reset_index(drop=True)
        
        return msg_df
    
    elif not os.path.exists(msg):
        print "ERROR when trying to open file, this file does not exist:" + msg    
        msg_df = [pd.DataFrame()]
        return msg_df



def MSG_to_txt(dates):       # see extract data function in MSG_read to define data dirs
    ## generate daily files of FRP-PIXEL product (met-8 & 11), excluding saturated pixels
    for date in dates:
        print "NOTE: Generating fires for "  + date
        times = pd.date_range(start='00:00',end='23:45', freq='15min')
        times=times.format(formatter=lambda x: x.strftime('%H%M'))
        header_1=True    
        for time in times:
            print "Pulling data for " + time
            msg_df = extract_data(datetime)
            
            if m11_df.empty: #from here just operates on met-11
                print "ERROR while opening file"
            else:
                mode = 'w' if header_1 else 'a'
                m11_df.to_csv(BaseDir +"/txt_data/SEVIRI/FRP_PIXELS_" +date, mode=mode,header=header_1,index=False) 
                header_1=False
              

def extract_QP(datetime):
    datetime_str = datetime.strftime('%Y%m%d%H%M')
    BaseDir = "C:/Users/Hannah.N/Documents/Earth_Observation/sandbox"
    m8_QP = BaseDir +"/Data/SEVIRI/NetCDF_LSASAF_MSG-IODC_FRP-PIXEL-QualityProduct_IODC-Disk_" + datetime_str
    m11_QP = BaseDir + "/Data/SEVIRI/HDF5_LSASAF_MSG_FRP-PIXEL-QualityProduct_MSG-Disk_" + datetime_str
    #dir_QP = 
    
    
    def get_QM(col,line,h5_QP):
        return m11_QP['QUALITYFLAG'][line][col]
    
    if  os.path.exists(m11_QP):
        #m8_QP = h5py.File(m8_QP,'r+')
        m11_QP = h5py.File(m11_QP,'r+') 
        
        lines = np.arange(0,3712,1)
        lines = np.tile(lines,3712)
        columns = np.arange(0,3712,1)
        columns = np.repeat(columns,3712)
        
        QP_df = pd.DataFrame({'LINE':lines,'COL':columns})
        
        for index, row in QP_df.iterrows():                
            Qproduct = np.array([get_QM(row['LINE'],row['COL'],m11_QP)])
            if index == 0:
                Qproduct_full = Qproduct
            else:
                Qproduct_full = np.concatenate((Qproduct_full,Qproduct)) 
                
        print Qproduct_full
        #value =  m11_QP['QUALITYFLAG'][3711][4]
        
        #print value


"""
dates = pd.date_range(start='2018-06-09',end='2018-06-26')
dates = dates.format(formatter=lambda x: x.strftime('%Y%m%d'))
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
            mode = 'w' if header_1 else 'a'
            m8_df.to_csv("./txt_data/"+date+"_m11_all_pixels", mode=mode,header=header_1,index=False) 
            header_1=False
            mode = 'w' if header_2 else 'a'
            m11_df.to_csv("./txt_data/"+date+"_m8_all_pixels", mode=mode,header=header_2,index=False) 
            header_2=False
"""            
     