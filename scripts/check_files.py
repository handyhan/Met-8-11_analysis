
from datetime import date,timedelta
import csv
import pandas as pd
import os.path
#from scp import SCPClient
import subprocess
import os




dates = pd.date_range(start='2018-06-09',end='2018-06-28')
dates = dates.format(formatter=lambda x: x.strftime('%Y%m%d'))

times = pd.date_range(start='00:00',end='23:45', freq='15min')
times=times.format(formatter=lambda x: x.strftime('%H%M'))

missing_list = []
for date in dates:
    n = 0
    print  date
    for time in times:
        m8 = "./NetCDF_LSASAF_MSG-IODC_FRP-PIXEL-ListProduct_IODC-Disk_" + date+time
        m11 = "./HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_MSG-Disk_" + date +time
        m8_QP = "./NetCDF_LSASAF_MSG-IODC_FRP-PIXEL-QualityProduct_IODC-Disk_" + date+time
        m11_QP = "./HDF5_LSASAF_MSG_FRP-PIXEL-QualityProduct_MSG-Disk_" + date +time
        products = [m8,m11,m11_QP,m8_QP ]
        for  product in products:
            
            if os.path.exists(product):
                pass
               
            else:
                n = n+1
                print "no file"

    if n > 0:
        missing_list.extend([date])
         

with open('./CheckList.txt', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for i in missing_list:
        writer.writerows([[i]])
                   
        
      