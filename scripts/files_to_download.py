# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 10:26:12 2018

@author: Hannah.
"""
from datetime import date,timedelta
import csv
import pandas as pd
import paramiko
#from scp import SCPClient
import subprocess
import os


yesterday = date.today() - timedelta(1)
day_before =  date.today() - timedelta(2)
day_before_b = date.today() - timedelta(3)
yesterday=  yesterday.strftime('%Y%m%d')
day_before = day_before.strftime('%Y%m%d')
day_before_b = day_before_b.strftime('%Y%m%d')

times = pd.date_range(start='00:00',end='23:45', freq='15min')
times=times.format(formatter=lambda x: x.strftime('%H%M'))

download_list=[]
for time in times:
    time_yesterday =  yesterday + time 
    time_day_before =  day_before + time 
    time_day_before_b = day_before_b + time 
    download_list.extend([time_yesterday,time_day_before,time_day_before_b])
    
    
with open('C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/scripts/DownloadList.txt', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for i in download_list:
        writer.writerows([[i]])
        
      
        
#print "hello"
#subprocess.call("/home/mobaxterm/Documents/Earth_Observation/sandbox/Data/daily_download.sh", shell=True)
#print "hello"

"""

ssh = paramiko.SSHClient()
ssh.load_system_host_keys()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh = ssh.connect(hostname='193.137.20.100',port=22, username='rhnguyen', password="rh#9!yen")

#stdin, stdout, stderr = ssh.exec_command('ls -1 /home/rhnguyen|head -n 5')
#print "STDOUT:\n%s\n\nSTDERR:\n%s\n" %( stdout.read(), stderr.read() )

scp = SCPClient(ssh.get_transport())
"""