from matplotlib import rcParams
import matplotlib.pyplot as plt
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



def extract_data(date, time):
    #print date
    m8 = "./Data/NetCDF_LSASAF_MSG-IODC_FRP-PIXEL-ListProduct_IODC-Disk_" + date+time
    m11 = "./Data/HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_MSG-Disk_" + date +time
    m8_QP = "./Data/NetCDF_LSASAF_MSG-IODC_FRP-PIXEL-QualityProduct_IODC-Disk_" + date+time
    m11_QP = "./Data/HDF5_LSASAF_MSG_FRP-PIXEL-QualityProduct_MSG-Disk_" + date +time
    
    
    if os.path.exists(m8) and os.path.exists(m11):
        m8_QP = h5py.File(m8_QP,'r+')
        m11_QP = h5py.File(m11_QP,'r+') 
    if os.path.exists(m8) and os.path.exists(m11):
        m8_f = h5py.File(m8,'r+')
        m11_f = h5py.File(m11,'r+') 
        m8_d = {'LATITUDE':m8_f['LATITUDE'][:],'LONGITUDE':m8_f['LONGITUDE'][:],'FRP':m8_f['FRP'][:],'FRP_UNCERTAINTY':m8_f['FRP_UNCERTAINTY'][:],'PIXEL_SIZE':m8_f['PIXEL_SIZE'][:],'PIXEL_VZA':m8_f['PIXEL_VZA'][:],'PIXEL_ROW':m8_f['ABS_LINE'][:],'PIXEL_COL':m8_f['ABS_PIXEL'][:]}
        m11_d = {'LATITUDE':m11_f['LATITUDE'][:],'LONGITUDE':m11_f['LONGITUDE'][:],'FRP':m11_f['FRP'][:],'FRP_UNCERTAINTY':m11_f['FRP_UNCERTAINTY'][:],'PIXEL_SIZE':m11_f['PIXEL_SIZE'][:],'PIXEL_VZA':m11_f['PIXEL_VZA'][:],'PIXEL_ROW':m11_f['ABS_LINE'][:],'PIXEL_COL':m11_f['ABS_PIXEL'][:]}
        m8_df = pd.DataFrame(data=m8_d,dtype=np.float32)
        m11_df = pd.DataFrame(data=m11_d,dtype=np.float32)
        m8_df['DATE']=date ; m11_df['DATE']=date; m8_df['TIME']=time; m11_df['TIME']=time
        cols = ['DATE','TIME','LATITUDE','LONGITUDE','FRP','FRP_UNCERTAINTY','PIXEL_SIZE','PIXEL_VZA','PIXEL_ROW','PIXEL_COL']
        m8_df = m8_df[cols]
        m11_df = m11_df[cols]    
        # account for scaling factors ect #
        m8_df['LATITUDE'] = m8_df['LATITUDE']/100 ; m11_df['LATITUDE'] = m11_df['LATITUDE']/100
        m8_df['LONGITUDE'] = m8_df['LONGITUDE']/100 ; m11_df['LONGITUDE'] = m11_df['LONGITUDE']/100
        m8_df['FRP'] = m8_df['FRP']/10 ; m11_df['FRP'] = m11_df['FRP']/10
        m8_df['FRP_UNCERTAINTY'] = m8_df['FRP_UNCERTAINTY']/100 ; m11_df['FRP_UNCERTAINTY'] = m11_df['FRP_UNCERTAINTY']/100
        m8_df['PIXEL_SIZE'] = m8_df['PIXEL_SIZE']/100 ; m11_df['PIXEL_SIZE'] = m11_df['PIXEL_SIZE']/100
        m8_df['PIXEL_VZA'] = m8_df['PIXEL_VZA']/100 ; m11_df['PIXEL_VZA'] = m11_df['PIXEL_VZA']/100
        # get rid of saturated pixels
        
        drop_list=[]
        for i in range(0,len(m8_df['PIXEL_ROW'] -1),1):
            pix_col = int(m8_df['PIXEL_COL'][i]) -1
            pix_row = int(m8_df['PIXEL_ROW'][i] - 1)
            Q_flag =  m8_QP['QUALITYFLAG'][pix_row][pix_col]
            if Q_flag == 2 :
                drop_list.append(i)
        m8_df = m8_df.drop(m8_df.index[drop_list])        
        m8_df = m8_df.reset_index(drop=True)
        drop_list=[]
        for i in range(0,len(m11_df['PIXEL_ROW'] -1),1):
            pix_col = int(m11_df['PIXEL_COL'][i]) -1
            pix_row = int(m11_df['PIXEL_ROW'][i]) - 1
            Q_flag =  m11_QP['QUALITYFLAG'][pix_row][pix_col]
            if Q_flag == 2 :
                drop_list.append(i)
        m11_df = m11_df.drop(m11_df.index[drop_list])        
        m11_df = m11_df.reset_index(drop=True)
        
        return m8_df, m11_df
    
    elif not os.path.exists(m8):
        print "ERROR when trying to open file, this file does not exist:" + m8    
        m8_df, m11_df = [pd.DataFrame(),pd.DataFrame()]
        return m8_df, m11_df
    elif not os.path.exists(m11):
        print "ERROR when trying to open file, this file does not exist:" + m11    
        m8_df, m11_df = [pd.DataFrame(),pd.DataFrame()]
        return m8_df, m11_df

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
            
def plot_pixel_size(plot_list):
    
    data_name= ['Meteosat-8','Metiosat-11']
    max_frp=np.max([np.max(m8_df.FRP),np.max(m11_df.FRP)])
    graph_max=roundup(max_frp)
    print graph_max

    fig = plt.figure()    
    m8_latitude = m8_df['LATITUDE'].values; m11_latitude = m11_df['LATITUDE'].values
    m8_longitude = m8_df['LONGITUDE'].values; m11_longitude = m11_df['LONGITUDE'].values
    m8_p_angle = m8_df['PIXEL_VZA'].values  ; m11_p_angle = m11_df['PIXEL_VZA'].values 
    m8_p_size = m8_df['PIXEL_SIZE'].values; m11_p_size = m11_df['PIXEL_SIZE'].values

    ax = fig.add_subplot(121)
    m8_z = np.polyfit(m8_p_angle,m8_p_size,3)
    m8_f = np.poly1d(m8_z)   
    m8_x_new = np.linspace(m8_p_angle[0], m8_p_size[-1], 50)
    m8_y_new = m8_f(m8_x_new)
    ax.plot(m8_p_angle,m8_p_size,'o', m8_x_new, m8_y_new)
    #ax.scatter(m8_p_angle,m8_p_size,s = 50, alpha = 0.8, cmap='hot',marker ='.')
    ax = fig.add_subplot(122)
    #ax.scatter(m11_p_angle,m11_p_size,s = 50, alpha = 0.8, cmap='hot',marker ='.') 
    
    m11_z = np.polyfit(m11_p_angle,m11_p_size,3)
    m11_f = np.poly1d(m11_z)   
    m11_x_new = np.linspace(m11_p_angle[0], m11_p_size[-1], 50)
    m11_y_new = m11_f(m11_x_new)
    print m8_z, m11_z
    ax.plot(m11_p_angle,m11_p_size,'o', m11_x_new, m11_y_new)



    plt.savefig('pixel_frp_plot' + date_time, dpi = 500)
    plt.close()
    
#plot_pixel_size(plot_list)
    

def plot_single_map(latitude,longitude,title):
    #formated_dt = date_time[0:4]+'-'+date_time[4:6]+'-'+date_time[6:8]+'_'+date_time[8:10]+':'+ date_time[10:12] +':00'    
    fig = plt.figure()
    #FRP = df['FRP'].values
      
    ax = fig.add_subplot(111)
    #m_m8 = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,
    #        llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    #(llcrnrlon=-35,llcrnrlat=-40,urcrnrlon=90,urcrnrlat=40,
    #                               resolution='i',projection='tmerc',lon_0=15,lat_0=0,ax=ax)
    m_m8 = Basemap(llcrnrlon=-35,llcrnrlat=-40,urcrnrlon=60,urcrnrlat=40,
                                   resolution='i',projection='tmerc',lon_0=15,lat_0=0,ax=ax) 
    m_m8.drawcoastlines(linewidth = 0.5)
    m_m8.drawmapboundary()
    #m.fillcontinents(lake_color='aqua',zorder=0)
    m_m8.drawparallels(np.arange(-30,60.,10), labels = [True],fontsize=5 )
    m_m8.drawmeridians(np.arange(-40.,100.,10))
    #m_m8.drawcountries(linewidth=0.25)
    #latitude = df['mean_lat'].values
    #longitude = df['mean_long'].values
    x, y = m_m8(longitude, latitude)
    #sc = plt.scatter(x, y, c='b', s = 8, alpha = 0.8, cmap='hot',marker =',',linewidth=0.0, label = 'VZA Binned FRP')
                            #zorder=1,norm=mpl.colors.SymLogNorm(linthresh=2, vmin=0, vmax=1000))
                            
    #latitude = df1['mean_lat'].values
    #longitude = df1['mean_long'].values    
    #sc = plt.scatter(x, y, c='g', s = 8, alpha = 0.8, cmap='hot',marker =',',linewidth=0.0, label = 'VZA Binned FRP')

    #latitude = df2['mean_lat'].values
    #longitude = df2['mean_long'].values    
    #sc = plt.scatter(x, y, c='k', s = 8, alpha = 0.8, cmap='hot',marker =',',linewidth=0.0, label = 'VZA Binned FRP')

    #latitude = df3['mean_lat'].values
    #longitude = df3['mean_long'].values    
    sc = plt.scatter(x, y, c='orange', s = 8, alpha = 0.8, cmap='hot',marker =',',linewidth=0.0, label = 'VZA Binned FRP')

    x,y = m_m8(-3.4,0) 
    m_m8.plot(x,y,"bo", ms = 3, label = 'Met-11 Sub-satellite loc')   ; m_m8.plot(x,y,"o", ms = 10, markerfacecolor = "none", markeredgecolor= "b")#,alpha = 0.5) #; m_m8.plot(x,y,"o", ms = 100,markerfacecolor = "none", markeredgecolor= "b") ; m_m8.plot(x,y,"o", ms = 150,markerfacecolor = "none", markeredgecolor= "b")  ;  m_m8.plot(x,y,"o", ms = 200,markerfacecolor = "none", markeredgecolor= "b")     
    x,y = m_m8(41.5,0) 
    m_m8.plot(x,y,"go", ms = 3, label = 'Met-8 Sub-satellite loc')   ;  m_m8.plot(x,y,"o", ms = 10, markerfacecolor = "none", markeredgecolor= "g")#, alpha = 0.5) #; m_m8.plot(x,y,"o", ms = 100,markerfacecolor = "none", markeredgecolor= "g") ; m_m8.plot(x,y,"o", ms = 150,markerfacecolor = "none", markeredgecolor= "g")  ;  m_m8.plot(x,y,"o", ms = 200,markerfacecolor = "none", markeredgecolor= "g")  
    plt.title(title,fontsize=9)
    
    plt.legend(loc='lower left', fontsize=7)
    #cb = m_m8.colorbar(sc, pad=0.001, ticks = [0, 1000] )#location='bottom')
    #cb.ax.tick_params(labelsize=8)
    #cb.ax.set_ylabel('FRP (MW)',fontsize=9)
    #cb.ax.set_yticklabels([])
    #cb.set_ticks([0,200,400,600,800,1000,1200,1400])         # need editting to be more case adaptive range
    #cb.ax.set_ylabel('FRP',fontsize=9)

    plt.savefig( ('C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/Plots/'+title), dpi = 500)
    plt.close()



def plot_scatter_fit(x,y,colour,z, perf_fit,fit_points, subtitle, title,textstr1,textstr2,c_max,fire_max): # 
    fig = plt.figure(figsize = [5,5])
    ax1 = fig.add_subplot(111)
    fire_max  = int(math.ceil(fire_max / 500.0)) * 500
    
    #print m1,c1,m2,c2
    #ax1.plot([w, z], [w, z],'k-', color = 'r')
    sc = ax1.scatter(x,y, c=colour, alpha = 0.5,vmin=0,vmax=c_max,linewidth=0.0)
    #ax1.scatter(outlier_fires['M11_summed_FRP'],outlier_fires['M11_summed_FRP'], c= 'r', alpha = 0.1, linewidth=0.3 )
    ax1.set_ylabel(' M8 Pixel count ',fontsize=9)
    ax1.set_xlabel('M11 Pixel count',fontsize=9)
    cb = plt.colorbar(sc)
    cb.ax.set_ylabel('Absolute VZA difference (Degrees)',fontsize=9)
    plt.subplots_adjust(top = 0.75)
    ax1.plot(fit_points,z,color='r' , label = 'OLS fit')
    ax1.plot(fit_points,perf_fit,color = 'b', label = '1:1 fit', ls = 'dashdot')
    plt.gca().set_aspect('equal', adjustable='box')
    axes = plt.gca()
    axes.set_xlim([0,fire_max])
    axes.set_ylim([0,fire_max])
    ax1.legend(loc='upper left', fontsize=8)
    plt.xticks(np.arange(0, fire_max+1, step=500),fontsize=7)
    plt.yticks(np.arange(0, fire_max+1, step=500),fontsize=7)
    #plt.xticks(np.arange(0, 20, step=2),fontsize=7)
    #plt.yticks(np.arange(0, 20, step=2),fontsize=7)
    plt.setp(ax1.get_xticklabels(), rotation=40, horizontalalignment='right')
    ax1.text(0.95,0.05,textstr1, transform=ax1.transAxes, fontsize=8,
        verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))# bbox=props)
    #ax1.text(0.5,0.95,textstr2, transform=ax1.transAxes, fontsize=8)
        #bbox=dict(facecolor='white', alpha = 0.8,  boxstyle='round,pad=1'))# bbox=props)
    plt.suptitle(textstr2,fontsize=10)
    plt.tight_layout()
    plt.savefig(title + subtitle, dpi = 800, bbox_inches='tight')
    plt.close()
    #ax1.scatter(m11_frp,m11_vza,c='b',)


def plot_ODR_scatter_residuals(m11_frp,m8_frp,m11_err,m8_err,fit_y, perf_fit,fit_points, residual,adjusted_err,params,vza_max,fire_max,fire_count,subtit,title):
    fire_max = int(math.ceil(fire_max / 500.0)) * 500
    fig = plt.figure(figsize = [5,5])
    ax1 = fig.add_subplot(111)
    plt.ylabel("M8 FRP [MW]")
    plt.xlabel("M11 FRP [MW]")
    plt.title("ODR Fit to Data")
    ax1.plot(m11_frp,m8_frp,'.',color = 'k',ms =4, linewidth=0.0)
    ax1.errorbar(m11_frp, m8_frp, xerr=m11_err, yerr=m8_err, fmt='none',ecolor='darkorange', alpha = 0.7, label = 'FRP uncertainty')
    ax1.plot(fit_points,perf_fit,color = 'b', label = '1:1 fit', ls = 'dashdot')
    ax1.plot(fit_points,fit_y,color = 'r', label = 'ODR fit')
    axes = plt.gca()
    axes.set_xlim([0,fire_max])
    axes.set_ylim([0,fire_max])
    plt.xticks(np.arange(0, fire_max+1, step=500),fontsize=7)
    plt.yticks(np.arange(0, fire_max+1, step=500),fontsize=7)
    plt.setp(ax1.get_xticklabels(), rotation=30, horizontalalignment='right')
    textstr1 ='Fire count=%d \n slope=%.3f \n intercept=%.3f'%(fire_count,params[1],params[0])
    ax1.text(0.95,0.95,textstr1, transform=ax1.transAxes, fontsize=8,va='top',ha='right', bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))
    ax1.legend(loc='upper left', fontsize=8)
    ax1.grid()    
    
    #residuals = fig.add_subplot(122) # 3 rows, 1 column, subplot 2
    #residuals.plot(m11_frp,residual,".", color = 'k',  ms =4, alpha =0.5)
    #residuals.errorbar(x=m11_frp,y=residual,yerr=adjusted_err, fmt="none", ecolor='darkorange',  label = "Residuals")
    #residuals.set_xlim(ax1.get_xlim())
    #plt.axhline(y=0, color='b')
    #plt.xlabel("M11 FRP [MW]")
    #plt.title("M8 Residuals with Adjusted Error")
    
    #plt.xticks(np.arange(0, fire_max+1, step=250),fontsize=7)
    #plt.yticks(np.arange(0, fire_max+1, step=250),fontsize=7)
    #plt.setp(ax1.get_xticklabels(), rotation=30, horizontalalignment='right')
    #plt.ticklabel_format(style='plain', useOffset=False, axis='x')
    #plt.ylabel("M8 FRP Residuals [MW]")
    #residuals.grid()
    
    plt.tight_layout()
    plt.savefig(title +subtit, dpi = 800,bbox_inches='tight')
    plt.close()

def plot_scatter_residuals(m11_frp,m8_frp,fit_y,residual,slope,intercept,vza_max,fire_count,subtit,title):
    fig = plt.figure(figsize = [8,4])
    ax1 = fig.add_subplot(121)
    plt.ylabel("M8 FRP [MW]")
    plt.xlabel("M11 FRP [MW]")
    plt.title("OLS Fit to Data")
    ax1.plot(m11_frp,m8_frp,'.',color = 'k',ms =4, linewidth=0.0)
    #ax1.errorbar(m11_frp, m8_frp, xerr=m11_err, yerr=m8_err, fmt='none',ecolor='darkorange', alpha = 0.7, label = 'FRP uncertainty')
    ax1.plot(m11_frp,fit_y,color = 'r', label = 'OLS fit')
    textstr1 ='Fire count=%d \n slope=%.3f \n intercept=%.3f'%(fire_count,slope,intercept)
    ax1.text(0.95,0.95,textstr1, transform=ax1.transAxes, fontsize=8,va='top',ha='right', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))
    ax1.legend(loc='upper left', fontsize=8)
    ax1.grid()    
    
    residuals = fig.add_subplot(122) # 3 rows, 1 column, subplot 2
    residuals.plot(m11_frp,residual,".", color = 'k',  ms =4, alpha =0.5)
    residuals.set_xlim(ax1.get_xlim())
    plt.axhline(y=0, color='b')
    plt.xlabel("M11 FRP [MW]")
    plt.title("M8 Residuals")
    plt.ticklabel_format(style='plain', useOffset=False, axis='x')
    plt.ylabel("M8 FRP Residuals [MW]")
    residuals.grid()
    
    plt.tight_layout()
    plt.savefig(title +subtit, dpi = 800,bbox_inches='tight')
    plt.close()


    
def plot_stats(x,y,z,w,title1, title2, labels): # 
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    label_1, label_2, label_3 = labels
    #print m1,c1,m2,c2
    #ax1.plot([w, z], [w, z],'k-', color = 'r')
    #sc = ax1.scatter(x,y, c=colour, alpha = 0.5,vmin=0,vmax=c_max)
    #ax1.set_ylabel('FRP (MW) ',fontsize=9)
    ax1.set_xlabel('Absolute VZA difference Bin',fontsize=9)
    #plt.yticks(np.arange(0, 1, step=0.1))
    #cb = plt.colorbar(sc)
    #cb.ax.set_ylabel('Absolute VZA difference (Degrees)',fontsize=9)
    #plt.subplots_adjust(top = 0.85)
    ax1.plot(x,y,color='dodgerblue', label = label_1)
    ax1.plot(x,z,color='darkorange', label = label_2)
    ax1.plot(x,w,color='forestgreen', label = label_3)
    #ax1.plot(x,w,color='palevioletred', label = 'Standard Error')
    ax1.legend()
    #ax1.text(0.3,0.95,textstr1, transform=ax1.transAxes, fontsize=8,
    #   va='top',ha='right')# bbox=props)
    plt.suptitle('Statistics on linear fit binned by their ' + title2 ,fontsize=10)
    #ax1 = fig.add_subplot(122)
    #ax1.plot(np.unique(w), np.poly1d(np.polyfit(w,z, 1))(np.unique(w)))
    #ax1.scatter(w,z, c='r')
    #ax1.set_ylabel('VZA M11',fontsize=9)
    #ax1.set_xlabel('VZA M8',fontsize=9)
    plt.savefig(title1+ title2 , dpi = 700)
    plt.close()
    #ax1.scatter(m11_frp,m11_vza,c='b',)

    'Slope','R$^2$','Standard Error'
    
    #ax1 = fire_stats.plot.scatter(x=m8_frp,y=m11_frp)
    


def plot_compare(plot_list,date_time):
    formated_dt = date_time[0:4]+'-'+date_time[4:6]+'-'+date_time[6:8]+'_'+date_time[8:10]+':'+ date_time[10:12] +':00'    
    m8_df,m11_df = plot_list
    fig = plt.figure()
    latitude = m8_df['LATITUDE'].values
    longitude = m8_df['LONGITUDE'].values
    FRP = m8_df['FRP'].values
      
    ax = fig.add_subplot(121)
    m_m8 = Basemap(llcrnrlon=-30,llcrnrlat=-35,urcrnrlon=60,urcrnrlat=40,
                                    resolution='i',projection='tmerc',lon_0=15,lat_0=0,ax=ax)
    #m_m8 = Basemap(llcrnrlon=28,llcrnrlat=8,urcrnrlon=32,urcrnrlat=10,
    #                                resolution='i',projection='tmerc',lon_0=30,lat_0=9,ax=ax)   
    m_m8.drawcoastlines(linewidth = 0.5)
    m_m8.drawmapboundary()
    #m.fillcontinents(lake_color='aqua',zorder=0)
    m_m8.drawparallels(np.arange(-35,60.,20), labels = [True],fontsize=5 )
    m_m8.drawmeridians(np.arange(-40.,100.,20))
    #m_m8.drawcountries(linewidth=0.25)
    x, y = m_m8(longitude, latitude)
    sc = plt.scatter(x, y, c=FRP, s = 10, alpha = 0.8, cmap='hot',marker =',',
                            zorder=1,norm=mpl.colors.SymLogNorm(linthresh=2, vmin=0, vmax=1000))
    plt.title('Meteosat-8',fontsize=9)
    cb = m_m8.colorbar(sc, pad=0.001, ticks = [0, 1000] )#location='bottom')
    cb.ax.tick_params(labelsize=8)
    cb.ax.set_yticklabels([])
    #cb.set_ticks([0,200,400,600,800])         # need editting to be more case adaptive range
    #cb.ax.set_ylabel('FRP',fontsize=9)


    latitude = m11_df['LATITUDE'].values
    longitude = m11_df['LONGITUDE'].values
    FRP = m11_df['FRP'].values
    
    ax = fig.add_subplot(122)
    m_m11 = Basemap(llcrnrlon=-30,llcrnrlat=-35,urcrnrlon=60,urcrnrlat=40,
                                    resolution='i',projection='tmerc',lon_0=15,lat_0=0,ax=ax)
    #m_m11 = Basemap(llcrnrlon=28,llcrnrlat=8,urcrnrlon=32,urcrnrlat=10,
    #                                resolution='i',projection='tmerc',lon_0=30,lat_0=9,ax=ax)    
    m_m11.drawcoastlines(linewidth = 0.5)
    m_m11.drawmapboundary()
    #m.fillcontinents(lake_color='aqua',zorder=0)
    m_m11.drawparallels(np.arange(-35,60.,20), labels = [True],fontsize=5 )
    m_m11.drawmeridians(np.arange(-40.,100.,20),labels = [True])
    #m_m11.drawcountries(linewidth=0.25)
    # m_ice.shadedrelief()
    x, y = m_m11(longitude, latitude)
    sc = plt.scatter(x, y, c=FRP, s = 10, alpha = 0.8, cmap='hot',marker =',',
                            zorder=1,norm=mpl.colors.SymLogNorm(linthresh=2, vmin=0, vmax=1000))
    plt.title('Meteosat-11',fontsize=9)
    plt.suptitle('FRP-PIXEL on ' + formated_dt ,fontsize=12 , y=0.09)    
    #fig.subplots_adjust(right=0.8)
    rcParams['figure.subplot.wspace'] = 0.2 # more width between subplots
    #cb_ax = fig.add_axes([0, 800, 0,800])
    cb = m_m11.colorbar(sc,  pad=0.001)#location='bottom')
    cb.ax.tick_params(labelsize=8)
    #cb.ax.set_yticklabels([0,10,50,100,200,500,1000])         # need editting to be more case adaptive range
    cb.ax.set_ylabel('FRP (MW)',fontsize=9)
    #fig.tight_layout()
    plt.savefig('FRP_cluster_' + date_time, dpi = 500)
    plt.close()


