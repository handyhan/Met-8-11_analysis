import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.colors as colors
import numpy as np
import pandas as pd



def plot_scatter_simple(x,y,z, perf_fit,fit_points, file_name ,textstr1,textstr2,fire_max,xaxis_lable,yaxis_lable,scaling_factor): # 
    
    fig = plt.figure(figsize = [5,7])
    ax1 = fig.add_subplot(111)
    ax1.scatter(x,y, c='orangered',s=15, alpha = 0.6,linewidth=0.0)
    ax1.set_ylabel( yaxis_lable ,fontsize=9)
    ax1.set_xlabel(xaxis_lable,fontsize=9)
    plt.subplots_adjust(top = 0.75)
    ax1.plot(fit_points,z,color='b' , label = 'OLS fit')
    ax1.plot(fit_points,perf_fit,color = 'k', label = '1:1 fit', ls = 'dashdot')
    plt.gca().set_aspect('equal', adjustable='box')
    axes = plt.gca()
    axes.set_xlim([0,fire_max])
    axes.set_ylim([0,fire_max])
    ax1.legend(loc='upper left', fontsize=8)        
    plt.xticks(np.arange(0, fire_max+1, step=scaling_factor),fontsize=7)
    plt.yticks(np.arange(0, fire_max+1, step=scaling_factor),fontsize=7)
    plt.setp(ax1.get_xticklabels(), rotation=40, horizontalalignment='right')
    ax1.text(0.95,0.05,textstr1, transform=ax1.transAxes, fontsize=8,
        verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))# bbox=props)
    ax1.text(0.95,0.30,textstr2, transform=ax1.transAxes, fontsize=8,
        verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))
    plt.tight_layout()
    plt.savefig(file_name , dpi = 800, bbox_inches='tight')
    plt.close()
       
    
def plot_single_map(latitude,longitude,variable3,title):
    fig = plt.figure()
    #FRP = df['FRP'].values
      
    ax = fig.add_subplot(111)
    ######### other map progections
    #m_m8 = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,
    #        llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    #(llcrnrlon=-35,llcrnrlat=-40,urcrnrlon=90,urcrnrlat=40,
    #                               resolution='i',projection='tmerc',lon_0=15,lat_0=0,ax=ax) #projection='ortho',lat_0=0,lon_0=20 )
    m_m8 = Basemap(llcrnrlon=-30,llcrnrlat=-35,urcrnrlon=60,urcrnrlat=40,   #projection='ortho',lat_0=0,lon_0=20 )#
                                   resolution='i',projection='tmerc',lon_0=10,lat_0=10,ax=ax) 
    #x, y = m_m8(longitude, latitude)
    #sc = plt.scatter(x, y, c='b', s = 8, alpha = 0.8, cmap='hot',marker =',',linewidth=0.0, label = 'VZA Binned FRP')
                            #zorder=1,norm=mpl.colors.SymLogNorm(linthresh=2, vmin=0, vmax=1000))
                             
    #x, y = m_m8(longitude2, latitude2)                      
    #sc = plt.scatter(x, y, c=colour, s = 8, alpha = 0.8, cmap='hot',marker =',',linewidth=0.0, label = 'VZA Binned FRP')

    x, y = m_m8(longitude, latitude)
    
    
    #x1, y1 = m_m8(longitude2, latitude2)
    m_m8.drawcoastlines(linewidth = 0.5)
    m_m8.drawmapboundary(fill_color='lightcyan')
    m_m8.drawparallels(np.arange(-40,60.,10), linewidth=0.5) #, labels = [True],fontsize=5
    m_m8.drawmeridians(np.arange(-40.,100.,10),linewidth=0.5)
    m_m8.drawcountries(linewidth=0.25)
    m_m8.fillcontinents(lake_color='powderblue',color='cornsilk')
    #m_m8.shadedrelief()
    sc = m_m8.scatter(x, y, c='r', s = (variable3/100), alpha = 0.2,marker ='o',linewidth=0.2,zorder=10)   # change colour 
    #x,y = m_m8(0,0) 
    #m_m8.plot(x,y,"go", ms = 3, label = 'Met-11 Sub-satellite loc')   ; m_m8.plot(x,y,"o", ms = 10, markerfacecolor = "none", markeredgecolor= "g")#,alpha = 0.5) #; m_m8.plot(x,y,"o", ms = 100,markerfacecolor = "none", markeredgecolor= "b") ; m_m8.plot(x,y,"o", ms = 150,markerfacecolor = "none", markeredgecolor= "b")  ;  m_m8.plot(x,y,"o", ms = 200,markerfacecolor = "none", markeredgecolor= "b")     
    #x,y = m_m8(0,41.5) 
    #m_m8.plot(x,y,"bo", ms = 3, label = 'Met-8 Sub-satellite loc')   ;  m_m8.plot(x,y,"o", ms = 10, markerfacecolor = "none", markeredgecolor= "b")#, alpha = 0.5) #; m_m8.plot(x,y,"o", ms = 100,markerfacecolor = "none", markeredgecolor= "g") ; m_m8.plot(x,y,"o", ms = 150,markerfacecolor = "none", markeredgecolor= "g")  ;  m_m8.plot(x,y,"o", ms = 200,markerfacecolor = "none", markeredgecolor= "g")  
    plt.title('MODIS fires correctly detected by MSG',fontsize=9)
    
    
    
    #plt.scatter(x1, y1, c='k', s = 8, alpha = 0.2, marker =',',linewidth=0.0)
    #if not isinstance(colour, basestring):
    #    cb = plt.colorbar(sc)
    #    print "Remember to set Coulor bar label and params in def(plot_single_map)"
    #    cb.ax.set_ylabel('Absolute VZA difference (Degrees)',fontsize=9)        
    #    #cb = m_m8.colorbar(sc, pad=0.001, ticks = [0, 1000] )#location='bottom')
    #    cb.ax.tick_params(labelsize=8)
    #    #cb.set_ticks([0,200,400,600,800,1000,1200,1400])         # need editting to be more case adaptive range
    
    
    plt.legend(loc='lower left', fontsize=7)

    plt.savefig( ('C:/Users/Hannah.N/Documents/Earth_Observation/sandbox/Plots/'+title), dpi = 500)
    plt.close()


