import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.colors as colors
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', 500)


       
    

def plot_single_map(latitude,longitude,colour,title):
    fig = plt.figure()
    #FRP = df['FRP'].values
      
    ax = fig.add_subplot(111)
    ######### other map progections
    #m_m8 = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,
    #        llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    #(llcrnrlon=-35,llcrnrlat=-40,urcrnrlon=90,urcrnrlat=40,
    #                               resolution='i',projection='tmerc',lon_0=15,lat_0=0,ax=ax) #projection='ortho',lat_0=0,lon_0=20 )
    m_m8 = Basemap(projection='ortho',lat_0=0,lon_0=20 )#llcrnrlon=20,llcrnrlat=-35,urcrnrlon=30,urcrnrlat=-25,  
                                   #resolution='i',projection='tmerc',lon_0=25,lat_0=-30,ax=ax) 
    m_m8.drawcoastlines(linewidth = 0.5)
    m_m8.drawmapboundary()
    m_m8.drawparallels(np.arange(-40,60.,10) ) #, labels = [True],fontsize=5
    m_m8.drawmeridians(np.arange(-40.,100.,10))
    #m_m8.drawcountries(linewidth=0.25)
    #x, y = m_m8(longitude, latitude)
    #sc = plt.scatter(x, y, c='b', s = 8, alpha = 0.8, cmap='hot',marker =',',linewidth=0.0, label = 'VZA Binned FRP')
                            #zorder=1,norm=mpl.colors.SymLogNorm(linthresh=2, vmin=0, vmax=1000))
                             
    #x, y = m_m8(longitude2, latitude2)                      
    #sc = plt.scatter(x, y, c=colour, s = 8, alpha = 0.8, cmap='hot',marker =',',linewidth=0.0, label = 'VZA Binned FRP')

    x, y = m_m8(longitude, latitude)
    plt.scatter(x, y, c=colour, s = 5, alpha = 0.5,marker ='.',linewidth=0.0)   # change colour 
    
    
    #x1, y1 = m_m8(longitude2, latitude2)
    #plt.scatter(x1, y1, c='k', s = 8, alpha = 0.2, marker =',',linewidth=0.0)
    #cb = plt.colorbar(sc)
    #cb.ax.set_ylabel('Absolute VZA difference (Degrees)',fontsize=9)
    #x,y = m_m8(latitude,longitude) 
    #m_m8.plot(x,y,"ro", ms = 3, label = 'Met-11 Sub-satellite loc')   ; m_m8.plot(x,y,"o", ms = 10, markerfacecolor = "none", markeredgecolor= "r")#,alpha = 0.5) #; m_m8.plot(x,y,"o", ms = 100,markerfacecolor = "none", markeredgecolor= "b") ; m_m8.plot(x,y,"o", ms = 150,markerfacecolor = "none", markeredgecolor= "b")  ;  m_m8.plot(x,y,"o", ms = 200,markerfacecolor = "none", markeredgecolor= "b")     
    #x,y = m_m8(23.960735,-29.0457) 
    #m_m8.plot(x,y,"bo", ms = 3, label = 'Met-8 Sub-satellite loc')   ;  m_m8.plot(x,y,"o", ms = 10, markerfacecolor = "none", markeredgecolor= "b")#, alpha = 0.5) #; m_m8.plot(x,y,"o", ms = 100,markerfacecolor = "none", markeredgecolor= "g") ; m_m8.plot(x,y,"o", ms = 150,markerfacecolor = "none", markeredgecolor= "g")  ;  m_m8.plot(x,y,"o", ms = 200,markerfacecolor = "none", markeredgecolor= "g")  
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



def plot_scatter_fit(vza_switch,df,x,y,colour,z, perf_fit,fit_points, subtitle, title,textstr1,textstr2,c_max,fire_max): # 
    
    x = df[x]
    y = df[y]
    colour = df[colour]
    if vza_switch == False:

        fig = plt.figure(figsize = [5,5])
        ax1 = fig.add_subplot(111)
        sc = ax1.scatter(x,y, c='orangered',s=15, alpha = 0.5,linewidth=0.0)
        #ax1.scatter(outlier_fires['M11_summed_FRP'],outlier_fires['M11_summed_FRP'], c= 'r', alpha = 0.1, linewidth=0.3 )
        ax1.set_ylabel(y,fontsize=9)
        ax1.set_xlabel(x,fontsize=9)
        plt.subplots_adjust(top = 0.75)
        ax1.plot(fit_points,z,color='b' , label = 'OLS fit')
        ax1.plot(fit_points,perf_fit,color = 'k', label = '1:1 fit', ls = 'dashdot')
        plt.gca().set_aspect('equal', adjustable='box')
        axes = plt.gca()
        axes.set_xlim([0,fire_max])
        axes.set_ylim([0,fire_max])
        ax1.legend(loc='upper left', fontsize=8)
        plt.xticks(np.arange(0, fire_max+1, step=2),fontsize=7)
        plt.yticks(np.arange(0, fire_max+1, step=2),fontsize=7)
        plt.setp(ax1.get_xticklabels(), rotation=40, horizontalalignment='right')
        ax1.text(0.95,0.05,textstr1, transform=ax1.transAxes, fontsize=8,
            verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))# bbox=props)
        #ax1.text(0.5,0.95,textstr2, transform=ax1.transAxes, fontsize=8)
            #bbox=dict(facecolor='white', alpha = 0.8,  boxstyle='round,pad=1'))# bbox=props)
        ax1.text(0.95,0.25,textstr2, transform=ax1.transAxes, fontsize=8,
            verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))
        plt.tight_layout()
        plt.savefig(title + subtitle, dpi = 800, bbox_inches='tight')
        plt.close()
       
        
    elif  vza_switch == True:
            
        fig = plt.figure(figsize = [5,5])
        ax1 = fig.add_subplot(111)
        
        sc = ax1.scatter(x,y, c=colour, alpha = 0.5,vmin=0,vmax=c_max,linewidth=0.0,cmap='bwr')
        #ax1.scatter(outlier_fires['M11_summed_FRP'],outlier_fires['M11_summed_FRP'], c= 'r', alpha = 0.1, linewidth=0.3 )
        ax1.set_ylabel(y,fontsize=9)
        ax1.set_xlabel(x,fontsize=9)
        cb = plt.colorbar(sc)
        cb.ax.set_ylabel(colour,fontsize=9)
        plt.subplots_adjust(top = 0.75)
        ax1.plot(fit_points,z,color='r' , label = 'OLS fit')
        ax1.plot(fit_points,perf_fit,color = 'b', label = '1:1 fit', ls = 'dashdot')
        plt.gca().set_aspect('equal', adjustable='box')
        axes = plt.gca()
        axes.set_xlim([0,fire_max])
        axes.set_ylim([0,fire_max])
        ax1.legend(loc='upper left', fontsize=8)
        plt.xticks(np.arange(0, fire_max+1, step=2),fontsize=7)
        plt.yticks(np.arange(0, fire_max+1, step=2),fontsize=7)
        plt.setp(ax1.get_xticklabels(), rotation=40, horizontalalignment='right')
        ax1.text(0.95,0.05,textstr1, transform=ax1.transAxes, fontsize=8,
            verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))# bbox=props)
        #ax1.text(0.5,0.95,textstr2, transform=ax1.transAxes, fontsize=8)
            #bbox=dict(facecolor='white', alpha = 0.8,  boxstyle='round,pad=1'))# bbox=props)
        plt.suptitle(textstr2,fontsize=10)
        plt.tight_layout()
        plt.savefig(title + subtitle, dpi = 800, bbox_inches='tight')
        plt.close()





def plot_ODR_scatter_residuals(m11_frp,m8_frp,m11_err,m8_err,fit_y, perf_fit,fit_points, residual,adjusted_err,params,vza_max,fire_max,fire_count,subtit,title,textstr1):
    fig = plt.figure(figsize = [5,5])
    ax1 = fig.add_subplot(111)
    plt.ylabel("M8 FRP [MW]")
    plt.xlabel("M11 FRP [MW]")
    plt.title("ODR Fit to Data")
    ax1.plot(m11_frp,m8_frp,'.',color = 'k',ms =4, linewidth=0.0, alpha = 0.4)
    ax1.errorbar(m11_frp, m8_frp, xerr=m11_err, yerr=m8_err, fmt='none',ecolor='darkorange', alpha = 0.7, label = 'FRP uncertainty')
    ax1.plot(fit_points,perf_fit,color = 'b', label = '1:1 fit', ls = 'dashdot')
    ax1.plot(fit_points,fit_y,color = 'r', label = 'ODR fit')
    axes = plt.gca()
    axes.set_xlim([0,fire_max])
    axes.set_ylim([0,fire_max])
    plt.xticks(np.arange(0, fire_max+1, step=500),fontsize=9)
    plt.yticks(np.arange(0, fire_max+1, step=500),fontsize=9)
    plt.setp(ax1.get_xticklabels(), rotation=30, horizontalalignment='right')
    #textstr1 ='Fire count=%d \n slope=%.3f \n intercept=%.3f'%(fire_count,params[1],params[0])
    ax1.text(0.95,0.15,textstr1, transform=ax1.transAxes, fontsize=8,va='top',ha='right', bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))
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

   

def plot_scatter_area(vza_switch,x,y,colour,z, perf_fit,fit_points, subtitle, title,textstr1,textstr2,c_max,fire_max, xaxis_lable,yaxis_lable,scaling_factor): # 
    
    #min_y = math.floor(min(y))
    if vza_switch == False:

        fig = plt.figure(figsize = [5,7])
        ax1 = fig.add_subplot(111)
        sc = ax1.scatter(x,y, c='orangered',s=15, alpha = 0.6,linewidth=0.0)
        #ax1.scatter(outlier_fires['M11_summed_FRP'],outlier_fires['M11_summed_FRP'], c= 'r', alpha = 0.1, linewidth=0.3 )
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
        #ax1.text(0.5,0.95,textstr2, transform=ax1.transAxes, fontsize=8)
            #bbox=dict(facecolor='white', alpha = 0.8,  boxstyle='round,pad=1'))# bbox=props)
        ax1.text(0.95,0.25,textstr2, transform=ax1.transAxes, fontsize=8,
            verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))
        plt.tight_layout()
        plt.savefig(title + subtitle, dpi = 800, bbox_inches='tight')
        plt.close()
       
        
    elif  vza_switch == True:
            
        fig = plt.figure(figsize = [5,5])
        ax1 = fig.add_subplot(111)
        
        sc = ax1.scatter(x,y, c=colour, alpha = 0.5,vmin=0,vmax=c_max,linewidth=0.0)
        #ax1.scatter(outlier_fires['M11_summed_FRP'],outlier_fires['M11_summed_FRP'], c= 'r', alpha = 0.1, linewidth=0.3 )
        ax1.set_ylabel(' M8 '+ axis_lable,fontsize=9)
        ax1.set_xlabel(' M11 '+ axis_lable,fontsize=9)
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
        plt.xticks(np.arange(0, fire_max+1, step=scaling_factor),fontsize=9)
        plt.yticks(np.arange(0, fire_max+1, step=scaling_factor),fontsize=9)
        plt.setp(ax1.get_xticklabels(), rotation=40, horizontalalignment='right')
        ax1.text(0.95,0.05,textstr1, transform=ax1.transAxes, fontsize=8,
            verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))# bbox=props)
        #ax1.text(0.5,0.95,textstr2, transform=ax1.transAxes, fontsize=8)
            #bbox=dict(facecolor='white', alpha = 0.8,  boxstyle='round,pad=1'))# bbox=props)
        plt.suptitle(textstr2,fontsize=10)
        plt.tight_layout()
        plt.savefig(title + subtitle, dpi = 800, bbox_inches='tight')
        plt.close()

def plot_scatter_log(vza_switch,x,y,colour,z, perf_fit,fit_points, subtitle, title,textstr1,textstr2,c_max,fire_max, xaxis_lable,yaxis_lable,scaling_factor): # 
    
    if vza_switch == False:

        fig = plt.figure(figsize = [5,6])
        ax1 = fig.add_subplot(111,aspect='equal')
        sc = ax1.scatter(x,y, c='orangered',s=15, alpha = 0.6,linewidth=0.0)
        #ax1.scatter(outlier_fires['M11_summed_FRP'],outlier_fires['M11_summed_FRP'], c= 'r', alpha = 0.1, linewidth=0.3 )
        ax1.set_ylabel( yaxis_lable ,fontsize=9)
        ax1.set_xlabel(xaxis_lable,fontsize=9)
        plt.subplots_adjust(top = 0.75)
        ax1.plot(fit_points,z,color='b' , label = 'OLS fit')
        ax1.plot(fit_points,perf_fit,color = 'k', label = '1:1 fit', ls = 'dashdot')
        plt.gca().set_aspect('equal', adjustable='box')
        axes = plt.gca()
        #axes.set_xlim([0,fire_max])
        #axes.set_ylim([0,fire_max])
        ax1.legend(loc='upper left', fontsize=8)
        plt.xscale('log')
        plt.yscale('log')
        ax1.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        #ax1.set_xticks([0,1])
        #ax1.set_yticks([0,1])
        plt.setp(ax1.get_xticklabels(), rotation=40, horizontalalignment='right')
        ax1.text(0.90,0.05,textstr1, transform=ax1.transAxes, fontsize=8,
            verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))# bbox=props)
        ax1.text(0.90,0.20,textstr2, transform=ax1.transAxes, fontsize=8,
            verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))
        plt.tight_layout()
        plt.savefig(title + subtitle, dpi = 800, bbox_inches='tight')
        plt.close()
       
        
    elif  vza_switch == True:
            
        fig = plt.figure(figsize = [5,5])
        ax1 = fig.add_subplot(111)
        
        sc = ax1.scatter(x,y, c=colour, alpha = 0.5,vmin=0,vmax=c_max,linewidth=0.0)
        #ax1.scatter(outlier_fires['M11_summed_FRP'],outlier_fires['M11_summed_FRP'], c= 'r', alpha = 0.1, linewidth=0.3 )
        ax1.set_ylabel(' M8 '+ axis_lable,fontsize=9)
        ax1.set_xlabel(' M11 '+ axis_lable,fontsize=9)
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
        plt.xticks(np.arange(0, fire_max+1, step=scaling_factor),fontsize=9)
        plt.yticks(np.arange(0, fire_max+1, step=scaling_factor),fontsize=9)
        plt.setp(ax1.get_xticklabels(), rotation=40, horizontalalignment='right')
        ax1.text(0.95,0.05,textstr1, transform=ax1.transAxes, fontsize=8,
            verticalalignment='bottom',ha='right',bbox=dict(facecolor='white', alpha = 0.8, edgecolor='black', boxstyle='round,pad=1'))# bbox=props)
        #ax1.text(0.5,0.95,textstr2, transform=ax1.transAxes, fontsize=8)
            #bbox=dict(facecolor='white', alpha = 0.8,  boxstyle='round,pad=1'))# bbox=props)
        plt.suptitle(textstr2,fontsize=10)
        plt.tight_layout()
        plt.savefig(title + subtitle, dpi = 800, bbox_inches='tight')
        plt.close()












"""


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

###################################################

   
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
    
###########################################################


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
"""

