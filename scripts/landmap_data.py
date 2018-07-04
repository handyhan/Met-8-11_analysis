import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import numpy.ma as ma
import pandas as pd
import h5py
import sys
import math
from MSG_read import extract_data,plot_compare,plot_scatter_fit, plot_stats,plot_single_map,plot_ODR_scatter_residuals,plot_scatter_residuals
from math import sin, cos, sqrt, atan2, radians
from pandas.plotting import table
import statsmodels.api as smapi
from statsmodels.formula.api import ols
import statsmodels.graphics as smgraphics
import scipy.odr as scodr 
import scipy.special, scipy.stats

pd.options.mode.chained_assignment = None  
pd.set_option('display.max_columns', 500)


hdf = h5py.File('landmap.h5','r+')
columns = ['X', 'Y', '11','14', '20', '30', '40', '50', '60', '70', '90', '100', '110', '120', '130', '140', '150', '160','170', '180', '190', '200', '210', '220']
df = pd.DataFrame()
for item in columns:
    df[item] =  hdf[item]['values'][:]
    
    
#df['11'] = df[df['11'] > 0 ]
df = df[df['11','14', '20', '30', '40', '50', '60', '70', '90', '100', '110', '120', '130', '140', '150', '160','170', '180', '190', '200', '210', '220'] !=
print df

fig = plt.figure()
ax1 = fig.add_subplot(111)
sc = ax1.scatter(df['X'],df['Y'], c=df['11'])
cb = plt.colorbar(sc)
plt.tight_layout()
plt.savefig('landcover_map', dpi = 800)
plt.close()
