# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 10:41:46 2018

@author: Hannah.N
"""
import statsmodels.api as smapi
from statsmodels.formula.api import ols
import statsmodels.graphics as smgraphics
import scipy.odr as scodr 
import scipy.special, scipy.stats

import numpy as np
import numpy.ma as ma
import pandas as pd
import sys
import math


def regression_log(df,x_name,y_name):
    df[y_name] = np.log10(df[y_name])
    df[x_name] = np.log10(df[x_name])
    mod = ols(formula = y_name +' ~ '+ x_name , data = df)
    res= mod.fit()
    results = (res.summary())    
    test = res.outlier_test()
    outlier_ix = test[test['bonf(p)'] < 0.05].index
    slope = res.params[x_name]
    intercept = res.params.Intercept
    r_squared = res.rsquared
    se = res.bse[x_name]
    residual = res.resid
    
    return results, outlier_ix,slope, intercept,r_squared, se, residual


def regression_free(df,x_name,y_name):
    mod = ols(formula = y_name +' ~ '+ x_name , data = df)
    res= mod.fit()
    results = (res.summary())    
    test = res.outlier_test()
    outlier_ix = test[test['bonf(p)'] < 0.05].index
    slope = res.params[x_name]
    intercept = res.params.Intercept
    r_squared = res.rsquared
    se = res.bse[x_name]
    residual = res.resid
    
    return results, outlier_ix,slope, intercept,r_squared, se, residual

def regression_c0(df,x_name,y_name):
    mod = ols(formula = y_name +' ~ '+ x_name + '+0' , data = df)
    res= mod.fit()
    results = (res.summary())    
    test = res.outlier_test()
    outlier_ix = test[test['bonf(p)'] < 0.05].index
    slope = res.params[x_name]
    r_squared = res.rsquared
    se = se = res.bse[x_name]
    intercept = 0.00
    residual = res.resid
    
    return results, outlier_ix, slope , intercept, r_squared, se, residual



def odr_regression(x,y,x_err,y_err):
        
    def linear(p,x) :        # A linear function with:      #   Constant Background          : p[0]       #   Slope                        : p[1]
        return p[0]+p[1]*x

    func=linear                           # set the function to an object
    p_guess = (0,1)    # can choose first guess intercept to be zero by changing linear (replate p[0] with 0, then only fitting slope)      
    data = scodr.RealData(x=x, y=y, sx=x_err, sy=y_err)     
    model = scipy.odr.Model(func)    
    print data                
    odr = scodr.ODR(data, model, p_guess, maxit=5000,job=10)      
    output = odr.run()
    params = output.beta 	# 'beta' is an array of the parameter estimates
    param_uncertainty = output.sd_beta # parameter standard uncertainties    
    delta = output.delta  # estimated x-component of the residuals
    eps   = output.eps    # estimated y-component of the residuals
    xstar = x_err*np.sqrt( ((y_err*delta)**2) / ( (y_err*delta)**2 + (x_err*eps)**2 ) )
    ystar = y_err*np.sqrt( ((x_err*eps)**2) / ( (y_err*delta)**2 + (x_err*eps)**2 ) )
    adjusted_err = np.sqrt(xstar**2 + ystar**2)
    residual = np.sign(y-func(params,x))*np.sqrt(delta**2 + eps**2)
    return params,param_uncertainty, adjusted_err, residual

