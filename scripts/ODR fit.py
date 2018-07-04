# Fitting with uncertainties in x and y 
#http://www.physics.utoronto.ca/~phy326/python/odr_fit_to_data.py
#Simulated experiment:
#Direct current through a diode as a function of the applied voltage  
import os
import inspect               # http://docs.python.org/library/inspect.html
import matplotlib            # http://matplotlib.sourceforge.net
import numpy  as np               # http://numpy.scipy.org/
import scipy                 # http://scipy.org/
import scipy.odr as scodr 
import scipy.special, scipy.stats
import matplotlib.pyplot as plt
from cluster_analysis import read_in_from_csv
from MSG_read import extract_data,plot_compare,plot_scatter_fit, plot_stats,plot_single_map




dates = ["20180515","20180516" ,"20180517"] #,"20180518","20180519","20180520",'20180521','20180522','20180523',"20180524",'20180525','20180526',"20180527"]
fire_data, fires_m8, fires_m11, pixels_m8, pixels_m11 = read_in_from_csv(dates)


# You can ignore x uncertainties if desired
x_uncertainties_ignored = False

## Choose function to fit:

def linear(p,x) :
    # A linear function with:
    #   Constant Background          : p[0]
    #   Slope                        : p[1]
    return p[0]+p[1]*x
func=linear
p_guess = (2,2)    # can choose first guess intercept to be zero by changing linear (replate p[0] with 0, then only fitting slope)

data_file = "linear_xy_errors.txt" # Data with large x uncertainties

# import in x,y and uncertainties
#x, dx, y_data, dy = numpy.loadtxt(data_file, comments='#', unpack = True)  ##!!!!
if x_uncertainties_ignored:
    # x uncertainties cannot be set exactly to zero with crashing the fit,
    #   but a tiny value seems to do the job.
    dx = len(x)*[1e-99]

# Load data for ODR fit
x = fire_data['M11_summed_FRP']
y_data = fire_data['M8_summed_FRP']
dx = fire_data['M11_FRP_uncertainty']
dy = fire_data['M8_FRP_uncertainty']
vza_diff = fire_data['VZA_diff']
colour = vza_diff
c_max = vza_diff.max()
print dx, dy

data = scipy.odr.RealData(x=x, y=y_data, sx=dx, sy=dy)             ## put data in scipy object
# Load model for ODR fit
model = scipy.odr.Model(func)                                     ## make an object which is the model

## Now fit model to data
#	job=10 selects central finite differences instead of forward differences
#		when doing numerical differentiation fo function
#	maxit is maximum number of iterations
fit = scipy.odr.ODR(data, model, p_guess, maxit=5000,job=10)      
output = fit.run()
p = output.beta 	# 'beta' is an array of the parameter estimates
cov = output.cov_beta   # parameter covariance matrix
uncertainty = output.sd_beta # parameter standard uncertainties

# Calculate initial residuals and the 'adjusted error' for each data point                    ###!!!!! look through this for res
delta = output.delta  # estimated x-component of the residuals
eps   = output.eps    # estimated y-component of the residuals
# (xstar,ystar) is the point where the 'residual line' (in black)
#   intersects the 'ellipse' created by xerr & yerr.
xstar = dx*numpy.sqrt( ((dy*delta)**2) / ( (dy*delta)**2 + (dx*eps)**2 ) )
ystar = dy*numpy.sqrt( ((dx*eps)**2) / ( (dy*delta)**2 + (dx*eps)**2 ) )
adjusted_err = numpy.sqrt(xstar**2 + ystar**2)
# residual is positive if the point lies above the fitted curve,
#             negative if below
residual = numpy.sign(y_data-func(p,x))*numpy.sqrt(delta**2 + eps**2)
# number of degrees of freedom for fit
dof = len(x) - len(p_guess)
#"Quasi-chi-squared" is defined to be the [total weighted sum of squares] / dof
#	i.e. same as numpy.sum((residual/adjusted_err)**2)/dof or
#       numpy.sum(((output.xplus-x)/dx)**2+((y_data-output.y)/dy)**2)/dof
#	This converges to the conventional chi-squared for zero x uncertainties.
quasi_chisq = output.res_var

print "\nQuasi Chi-Squared/dof   = {0:10.5f}, Chi-Squared CDF = {1:10.5f}%".\
    format(quasi_chisq, 100.*float(scipy.special.chdtrc(dof,dof*quasi_chisq)))
print "   WARNING:Above CDF is not valid for large x uncertainties!"

print "\nTo use Monte Carlo simulation to more accurately estimate CDF for"
print '      large x uncertainties, re-run program with '
print '     "Run_Monte_Carlo_CDF_Estimator = True" and'
print '     "Number_of_MC_iterations >= 1000." This may take some time\n'
## Plot

fig = fig = plt.figure()
# 3 rows, 1 column, subplot 1
#   3 rows are declared, but there are only 2 plots; this leaves room for text
#       in the empty 3rd row
ax1 = fig.add_subplot(311)
# remove tick labels from upper plot (for clean look)
ax1.set_xticklabels( () ) 
fire_max = round((max(x.max(),y.max())) / 500.0) * 500.0
plt.ylabel("M8 FRP (MW)")

plt.title("Orthogonal Distance Regression Fit to Data")
# Plot data as red circles, and fitted function as (default) line.
#   For a smooth look,generate many x values for plotting the model 
sc = ax1.scatter(x,y_data, c=colour, alpha = 0.5,vmin=0,vmax=c_max,linewidth=0.0)
ax1.errorbar(x,y,xerr=dx, yerr=dy,fmt="r+")
#fit.plot(x,y,'ro', x_model, func(p,x_model))
# Add error bars on data as red crosses.
ax1.set_yscale('linear')
#   draw starting guess as dashed green line ('r-')

# separate plot to show residuals
residuals = fig.add_subplot(312) # 3 rows, 1 column, subplot 2
residuals.errorbar(x=x,y=residual, yerr=adjusted_err,
                   			fmt="r+", label = "Residuals")
# make sure residual plot has same x axis as fit plot
residuals.set_xlim(ax1.get_xlim())
# Draw a horizontal line at zero on residuals plot
plt.axhline(y=0, color='b')
# Label axes
plt.xlabel("FRP [MW]")
# These data look better if 'plain', not scientific, notation is used,
#   and if the tick labels are not offset by a constant (as is done by default).
plt.ticklabel_format(style='plain', useOffset=False, axis='x')
plt.ylabel("Residuals")

residuals.grid()

plt.savefig('test_data', dpi = 800)
plt.close()




############################## PRINT THE RESULTS ##############################
print "***********************************************************"
print "               ORTHOGONAL DISTANCE REGRESSION"
print "***********************************************************\n"
print "ODR algorithm stop reason: " + output.stopreason[0]
print "\nFit {0} Data points from file: {1}".format(len(x),data_file)
print "To Model :"
print " ",inspect.getsource(func)
if x_uncertainties_ignored:
    print "** WARNING: x uncertainties set to zero in fit. **\n"

print "Estimated parameters and uncertainties"
for i in range(len(p)) :
    print ("   p[{0}] = {1:10.5g} +/- {2:10.5g}"+
           "          (Starting guess: {3:10.5g})").\
            format(i,p[i],uncertainty[i],p_guess[i])

print "\nCorrelation Matrix :"
for i,row in enumerate(cov):
    for j in range(len(p)) :
        print "{0:< 8.3g}".format(cov[i,j]/numpy.sqrt(cov[i,i]*cov[j,j])),
            # Newbie Notes: "{0:< 8.3g}" left justifies output with space in
            #   front of positive numbers, with 3 sig figs;
            #   the comma at end of print statement suppresses new line.
    print 
    
##### END of odr_fit_to_data.py