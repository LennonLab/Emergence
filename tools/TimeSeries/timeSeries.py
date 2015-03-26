from __future__ import division
import  matplotlib.pyplot as plt
#from pandas.tools.plotting import autocorrelation_plot
import numpy as np
from numpy import cumsum, log, polyfit, sqrt, std, subtract, random, square, pi, exp
from numpy.random import randn 
from numpy import *  
import sys



def get_lag(x):
    
    #plt.close()
    #plt.figure()
    #autocorrelation_plot(x)
    
    t = 0
    result = []
    while t < 100:
        corcof = np.corrcoef(np.array([x[0:len(x)-t], x[t:len(x)]]))
        result.append(corcof[1][0])
        t+=1
        
    #plt.plot(result, color = 'm')
    #plt.xlim(0, len(result))
    #plt.ylim(-1.0,1.0)
    #plt.show()
    
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    lagVal = 0
    autocorrVal = 0
    
    for i, autocorrVal in enumerate(result):
        if max(result[i:]) < 0.2:
            return i
    
    return 0



def get_uncorrelated(_list):
    
    if len(_list) < 100:
        return 0
    
    elif np.var(_list) == 0.0:
        return 1
    
    else:    
        lag = get_lag(_list)
    
    if lag == 0:
        return 0        
    else:                    
        return lag
      
   
 
def normcdf(X):
    (a1,a2,a3,a4,a5) = (0.31938153, -0.356563782, 1.781477937, -1.821255978, 1.330274429)
    L = abs(X)
    K = 1.0 / (1.0 + 0.2316419 * L)
    w = 1.0 - 1.0 / sqrt(2*pi)*exp(-L*L/2.) * (a1*K + a2*K*K + a3*pow(K,3) + a4*pow(K,4) + a5*pow(K,5))
    if X<0:
        w = 1.0-w
    return w
 
 
def vratio(a, lag = 2, cor = 'hom'):
    t = (std((a[lag:]) - (a[1:-lag+1])))**2
    b = (std((a[2:]) - (a[1:-1]) ))**2
 
    n = float(len(a))
    mu  = sum(a[1:n]-a[:-1])/n
    m=(n-lag+1)*(1-lag/n)
    b=sum(square(a[1:n]-a[:n-1]-mu))/(n-1)
    t=sum(square(a[lag:n]-a[:n-lag]-lag*mu))/m
    vratio = t/(lag*b)
 
    la = float(lag)
     
 
    if cor == 'hom':
        varvrt=2*(2*la-1)*(la-1)/(3*la*n)
 
 
    elif cor == 'het':
          varvrt=0;
          sum2=sum(square(a[1:n]-a[:n-1]-mu)); 
          for j in range(lag-1):
             sum1a=square(a[j+1:n]-a[j:n-1]-mu); 
             sum1b=square(a[1:n-j]-a[0:n-j-1]-mu)
             sum1=dot(sum1a,sum1b); 
             delta=sum1/(sum2**2);
             varvrt=varvrt+((2*(la-j)/la)**2)*delta
 
    zscore = (vratio - 1) / sqrt(float(varvrt))
    pval = normcdf(zscore);
 
    return  vratio, zscore, pval
 


def hurst(ts):
    """Returns the Hurst Exponent of the time series vector ts"""
    # Create the range of lag values
    lags = range(2, 100)

    # Calculate the array of the variances of the lagged differences
    tau = [sqrt(std(subtract(ts[lag:], ts[:-lag]))) for lag in lags]
    
    for index, val in enumerate(tau): # remove 0's
        if val <= 0.0:
            tau.pop(index)
            lags.pop(index)
    
    # Use a linear fit to estimate the Hurst Exponent
    lags = [0.00000001 if e == 0 else e for e in lags]
    tau = [0.00000001 if e == 0 else e for e in tau]

    poly = polyfit(log(lags), log(tau), 1)
    
    # Return the Hurst exponent from the polyfit output
    return poly[0]*2.0

        

def get_burnin(Nlist):
    
    if len(Nlist) < 100:
        return 0
    
    Nlist = Nlist[-100:]
        
    Hurst = hurst(Nlist)
    Vratio, zscore, pval = vratio(np.array(Nlist))
    
    return [Hurst, pval]
    


"""
# Create a Geometric Brownian Motion, Mean-Reverting and Trending Series
gbm = log(cumsum(randn(100000))+1000)
print 'gbm:', gbm
mr = log(randn(100000)+1000)
tr = log(cumsum(randn(100000)+1)+1000)

# Output the Hurst Exponent for each of the above series
# and the price of Google (the Adjusted Close price) for 
# the ADF test given above in the article
print "Hurst(GBM):   %s" % hurst(gbm),
gbm_vr = vratio(gbm)
print gbm_vr[2]

print "Hurst(MR):    %s" % hurst(mr),
mr_vr = vratio(mr)
print mr_vr[2]

print "Hurst(TR):    %s" % hurst(tr),
tr_vr = vratio(tr)
print tr_vr[2]
"""
