# -*- coding: utf-8 -*-

""" The future site for a function that generates the taxa-time relationship


The first step in constructing a species–time relationship is to estimate
the number of species found in a given time span of observation. 

We used a moving window approach where species richness was determined 
for every possible window of each time span (e.g. in a 20 year time series,
 20 one-year windows, 19 two-year windows, etc.). 
 
These values were then averaged within each time span (we used arithmetic means,
 but using geometric means did not affect the results).
 
 Next, we fit the two statistical models (power and logarithmic functions) to
 these averages using two methods: ordinary least squares regression (OLS) on
 log-transformed data and non-linear regression on untransformed data. 
 
 For the non-linear fitting we used the parameter estimates from the OLS as 
 the initial parameter values. OLS and non-linear fitting produced nearly 
 identical parameter estimates (Pearson's r >= 0.99 for both functions). 
 We compared fits of the non-linear power and logarithmic regressions using r2
 values (because there are two parameters in both models there is no need for
 more complex model comparison techniques).
"""

