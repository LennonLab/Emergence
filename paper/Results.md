#Results
##P1
**S1:** Across microbial datasets, SAD predictions from the maximum entropy theory of ecology (METE) generally failed to explain more than 60% of variation in abundance among species and often produced negative r-squared values, which have no straigtforward interpretation (Figure 1; Table 1).  
**S2:** Previous tests of METE on communities of macroscopic plants and animals produced much greater success, often explaining 90% or greater variation in abundance (Harte et al. 2008, 2009, White et al. 2012, Xiao et al. 2014).  
**S3:** As expected, due to its relatively even form, the broken-stick model (effectively the geometric distribution) performed considerably worse than METE and generally produced negative r-square values.   
**S4:** Negative values were possible because the relationship is not fitted, i.e., estimating variation around a line with a slope of 1.0 and intercept of zero (White et al. 2012, Locey and White 2013, Xiao et al. 2014).  
**S5:** While the log-series (METE) characterizes the form of the SAD better than the broken-stick, microbial SADs are still characterized by disparities in abundance that METE fails to capture.  
**S6:** Both METE and the broken-stick under-predict the abundance of the most abundant species and over-predict the abundance of the rarest species.

##P2
**S1:** We found that the success of METE and the broken-stick were influenced by the two primary state-variables (*N* and *S*) and the primary constraint of average abundance (*N*/*S*) (Table 1). Across each dataset (EMP, HMP, MG-RAST) increasing *N* led to decreasing fits of each model while 




Table 1. 

| Dataset | Model | Variable |  r^2  | p-value |
|:-------:|:-----:|:--------:|:-----:|:-------:|
|   HMP   |   BS  |     N    | -0.39 |   0.0   |
|   HMP   |  METE |     N    | -0.19 |   0.0   |
|   HMP   |   BS  |     S    |       |         |
|   HMP   |  METE |     S    |       |         |
|   HMP   |   BS  |    N/S   |       |         |
|   HMP   |  METE |    N/S   |       |         |
|   EMP   |   BS  |     N    |       |         |
|   EMP   |  METE |     N    |       |         |
|   EMP   |   BS  |     S    |       |         |
|   EMP   |  METE |     S    |       |         |
|   EMP   |   BS  |    N/S   |       |         |
|   EMP   |  METE |    N/S   |       |         |
