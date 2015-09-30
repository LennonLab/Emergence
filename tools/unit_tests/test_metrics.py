from __future__ import division

import os
import sys
import numpy as np
#import scipy
#from scipy import stats
#from scipy.stats import gaussian_kde
#from scipy.optimize import fsolve
import math


mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/tools/metrics")
import metrics

# content of test1.py

def test_percent_ones():

    """ tests for the percent of species in a list called sad,
    represented by a single individual

    where the function returns 100 * sad.count(1)/len(sad) """

    assert metrics.percent_ones([1,10,1,11,1.1,0.1,1,1,111,1.1]) == 100*4.0/10.0
    assert metrics.percent_ones([1]) == 100
    assert metrics.percent_ones([2, 3, 4, 11, 111, 0.1, 1.1]) == 0.0


def test_percent_pt_one():

    """ test for percent taxa with less than 0.1 percent relative abundance
    in a list called sad """

    assert metrics.percent_pt_one([1]) == 0
    assert metrics.percent_pt_one([1000, 1]) == 50
    assert metrics.percent_pt_one([999, 1]) == 0
