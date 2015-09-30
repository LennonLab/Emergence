from __future__ import division
import os
import sys

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/simplex/tools/metrics")
import metrics


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


def test_Rlogskew():

    """ test for log-skew of a frequency distribution """

    assert metrics.Rlogskew([1,1]) == 'S < 2, cannot compute log-skew'
    assert metrics.Rlogskew([1.0,123]) == 'S < 2, cannot compute log-skew'
    assert metrics.Rlogskew([5,5,5,5,5]) == '0 variance, cannot compute log-skew'
    assert metrics.Rlogskew([5,5,5,5,5,5,5,5,1]) == -3.0



def test_Preston():

    """ test for the value of alpha for Preston's lognormal distribution """
    assert metrics.Preston([0]) == 'sum <= 0, cannot compute'
    assert metrics.Preston([1]) == (0.72823274985908792, 3.0525896286614604)
    assert metrics.Preston([1,2,3,4,5,6,7,8,9,10]) == (0.16113216807474443, 1123.4019840427286)


def test_Berger_Parker():

    """ test for the value of the Berger-Parker index of dominance """
    assert metrics.Berger_Parker([100]) == 1.0
    assert metrics.Berger_Parker([0]) == 'sum <= 0, cannot compute'
    assert metrics.Berger_Parker([2, 6, 76, 1, -1]) == 'all elements must be < 0'
