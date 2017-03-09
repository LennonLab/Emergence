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

    assert metrics.Rlogskew([1,1]) == 'NaN'
    assert metrics.Rlogskew([1.0,123]) == 'NaN'
    assert metrics.Rlogskew([5,5,5,5,5]) == 'NaN'
    assert metrics.Rlogskew([5,5,5,5,5,5,5,5,1]) == -3.0



def test_Preston():

    """ test for the value of alpha for Preston's lognormal distribution """
    assert metrics.Preston([0]) == 'NaN'
    assert metrics.Preston([1]) == (0.72823274985908792, 3.0525896286614604)
    assert metrics.Preston([1,2,3,4,5,6,7,8,9,10]) == (0.16113216807474443, 1123.4019840427286)


def test_Berger_Parker():

    """ test for the value of the Berger-Parker index of dominance, i.e.,
    relative abundance of the most abundant species """

    assert metrics.Berger_Parker([100]) == 1.0
    assert metrics.Berger_Parker([1,1,1,1,1,1,1,1,1,1]) == 0.1
    assert metrics.Berger_Parker([0]) == 'NaN'
    assert metrics.Berger_Parker([2, 6, 76, 1, -1]) == 'NaN'


def test_McNaughton():

    """ test for the value of the McNaughton index of dominance, i.e.,
    perent relative abundance of the two most abundant species """

    assert metrics.McNaughton([100]) == 'NaN'
    assert metrics.McNaughton([1,1,1,1,1,1,1,1,1,1]) == 20
    assert metrics.McNaughton([0]) == 'NaN'
    assert metrics.McNaughton([2, 6, 76, 1, -1]) == 'NaN'


def test_Shannons_H():

    """ test for the value of Shannon's information entropy """

    assert metrics.Shannons_H([100]) == 0.0
    assert metrics.Shannons_H([1,1,1,1,1,1,1,1,1,1]) == 2.302585

    assert metrics.Shannons_H([0]) == 'NaN'
    assert metrics.Shannons_H([2, 6, 76, 1, -1]) == 'NaN'


def test_simpsons_dom():

    """ test for the value of Shannon's information entropy """

    assert metrics.simpsons_dom([100]) == 0.0
    assert metrics.simpsons_dom([1,1,1,1,1,1,1,1,1,1]) == 0.9
    assert metrics.simpsons_dom([0]) == 'NaN'
    assert metrics.simpsons_dom([2, 6, 76, 1, -1]) == 'NaN'


def test_e_shannon():

    """ test for the value of Pielou's evenness index, which is based on
    Shannon's information entropy """

    assert metrics.e_shannon([100, 100]) == 1.0
    assert metrics.e_shannon([1,1,1,1,1,1,1,1,1,1]) == 1.0
    assert metrics.e_shannon([0]) == 'NaN'
    assert metrics.e_shannon([2, 6, 76, 1, -1]) == 'NaN'


def test_simplest_gini():

    """ test for the value of Gini's inequality index """

    assert metrics.simplest_gini([100, 100]) == 0.0
    assert metrics.simplest_gini([1,1,1,1,1,1,1,1,1,1]) == 0.0
    assert metrics.simplest_gini([1000,1]) == 0.499001
    assert metrics.simplest_gini([0]) == 'NaN'
    assert metrics.simplest_gini([2, 6, 76, 1, -1]) == 'NaN'


def test_e_Mcintosh():

    """ tests for the value of McIntosh's evenness metric """
    assert metrics.e_Mcintosh([100, 100]) == 1.0
    assert metrics.e_Mcintosh([1,1,1,1,1,1,1,1,1,1]) == 1.0
    assert metrics.e_Mcintosh([1000,1]) == 0.003409
    assert metrics.e_Mcintosh([0]) == 'NaN'
    assert metrics.e_Mcintosh([2, 6, 76, 1, -1]) == 'NaN'
