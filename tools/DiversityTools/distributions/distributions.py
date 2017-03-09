from __future__ import division
import sys
import os
import numpy as np
from scipy import stats
import  matplotlib.pyplot as plt

import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd

########### PATHS ##############################################################

tools = os.path.expanduser("~/tools")
sys.path.append(tools + "/macroeco_distributions")
import macroeco_distributions as md

######### CLASS ################################################################

class zipf:

    """ A class to obtain a zipf object with inherited mle shape parameter,
    mle form of the rank-abundance distribution, and a rank-abundance curve
    based on fitting the zipf to the observed via a generalized linear model."""

    def __init__(self, obs):
        self.obs = obs


    def from_cdf(self):
        """ Obtain the maximum likelihood form of the Zipf distribution, given
        the mle value for the Zipf shape parameter (a). Using a, this code
        generates a rank-abundance distribution (RAD) from the cumulative
        density function (cdf) using the percent point function (ppf) also known
        as the quantile function.
        see: http://www.esapubs.org/archive/ecol/E093/155/appendix-B.htm

        This is an actual form of the Zipf distribution, obtained from getting
        the mle for the shape parameter.
        """

        p = md.zipf_solver(self.obs)
        S = len(self.obs)
        rv = stats.zipf(a=p)
        rad = []
        for i in range(1, S+1):
            val = (S - i + 0.5)/S
            x = rv.ppf(val)
            rad.append(int(x))

        return rad


    def from_glm(self):

        """ Fit the Zipf distribution to the observed vector of integer values
        using a generalized linear model.

        Note: This is a fitted curve; not an actual form of the Zipf distribution

        This method was inspired by the vegan
        package's open source code on vegan's public GitHub repository:
        https://github.com/vegandevs/vegan/blob/master/R/rad.zipf.R
        on Thursday, 19 Marth 2015 """

        ranks = np.log(range(1, len(self.obs)+1))
        off = [np.log(sum(self.obs))] * len(self.obs)

        d = pd.DataFrame({'ranks': ranks, 'off': off, 'x':self.obs})

        lm = smf.glm(formula='x ~ ranks', data = d, family = sm.families.Poisson()).fit()
        pred = lm.predict()

        return pred


# zipf_pred = zipf(ad)
# zipf_mle = zipf_pred.from_cdf()
# zipf_glm = zipf_pred.from_glm()
