from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import os
import sys

import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table


mydir = os.path.expanduser('~/GitHub/ResTime')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

dat = pd.read_csv(mydir + '/results/simulated_data/examples/2015_12_08/SimData.csv')
dat = dat[np.isfinite(dat['total.abundance'])]
#dat[dat['total.abundance'] == 0] = 1
dat = dat[dat['motion'] == '\'fluid\'']


#N = np.log10(dat['total.abundance']).tolist()
N = dat['total.abundance'].tolist()
AvgGrow = dat['avg.per.capita.growth'].tolist()

MaxRes = dat['max.res.val'].tolist()
ResPart = dat['resource.particles'].tolist()
ResIn = dat['res.inflow'].tolist()
ResConc = dat['resource.concentration'].tolist()

MaxActDisp = dat['max.active.dispersal'].tolist()
AvgActDisp = dat['avg.per.capita.active.dispersal'].tolist()

AvgMaint = dat['avg.per.capita.maint'].tolist()
MaxMaint = dat['max.met.maint'].tolist()

Phase = dat['phase'].tolist()
Freq = dat['frequency'].tolist()
Amp = dat['amplitude'].tolist()

AvgPEff = dat['avg.per.capita.P.efficiency'].tolist()
AvgNEff = dat['avg.per.capita.N.efficiency'].tolist()
AvgCEff = dat['avg.per.capita.C.efficiency'].tolist()

tau = np.log10(dat['width']/dat['flow.rate']).tolist()
tau = np.log10(dat['individual.tau'])
dat['tau'] = list(tau)

d = pd.DataFrame({'N': list(N)})
d['tau'] = list(tau)

#### plot figure ###############################################################
#xlab = r"$log_{10}$"+'(' + r"$\tau$" +')'
fs = 8 # fontsize
fig = plt.figure()

#### N vs. Tau #################################################################

f2 = smf.ols('N ~ tau + I(tau ** 2.0)', d).fit()
f2 = smf.ols('N ~ tau', d).fit()
print f2.summary()
#sys.exit()

st, data, ss2 = summary_table(f2, alpha=0.05)
fitted = data[:,2]

aboveMaxMaint = []
belowMaxMaint = []

aboveAvgGrow = []
belowAvgGrow = []

aboveMaxRes = []
belowMaxRes = []

aboveResPart = []
belowResPart = []

aboveMaxActDisp = []
belowMaxActDisp = []

aboveAvgActDisp = []
belowAvgActDisp = []

aboveAvgMaint = []
belowAvgMaint = []

aboveResIn = []
belowResIn = []

abovePhase = []
belowPhase = []

aboveFreq = []
belowFreq = []

aboveAmp = []
belowAmp = []

aboveResConc = []
belowResConc = []

aboveAvgPEff = []
belowAvgPEff = []

aboveAvgNEff = []
belowAvgNEff = []

aboveAvgCEff = []
belowAvgCEff = []



for i, fval in enumerate(fitted):
    if fval > N[i]:
        belowAvgGrow.append(AvgGrow[i])
        belowResPart.append(ResPart[i])

        belowAvgPEff.append(AvgPEff[i])
        belowAvgNEff.append(AvgNEff[i])
        belowAvgCEff.append(AvgCEff[i])

        belowMaxRes.append(MaxRes[i])
        belowMaxActDisp.append(MaxActDisp[i])

        belowAvgActDisp.append(AvgActDisp[i])
        belowAvgMaint.append(AvgMaint[i])
        belowMaxMaint.append(MaxMaint[i])

        belowResIn.append(ResIn[i])
        belowPhase.append(Phase[i])
        belowFreq.append(Freq[i])
        belowAmp.append(Amp[i])
        belowResConc.append(ResConc[i])


    elif fval < N[i]:
        aboveAvgGrow.append(AvgGrow[i])
        aboveResPart.append(ResPart[i])

        aboveAvgPEff.append(AvgPEff[i])
        aboveAvgNEff.append(AvgNEff[i])
        aboveAvgCEff.append(AvgCEff[i])

        aboveMaxRes.append(MaxRes[i])
        aboveMaxActDisp.append(MaxActDisp[i])

        aboveAvgActDisp.append(AvgActDisp[i])
        aboveAvgMaint.append(AvgMaint[i])
        aboveMaxMaint.append(MaxMaint[i])

        aboveResIn.append(ResIn[i])
        abovePhase.append(Phase[i])
        aboveFreq.append(Freq[i])
        aboveAmp.append(Amp[i])
        aboveResConc.append(ResConc[i])



fig.add_subplot(2, 2, 1)
plt.boxplot([aboveAvgPEff, belowAvgPEff])

fig.add_subplot(2, 2, 2)
plt.boxplot([aboveAvgCEff, belowAvgCEff])

fig.add_subplot(2, 2, 3)
plt.boxplot([aboveAvgNEff, belowAvgNEff])


#fig.add_subplot(2, 2, 1)
#plt.boxplot([aboveAvgGrow, belowAvgGrow])

#fig.add_subplot(2, 2, 1)
#plt.boxplot([aboveMaxActDisp, belowMaxActDisp])

#fig.add_subplot(2, 2, 2)
#plt.boxplot([aboveAvgActDisp, belowAvgActDisp])

#fig.add_subplot(3, 3, 4)
#plt.boxplot([aboveAvgMaint, belowAvgMaint])

#fig.add_subplot(3, 3, 5)
#plt.boxplot([aboveResIn, belowResIn])

#fig.add_subplot(2, 2, 3)
#plt.boxplot([aboveMaxMaint, belowMaxMaint])

#fig.add_subplot(3, 3, 7)
#plt.boxplot([aboveFreq, belowFreq])

#fig.add_subplot(3, 3, 8)
#plt.boxplot([aboveMaxRes, belowMaxRes])

#fig.add_subplot(2, 2, 4)
#plt.boxplot([aboveResConc, belowResConc])
#plt.ylim(0,2)

#plt.ylabel(r"$log_{10}$"+'(' + r"$N$" +')', fontsize=fs+6)
#plt.xlabel(xlab, fontsize=fs+6)
#plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.ylim(-2, 5)

#### Final Format and Save #####################################################
#plt.subplots_adjust(wspace=0.4, hspace=0.4)
#plt.savefig(mydir2 + 'Desktop/ResTime/figures/Fig3.png', dpi=600, bbox_inches = "tight")
plt.show()
