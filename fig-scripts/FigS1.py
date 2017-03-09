from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

mydir = os.path.expanduser('~/GitHub/residence-time')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df2 = pd.DataFrame({'width' : df['width']})
df2['flow'] = df['flow.rate']
df2['tau'] = np.log10((df2['width']**3)/df2['flow'])

df2['R'] = np.log10(df['resource.particles'])
df2['RDens'] = np.log10(df['resource.concentration'])

#### plot figure ###############################################################
gd = 20
mnct = 0
xlab = r"$log_{10}$"+'(' + r"$\tau$" +')'
fs = 6 # fontsize
fig = plt.figure()


#### R vs. Tau #################################################################
fig.add_subplot(2, 2, 1)

plt.hexbin(df2['tau'], df2['R'], mincnt=mnct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(r"$log_{10}$"+'(' + r"$R$" +')', fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(1.1, 2, 'A', color = 'y', fontweight='bold')

#### RDens vs. Tau ########################################################
#dat = dat.convert_objects(convert_numeric=True).dropna()
fig.add_subplot(2, 2, 2)

plt.hexbin(df2['tau'], df2['RDens'], mincnt=mnct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(r"$log_{10}$"+'(Resource concentration)', fontsize=fs+3)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.text(4.2, -0.25, 'B', color = 'y', fontweight='bold')

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/FigS1.png', dpi=600, bbox_inches = "tight")
#plt.show()
plt.close()
