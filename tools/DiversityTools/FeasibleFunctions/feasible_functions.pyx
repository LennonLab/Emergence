from __future__ import division

import sys
import os
from os import path, access, R_OK  # W_OK for write permission

mydir = os.path.expanduser("~/GitHub/rare-bio/")
mydir2 = os.path.expanduser("~/Desktop/")

sys.path.append(mydir + "/tools/macroecotools")
import macroecotools
sys.path.append(mydir + "/tools/partitions")
import partitions as parts

import  matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde

import random
from random import choice
import decimal
import re

import math

########################################################################################################
######   A Section devoted to evenness indices and descriptive statistical functions ###################

"""
    Modified from Locey and White (2013)
    (insert reference)

"""
def min_max(n,s):

    _min = int(math.floor(float(n)/float(s)))
    if int(n%s) > 0:
        _min +=1

    return _min


def get_modes(_list,which):

    modes = []
    for i in _list:
        _mode = np.mode(i)
        if which is 'high':
            modes.append(int(max(_mode)))
        if which is 'low':
            modes.append(int(min(_mode)))

    return modes


def get_modal(_list):

    """ Finds the mode from a kernel density function across a sample """
    exp_mode = 0.0
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list),n)
    density.covariance_factor = lambda : .001
    density._compute_covariance()
    D = [xs,density(xs)]
    d = 0
    maxd = 0.0
    while d < len(D[1]):
        if D[1][d] > maxd:
            maxd = D[1][d]
            exp_mode = D[0][d]
        d += 1
    return exp_mode

def get_kdens_choose_kernel(_list,kernel):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D

def get_kdens(_list):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    #xs = np.linspace(min(_list),max(_list),n)
    xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : 0.5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D


#######################################################################################################
#### A section devoted to finding constraint combinations and empirical DOWs/SADs from datasets #######

def GetObsPred(dataset, kind):

    PATH = mydir2 + 'data/' + kind + '/' + dataset + '/' + dataset + '_obs_pred.txt'
    if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
        DATA = open(PATH,'r')

        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group()
        m2 = int(re.findall(r'\d*\S$',d)[0])
        SAD = [int(m2)]
        SADs = []

        for d in DATA:
            ct1+=1
            m1 = re.match(r'\S*',d).group()
            if m1 == m0:
                m2 = int(re.findall(r'\d*\S$',d)[0])
                if m2 > 0:SAD.append(m2)
            else:
                site_name = m0
                m0 = m1
                if len(SAD) > 1 and sum(SAD) <= 100000:
                    SAD.sort()
                    SAD.reverse()
                    SADs.append(SAD) # can also append, site_name, len(SAD), and sum(SAD)
                    ct2+=1
                SAD = []
                abundance = int(re.findall(r'\d*\S$',d)[0])
                if abundance > 0:SAD.append(abundance)
        DATA.close()
        return(SADs)


def get_SADs(dataset):

    DATA = open('/home/kenlocey/data1/' + dataset + '/' + dataset + '-data.txt','r')
    ct1 = 0
    ct2 = 0
    d = DATA.readline()
    m0 = re.match(r'\A\S*',d).group()
    m2 = int(re.findall(r'\d*\S$',d)[0])
    SAD = [int(m2)]
    SADs = []

    for d in DATA:
        ct1+=1
        m1 = re.match(r'\A\S*',d).group()
        if m1 == m0:
            m2 = int(re.findall(r'\d*\S$',d)[0])
            if m2 > 0:
                SAD.append(m2)

        else:
            site_name = m0
            m0 = m1
            if len(SAD) > 9 and sum(SAD) <= 100000:
                SAD.sort()
                SAD.reverse()
                SADs.append(SAD) # can also append, site_name, len(SAD), and sum(SAD)  !!THIS NEEDS TO BE DEALT WITH!!
                ct2+=1
            SAD = []
            abundance = int(re.findall(r'\d*\S$',d)[0])
            if abundance > 0:SAD.append(abundance)
    DATA.close()
    return(SADs)





########################################################################################################
######   A Section devoted to finding macrostates/integer partitions ############################

def get_random_macrostates(NS_combo): # This function call rand_parts (derived by KJL), a function
    # that finds random macrostates much faster than Sage. Not used for Locey and White (2013).
    N = int(NS_combo[0])
    S = int(NS_combo[1])
    sample_size = 63
    return parts.rand_parts1(N,S,sample_size)



########################################################################################################
##### A Section devoted to examining randoms samples of feasibles sets and empirical data ##############


def get_hottest_SAD(unique_SADs):
    """ Find the SAD in a random sample with the greatest average commonness
        among its ranked abundance states. This SAD is taken to represent the
        central tendency of the set, based on the SAD shape. """

    if len(unique_SADs) > 500:
        unique_SADs = random.sample(unique_SADs,500)

    N = sum(unique_SADs[0])
    S = len(unique_SADs[0])
    a1 = 0 # SAD mean
    v1 = 0 # SAD variance
    for rad in unique_SADs:
        in_common = []
        ct1 = 0
        for a in rad: # for each rank
            c = 0
            for sad in unique_SADs:
                if a == sad[ct1]:
                    c += 1
            in_common.append(np.log(c))
            ct1 += 1
        a2 = np.mean(in_common)
        v2 = np.var(in_common)
        if a2 > a1:
            a1 = a2
            v1 = v2
            xRAD = rad
        elif a2 == a1:
            if v2 < v1:
                a1 = a2
                v1 = v2
                xRAD = rad
    #percentile_evar = stats.percentileofscore(sample_evar,obs_evar)
    return xRAD



def import_obs_pred_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    data = np.genfromtxt(input_filename, dtype = "S15,f8,f8", names = ['site','obs','pred'], delimiter = " ")
    return data


def hist_mete_r2(sites, obs, pred):  # TAKEN FROM Macroecotools or the mete_sads.py script used for White et al. (2012)
    """Generate a kernel density estimate of the r^2 values for obs-pred plots"""
    r2s = []
    for site in np.unique(sites):
        obs_site = obs[sites==site]
        pred_site = pred[sites==site]
        r2 = macroecotools.obs_pred_rsquare(obs_site, pred_site)
        #print site,r2
        r2s.append(r2)
    hist_r2 = np.histogram(r2s, range=(0, 1))
    xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
    xvals = xvals[0:len(xvals)-1]
    yvals = hist_r2[0]
    plt.plot(xvals, yvals, 'k-', linewidth=2)
    plt.axis([0, 1, 0, 1.1 * max(yvals)])


def obs_pred_r2_multi(datasets, data_dir='/home/kenlocey/data1/'): # TAKEN FROM THE mete_sads.py script
    print 'generating 1:1 line R-square values for dataset(s)'
    for i, dataset in enumerate(datasets):
        obs_pred_data = import_obs_pred_data(data_dir + dataset + '/' + dataset + '_obs_pred.txt')
        obs = ((obs_pred_data["obs"]))
        pred = ((obs_pred_data["pred"]))
        print dataset,' ',macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))


def plot_obs_pred_sad(datasets, data_dir='/home/kenlocey/data1/', radius=2): # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    # Figure 3 Locey and White (2013)        ##########################################################################################################

    """Multiple obs-predicted plotter"""
    fig = plt.figure()
    I = [1,2,5,6,9,10,13,14]
    xs = [[60,1], [100,1], [20,1], [60,1], [40,1], [200,1], [800,1.5], [200,1.5]]
    rs = ['0.93','0.77','0.84','0.81','0.78','0.83','0.58','0.76']
    for i, dataset in enumerate(datasets):
        i = I[i]
        print dataset
        obs_pred_data = import_obs_pred_data(data_dir + dataset + '/' + dataset + '_obs_pred.txt')
        site = ((obs_pred_data["site"]))
        obs = ((obs_pred_data["obs"]))
        pred = ((obs_pred_data["pred"]))

        axis_min = 0.5 * min(obs)
        axis_max = 2 * max(obs)
        ax = fig.add_subplot(4,4,i+1)
        macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1,
                                            plot_obj=plt.subplot(4,4,i+1))
        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        plt.xlim(axis_min, axis_max)
        plt.ylim(axis_min, axis_max)
        plt.tick_params(axis='both', which='major', labelsize=8)
        r2 = 0.8
        plt.subplots_adjust(wspace=0.5, hspace=0.3)
        plt.text(xs[0][1],xs[0][0],dataset+'\n'+rs[0],fontsize=8)
        xs.pop(0)
        rs.pop(0)
        # Create inset for histogram of site level r^2 values
        axins = inset_axes(ax, width="30%", height="30%", loc=4)
        hist_mete_r2(site, np.log10(obs), np.log10(pred))
        plt.setp(axins, xticks=[], yticks=[])

    plt.text(-8,-80,'Rank-abundance at the centre of the feasible set',fontsize=10)
    plt.text(-8.5,500,'Observed rank-abundance',rotation='90',fontsize=10)
    plt.savefig('obs_pred_plots.png', dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)




def histOutline(dataIn, *args, **kwargs): # found at http://www.scipy.org/Cookbook/Matplotlib/UnfilledHistograms
    (histIn, binsIn) = np.histogram(dataIn, *args, **kwargs)

    stepSize = binsIn[1] - binsIn[0]
    bins = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]
    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0

    return (bins, data)



def next_partition(p):
  if max(p) == 1:
    #return [sum(p)]
    return [0]
  p.sort()
  p.reverse()
  q = [ p[n] for n in range(len(p)) if p[n] > 1 ]
  q[-1] -= 1
  if (p.count(1)+1) % q[-1] == 0:
    return q + [q[-1]]*((p.count(1)+1) / q[-1])
  else:
    return q + [q[-1]]*((p.count(1)+1) / q[-1]) + [(p.count(1)+1) % q[-1]]
