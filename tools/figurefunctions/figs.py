from __future__ import division

import  matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

import numpy as np
import random
import os
import sys

mydir = os.path.expanduser('~/Desktop/Repos/HYDRO-BIDE')
sys.path.append(mydir+'/tools')


def obs_pred_rsquare(obs, pred):
    """Determines the prop of variability in a data set accounted for by a model
        In other words, this determines the proportion of variation explained by
        the 1:1 line in an observed-predicted plot. """
    
    return 1 - sum((obs - pred) ** 2) / sum((obs - np.mean(obs)) ** 2)
    


def figSet(fig, i, var_to_vary, varlist, lists, TauList, x_label, y_label, colors):
    
    fig.add_subplot(2,2,i+1)
    fs = 14 # fontsize
    
    ymax = 0
    
    for j, _list in enumerate(lists):
         
        plt.plot(TauList[j], _list, linewidth=2, c=colors[j], alpha=0.8)
        if max(_list) > ymax: ymax = max(_list)
        
    #plt.xlim(1.0, max(TauList[j]))
    plt.ylim(0.0, ymax*1.25)
                                                                            
    if y_label == 'Avg abundance, N/S':# or y_label == 'Richness, S':
        plt.ylim(1.0, ymax*1.25)
        plt.yscale('log')
        
    plt.xscale('log')
    plt.xlabel(x_label, fontsize=fs+2)
    plt.ylabel(y_label, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-2)
    plt.subplots_adjust(wspace=0.35, hspace=0.35)
    
    return fig


def figPred(fig, i, var_to_vary, varlist, lists, TauList, x_label, y_label, tauAtQmaxList, tauAtAvgMuMaxList, colors):
    
    fig.add_subplot(2,2,i+1)
    fs = 14 # fontsize
    ymax = 0
    xmax = 0
    for j, _list in enumerate(lists):
        
        if xmax < max(TauList[j]): xmax = max(TauList[j])
        
        plt.plot(TauList[j], _list, linewidth=2, c=colors[j], alpha=0.8)
        if i == 0: 
            x = tauAtAvgMuMaxList[j]
            y = TauList[j].index(x)
            if x != TauList[j][y]:
                print 'AvgMuMax plot error: x != TauList[j]'
            plt.scatter(x, _list[y], s=60, facecolor=colors[j], edgecolor=colors[j], alpha=0.8)                        
            plt.xlim(1.0, xmax)
            plt.yscale('log')
    
        elif i == 1:
            x = tauAtQmaxList[j]
            y = TauList[j].index(x)
            if x != TauList[j][y]:
                print 'tauAtQmax plot error: x != TauList[j]'
            plt.scatter(x, _list[y], s=60, facecolor=colors[j], edgecolor=colors[j], alpha=0.8)                        
            plt.xlim(1.0, xmax)
            
        if max(_list) > ymax: ymax = max(_list)
        
    plt.ylim(0.0, ymax*1.25)                                                                            
    plt.xscale('log')
    
    plt.xlabel(x_label, fontsize=fs+2)
    plt.ylabel(y_label, fontsize=fs)
    
    plt.tick_params(axis='both', which='major', labelsize=fs-2)
    plt.subplots_adjust(wspace=0.35, hspace=0.35)
    
    return fig



def ObsPredAvgMuMax(fig, ObstauAtAvgMuMaxList, tauAtAvgMuMaxList, colors):

    fig.add_subplot(2,2,3)
    fs = 14 # fontsize
    axesMax = 0
    
    for j, obs in enumerate(ObstauAtAvgMuMaxList):
        x = np.log(ObstauAtAvgMuMaxList[j])
        y = np.log(tauAtAvgMuMaxList[j])
        plt.scatter(x, y, s=60, facecolor=colors[j], edgecolor=colors[j], alpha=0.8)                        
        axesMax = max([x, y, axesMax])
    
        axesMax = max([x, y, axesMax])     
    
    x2 = [0, axesMax*2]
    y2 = [0, axesMax*2]
    y2 = [axesMax*2, 0]
                            
    plt.plot(x2, y2, color='k', lw=1)
    plt.xlim(0, axesMax*2)
    plt.ylim(0, axesMax*2)
    
    plt.xlabel('log(obs greatest avg gmax)', fontsize=fs-2)
    plt.ylabel('log(exp greatest avg gmax)', fontsize=fs-2)
    plt.tick_params(axis='both', which='major', labelsize=fs-2)
    plt.subplots_adjust(wspace=0.35, hspace=0.35)
    
    return fig



def NvsTau(fig, r, c, i, wsp, hsp, Nlist, REStimes, mean, tauAtQmaxList, colorBand=False):
    
    #y_label = 'Total abundance'
    y_label = 'Biomass turnover'
    x_label = r"$\tau$"
    
    fig.add_subplot(r, c, i)
    fs = 14 # fontsize
    colors = ['0.0','0.3','0.6']
    
    ymax = 0
    x3max = 1.0
    
    for j, _list in enumerate(Nlist):  
        plt.scatter(REStimes[j], _list, linewidth=0.0, c=colors[j], alpha=0.8)
        
        if max(_list) > ymax: 
            ymax = max(_list)
        if max(REStimes[j]) > x3max:
            x3max = max(REStimes[j])
            
    if i+1 == 1:
        plt.legend(bbox_to_anchor=(-0.03, 1.1, 2.57, .2), loc=10, ncol=3, mode="expand",prop={'size':12})        
    
    if colorBand == True:
        x1 = [1.0, np.mean(tauAtQmaxList)]
        plt.fill_between(x1, 0, 10**10, facecolor='red', interpolate=True, alpha = 0.1)                        
    
        x3 = [mean**-1, 10**5]
        plt.fill_between(x3, 0, 10**10, facecolor='yellow', interpolate=True, alpha = 0.1)                        
    
    plt.xlim(1.0, x3max)
    plt.ylim(0.0, ymax*1.25)
                                                                                
    plt.xscale('log')
    plt.xlabel(x_label, fontsize=fs)
    plt.ylabel(y_label, fontsize=fs)
    
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.subplots_adjust(wspace= wsp, hspace= hsp)
    
    return fig
    

def MCTvsTau(fig, r, c, i, wsp, hsp, MCTlist, REStimes, tauAtQmaxList):
    
    # plot mean cell residence time vs. ecosystem residence time
    
    y_label = 'MCT'
    x_label = r"$\tau$"
    
    fig.add_subplot(r, c, i)
    fs = 14 # fontsize
    
    ymax = 0
    x3max = 1.0
    
    obs = []
    pred = []
    
    for j, _list in enumerate(MCTlist):  
        
        for ii, dot in enumerate(_list):
            
            random.seed()    
            r = lambda: random.randint(0,255)
            color = '#%02X%02X%02X' % (r(),r(),r())
            plt.scatter(REStimes[j][ii], dot, linewidth=0.0, c= color, s=25, alpha=0.85)
            
            obs.append(REStimes[j][ii])
            pred.append(dot)    
            
        if max(_list) > ymax: 
            ymax = max(_list)
        if max(REStimes[j]) > x3max:
            x3max = max(REStimes[j])
            
    plt.plot([1,x3max],[1,x3max], color='gray')
    R2 = obs_pred_rsquare(np.array(obs), np.array(pred))
    
    textymax = np.exp(0.8 * float(np.log(x3max)))
    
    print 'r-square:',R2
    plt.text(1.4, textymax, r'$R^2$' + '= '+str(round(R2,3)), fontsize=fs)
    
    plt.xlim(1.0, x3max)
    plt.ylim(1.0, x3max)
                                                                                
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(x_label, fontsize=fs)
    plt.ylabel(y_label, fontsize=fs)
    
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.subplots_adjust(wspace= wsp, hspace= hsp)
    
    return fig
    


def MUvsD(fig, r, c, i, wsp, hsp, Plist, REStimes, tauAtQmaxList):
    
    y_label = 'Growth rate'
    x_label = r"$\tau^{-1}$"
    
    fig.add_subplot(r, c, i)
    
    fs = 14 # fontsize
    
    # plot 3: growth rate vs. dilution rate
    obs = []
    pred = []
    
    for j, _list in enumerate(Plist):  
        
        Ds = []
        for tau in REStimes[j]:
            Ds.append(tau**-1)
        
        for ii, dot in enumerate(_list):        
            random.seed()    
            r = lambda: random.randint(0,255)
            color = '#%02X%02X%02X' % (r(),r(),r())
            plt.scatter(Ds[ii], dot, linewidth=0.0, c= color, s=25, alpha=0.85)
            
            obs.append(Ds[ii])
            pred.append(dot)
            
    plt.plot([0,1],[0,1], color='gray')
    R2 = obs_pred_rsquare(np.array(obs), np.array(pred))
    
    print 'r-square:',R2
    plt.text(0.1, 0.8, r'$R^2$' + '= '+str(round(R2,3)), fontsize=fs)
    
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
                                                                                
    plt.xlabel(x_label, fontsize=fs)
    plt.ylabel(y_label, fontsize=fs)
    
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.subplots_adjust(wspace= wsp, hspace= hsp)
    
    return fig
    
    
    
def BTvsTau(fig, r, c, i, wsp, hsp, BTlist, REStimes, tauAtQmaxList):

    # plot mean cell residence time vs. ecosystem residence time
    
    y_label = 'Biomass turnover time'
    x_label = r"$\tau$"
    
    fig.add_subplot(r, c, i)
    fs = 14 # fontsize
    
    ymax = 0
    x3max = 1.0
    
    obs = []
    pred = []
    
    for j, _list in enumerate(BTlist):  
        
        for ii, dot in enumerate(_list):
            
            dot = 1/dot
            random.seed()    
            r = lambda: random.randint(0,255)
            color = '#%02X%02X%02X' % (r(),r(),r())
            plt.scatter(REStimes[j][ii], dot, linewidth=0.0, c= color, s=40, alpha=0.9)
            
            obs.append(REStimes[j][ii])
            pred.append(dot)    
            
        if max(_list) > ymax: 
            ymax = min(_list)
        if max(REStimes[j]) > x3max:
            x3max = max(REStimes[j])
            
    #plt.plot([1,x3max],[1,x3max], color='gray')
    #R2 = obs_pred_rsquare(np.array(obs), np.array(pred))
    
    #textymax = np.exp(0.8 * float(np.log(x3max)))
    
    #print 'r-square:',R2
    #plt.text(1.4, textymax, r'$R^2$' + '= '+str(round(R2,3)), fontsize=fs)
                                                                                
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(x_label, fontsize=fs)
    plt.ylabel(y_label, fontsize=fs)
    
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.subplots_adjust(wspace= wsp, hspace= hsp)
    
    return fig
    
    

def TTvsTau(fig, r, c, i, wsp, hsp, BTlist, REStimes, tauAtQmaxList):  # Turnover time

    # plot mean cell residence time vs. ecosystem residence time
    
    y_label = 'Turnover time'
    x_label = r"$\tau$"
    
    fig.add_subplot(r, c, i)
    fs = 14 # fontsize
    
    ymax = 0
    x3max = 1.0
    
    obs = []
    pred = []
    
    for j, _list in enumerate(BTlist):  
        
        for ii, dot in enumerate(_list):
            
            random.seed()    
            r = lambda: random.randint(0,255)
            color = '#%02X%02X%02X' % (r(),r(),r())
            plt.scatter(REStimes[j][ii], dot, linewidth=0.0, c= color, s=25, alpha=0.85)
            
            obs.append(REStimes[j][ii])
            pred.append(dot)    
            
        if max(_list) > ymax: 
            ymax = max(_list)
        if max(REStimes[j]) > x3max:
            x3max = max(REStimes[j])
            
    #plt.plot([1,x3max],[1,x3max], color='gray')
    #R2 = obs_pred_rsquare(np.array(obs), np.array(pred))
    
    #textymax = np.exp(0.8 * float(np.log(x3max)))
    
    #print 'r-square:',R2
    #plt.text(1.4, textymax, r'$R^2$' + '= '+str(round(R2,3)), fontsize=fs)
                                                                                
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(x_label, fontsize=fs)
    plt.ylabel(y_label, fontsize=fs)
    
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.subplots_adjust(wspace= wsp, hspace= hsp)
    
    return fig
    

def CTvsTau(Slists, REStimes):
    
    """
    Slists  :  A list of lists (subslists).  
    
    
    """
    
    #for Slists 
    #taxa_time_relationship = get_STR(Slist)
    return
