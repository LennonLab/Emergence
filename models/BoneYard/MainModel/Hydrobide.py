from __future__ import division
import os
import sys
import numpy as np
from random import randrange, randint,  choice, sample

mydir = os.path.expanduser('~/Desktop/Repos/HYDRO-BIDE')
sys.path.append(mydir+'/tools')
import metrics
import timeSeries

sys.path.append(mydir+'/models/NoHydro/bide')
import bide






def get_params(system, REStimes):

        V = 1
        lgp = choice([0.9, 0.92, 0.94, 0.96, 0.98, 0.999])

        Kq = randrange(1,10)/1000
        PropQ = choice([0.05, 0.1, 0.15, 0.2])

        uptake_mu = randrange(500, 1000)/1000
        uptake_std = randrange(1, 100)/1000

        Pcons = 0
        Rcons = randrange(10, 100)

        env_noise = 0
        spatial_het = 0

        dkern_std = (float(randrange(10**5))/10**5)*V
        dkern_mu = 0

        res_pulse, im_pulse = [1, 1]
        maint = 0.0 #(randrange(1, 10)/1000)

        if system != 'ideal':

            Pcons = randrange(1, 50)

        if system == 'complex_ecosystem':
            res_pulse = randint(1,5)
            im_pulse = randint(1,5)

            env_noise = 1
            spatial_het = 1

        return [V, lgp, Kq, PropQ, Pcons, Rcons, uptake_mu, uptake_std, maint, res_pulse, im_pulse, env_noise, spatial_het, dkern_mu, dkern_std]





def Hbide(V, REStime, Rcons, Pcons, PropQ, time, mean, std, Kq, lgp, maint, res_pulse, im_pulse, env_noise, spatial_het, basic, system, disp_mu, disp_std):

    """
    V         :  Volume

    Rcons     :  resource concentration

    Pcons     :  propagule concentration

    maint     :  resource loss due to cell maintenance

    REStime   :  residence time

    time      :  number of time steps to run simulation. Needs to be large
    enough to allow the community to reach a relatively stable state

    dormancy  :  designate whether dormancy is allowed. If dormancy is True,
    then cell quotas of 0 result in dormancy, if dormancy = False, then cell
    quotas of 0 result in death and immediate 'removal'

    """

    r = V/REStime # rate of flow
    TR = Rcons * V # total resources in the local environment


    num_in = int(round(r * Pcons))
    if num_in < 1 and Pcons > 1:
        return ['Continuously empty system']

    Nlist = [] # total abundance
    Glist = [] # average growth rate
    Slist = []  # richness
    Prodlist = [] # productivity list
    Rlist = [] # total resources
    Qlist = [] # cell quota
    MAXUlist = []
    Tlist = [] # number of generations in the system
    MCT = []

    if basic == False:
        SoNlist = [] # average species abundance
        TOlist = [] # biomass turnover
        Evarlist = [] # evenness
        Hlist = [] # Shannon's diversity
        CTlist = [] # compositional turnover
        RADlist = [] # list of RADs

    maxUptake_list = []
    maxUptake_dict = {}

    xcoords =  []
    ycoords =  []
    zcoords =  []

    DispParamsDict = {}

    Qs1 = []
    Qs2 = list(Qs1)
    N = len(Qs1)

    Slist1 = []
    Slist2 = list(Slist1)

    inds1 = []
    inds2 = list(inds1)

    ID = 0
    lag = 0

    """ Initial bout of immigration """
    init_n = 1000
    Qs2, Slist2, inds2, Tlist, ID, maxUptake_dict, maxUptake_list, DispParamsDict, xcoords, ycoords, zcoords = bide.immigration(system, V, Qs2, Slist2, inds2, Tlist, ID, maxUptake_dict, maxUptake_list, DispParamsDict, mean, std, lgp, init_n, PropQ, xcoords, ycoords, zcoords, disp_mu, disp_std)

    Nlist2 = []
    pastBurnIn = 'no'

    Hurst = float()
    pval = float()

    for t in range(time):

        """ Resource resupply """
        TR += r * Rcons # Total Resources increase due to inflow

        """ Immigration """
        if system != 'ideal':
            Qs2, Slist2, inds2, Tlist, ID, maxUptake_dict, maxUptake_list, DispParamsDict, xcoords, ycoords, zcoords = bide.immigration(system, V, Qs2, Slist2, inds2, Tlist, ID, maxUptake_dict, maxUptake_list, DispParamsDict, mean, std, lgp, num_in, PropQ, xcoords, ycoords, zcoords, disp_mu, disp_std)


        """ Growth """
        Qs2, Slist2, inds2, maxUptake_list, Tlist, TR, glist, xcoords, ycoords, zcoords = bide.growth(Qs2, Slist2, inds2, maxUptake_list, Tlist, TR, Kq, maint, xcoords, ycoords, zcoords)

        N = len(Qs2)
        if N == 0:
            return ['N = 0']

        """ Reproduction """
        Qs2, Slist2, inds2, Tlist, maxUptake_list, Kq, ID, prod, xcoords, ycoords, zcoords = bide.reproduction(Qs2, Slist2, inds2, Tlist, maxUptake_list, Kq, ID, maint, xcoords, ycoords, zcoords)


        """ Check community variables and list lengths"""
        check = checks(Qs2, Slist2, inds2, Tlist, xcoords, ycoords, zcoords)
        if check[0] != 'pass':
            return check


        """ Emigration """ # will eventually be modified to allow individuals to naturally flow out
        Qs2, Slist2, inds2, Tlist, maxUptake_list, xcoords, ycoords, zcoords = bide.emigration(Qs2, Slist2, inds2, Tlist, maxUptake_list, r, V, xcoords, ycoords, zcoords, system)

        N = len(Qs2)
        if N == 0:
            return ['N = 0 after emigration']


        """ Dispersal """
        Qs2, Slist2, inds2, Tlist, maxUptake_list, DispParamsDict, Kq, ID, xcoords, ycoords, zcoords = bide.dispersal(V, Qs2, Slist2, inds2, Tlist, maxUptake_list, DispParamsDict, Kq, ID, xcoords, ycoords, zcoords, system)

        """ Outflow """
        TR, V = bide.outflow(TR, r, V)

        N = len(Qs2)
        if N == 0:
            return ['N = 0']


        if pastBurnIn == 'no':

            Nlist2.append(N)
            if len(Nlist2) > 500:

                Hurst, pval = timeSeries.get_burnin(Nlist2)
                if Hurst < 0.5 and pval < 0.01:
                    pastBurnIn = 'yes'
                    Nlist2 = []
                    burnIn = t
                    time += time # let the simulation go longer

                else:
                    pastBurnIn = 'no'


        if pastBurnIn == 'yes':

            """ Record community info """

            Glist.append(np.mean(glist))
            Nlist.append(prod)
            Slist.append(len(set(Slist2)))
            Prodlist.append(prod)
            Rlist.append(float(TR))
            Qlist.append(float(sum(Qs2)/len(Qs2)))
            MAXUlist.append(float(np.mean(maxUptake_list)))

            if min(Tlist) == 0 and max(Tlist) == 0:
                MCT.append(0)
            else:
                Tlist_noZeros = filter(lambda a: a != 0, Tlist)
                MCT.append(float(np.mean(Tlist_noZeros)))

            if basic == False:

                S1 = set(Slist1)
                s1 = len(S1)
                S2 = set(Slist2)
                s2 = len(S2)
                c = len(S1 & S2)
                CT = ((s1 - c) + (s2 - c))/np.mean([s1, s2])
                # Whittaker's tunover (species)

                I1 = set(inds1)
                i1 = len(I1)
                I2 = set(inds2)
                i2 = len(I2)
                c = len(I1 & I2)
                TO = ((i1 - c) + (i2 - c))/np.mean([i1, i2])
                # Whittaker's tunover (individuals)

                RAD = [] # Generate the rank-abundance distribution (RAD)
                for i in S2:
                    ab = Slist2.count(i)
                    RAD.append(ab)

                RAD.sort()
                RAD.reverse()

                RADtaulist = [10, 100, 1000]
                ct = RADtaulist.count(REStime)
                if ct > 0 and t == time-1:
                    RADlist = list([V, REStime, RAD])

                TOlist.append(TO)
                SoNlist.append(float(N/s2))

                Hlist.append(float(metrics.Shannons_even(RAD)))
                Evarlist.append(float(metrics.e_var(RAD)))
                CTlist.append(float(CT))


            if len(Nlist) >= 500:

                #lag = timeSeries.get_uncorrelated(Nlist)
                if REStime <= 10:
                    lag = 10
                else: lag = REStime

                if lag > 0:
                    samp_size = int(round(len(Nlist)/lag))

                    if samp_size >= 10:

                        Glist = sample(Glist, samp_size)
                        Nlist = sample(Nlist, samp_size)
                        Slist = sample(Slist, samp_size)
                        Prodlist = sample(Prodlist, samp_size)
                        MAXUlist = sample(MAXUlist, samp_size)
                        Rlist = sample(Rlist, samp_size)
                        Qlist = sample(Qlist, samp_size)
                        MCT = sample(MCT, samp_size)

                        if basic == True:
                            return [np.mean(Glist), np.mean(Nlist), np.mean(Slist),
                                    np.mean(Rlist), np.mean(Qlist), np.mean(MAXUlist),
                                    burnIn, lag, np.mean(Prodlist), np.mean(MCT),
                                    Hurst, pval]

                        elif basic == False:
                            TOlist = sample(TOlist, samp_size)
                            SoNlist = sample(SoNlist, samp_size)
                            Hlist = sample(Hlist, samp_size)
                            Evarlist = sample(Evarlist, samp_size)
                            CTlist = sample(CTlist, samp_size)

                            return [np.mean(Glist), np.mean(Nlist), np.mean(Rlist), np.mean(TOlist),
                                    np.mean(Qlist), np.mean(CTlist),
                                    np.mean(Evarlist), np.mean(Hlist), RADlist,
                                    np.mean(SoNlist), np.mean(MAXUlist), burnIn,
                                    lag, np.mean(Prodlist), np.mean(MCT), Hurst, pval]


        """ saving this generations lists """
        Slist1 = list(Slist2)
        inds1 = list(inds2)
        Qs1 = list(Qs2)

        if t == time-1 and pastBurnIn == 'no':
            return ['failed to reach stationarity in '+str(time)+' generations, H = '+str(round(Hurst, 3))+', p = '+str(round(pval, 3))]


    #lag = timeSeries.get_uncorrelated(Nlist)
    if REStime <= 10:
        lag = 10
    else: lag = REStime


    if lag == 0:
        return ['no sufficient time lag in '+str(len(Nlist))+' generations']

    else:
        samp_size = int(round(len(Nlist)/lag))
        if samp_size < 10:
            return ['insufficient sample of unbiased time points from '+str(len(Nlist))+' generations']

        Glist = sample(Glist, samp_size)
        Nlist = sample(Nlist, samp_size)
        Slist = sample(Slist, samp_size)
        Prodlist = sample(Prodlist, samp_size)
        MAXUlist = sample(MAXUlist, samp_size)
        Rlist = sample(Rlist, samp_size)
        Qlist = sample(Qlist, samp_size)
        MCT = sample(MCT, samp_size)

        if basic == True:
            return [np.mean(Glist), np.mean(Nlist), np.mean(Slist), np.mean(Rlist),
                    np.mean(Qlist), np.mean(MAXUlist), burnIn,
                    lag, np.mean(Prodlist), np.mean(MCT), Hurst, pval]

        elif basic == False:
            TOlist = sample(TOlist, samp_size)
            SoNlist = sample(SoNlist, samp_size)
            Hlist = sample(Hlist, samp_size)
            Evarlist = sample(Evarlist, samp_size)
            CTlist = sample(CTlist, samp_size)

            return [np.mean(Glist), np.mean(Nlist), np.mean(Rlist),
                    np.mean(TOlist), np.mean(Qlist), np.mean(CTlist),
                    np.mean(Evarlist), np.mean(Hlist), RADlist,
                    np.mean(SoNlist), np.mean(MAXUlist), burnIn,
                    lag, np.mean(Prodlist), np.mean(MCT), Hurst, pval]


    return ['no result']






def run_hydrobide(use_picloud, system, Tau, basic, replicates, time):


    """ Properties:
    1. Constant Immigration
    2. Species differences in physiological parameters
    """

    REStimes = [1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.4,
                3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.6, 5.8, 6.0,
                6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 20.0, 40.0]

    REStimes = [1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.4,
                3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.6, 5.8, 6.0,
                6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 16.0,
                18.0, 20.0, 24.0, 28.0, 32.0, 40.0, 60.0, 70.0, 80.0, 90.0,
                100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 1000.0, 2000.0]

    REStimes = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.4,
                3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.6, 5.8, 6.0,
                6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 20.0, 40.0, 80.0,
                100, 128, 254, 508, 1016]


    #REStimes = [1, 2, 4, 8, 16, 32, 64, 128, 264]

    Qmax = 0.0
    tauAtQmax = 0.0

    AvgGlist = []
    AvgNlist = []   # mean total abundance
    AvgPlist = []   # meanproductivity

    AvgRlist = []  # mean total resources
    AvgQlist = []   # mean cell quota
    AvgMCTlist = [] # mean cell residence time

    if basic == False:
        AvgTOlist = []
        AvgSlist = []
        AvgAblist = []
        AvgHlist = []
        AvgEvarlist = []
        AvgCTlist = []


    tauAtQmaxList = []
    tauAtAvgMuMaxList = []

    ObstauAtAvgMuMaxList = []
    TauList = []

    V, lgp, Kq, PropQ, Pcons, Rcons, uptake_mean, uptake_std, maint, res_pulse, im_pulse, env_noise, spatial_het, disp_mu, disp_std = get_params(system, REStimes)

    iteration = 0
    while iteration < replicates: # number of replicates

        iteration += 1

        REStimes2 = []

        AvgG = []
        AvgN = [] # mean total abundance
        AvgS = []
        AvgP = [] # mean productivity

        AvgR = [] # mean total resources
        AvgQ = [] # mean cell quota

        AvgMCT = [] # mean mean cell residence time
        AvgMaxU = [] # mean maximum uptake

        if basic == False:
            AvgTO = []
            AvgAb = []
            AvgH = []
            AvgEvar = []
            AvgCT = []


        for ri, REStime in enumerate(REStimes):

            if Tau == 'random':
                #REStime = choice(REStimes)
                V, lgp, Kq, PropQ, Pcons, Rcons, uptake_mean, uptake_std, maint, res_pulse, im_pulse, env_noise, spatial_het, disp_mu, disp_std = get_params(system, REStimes)

            print system,': ', iteration,'of', replicates,

            avg_vals = []
            basic = True

            if use_picloud == 'n':
                avg_vals = Hbide(V, REStime, Rcons, Pcons, PropQ, time, uptake_mean, uptake_std, Kq, lgp, maint, res_pulse, im_pulse, env_noise, spatial_het, basic, system, disp_mu, disp_std)

            elif use_picloud == 'y':
                job_id = cloud.call(Hbide, V, REStime, Rcons, Pcons, PropQ, time, uptake_mean, uptake_std, Kq, lgp, maint, res_pulse, im_pulse, env_noise, spatial_het, basic, system, disp_mu, disp_std, _type='m1')
                avg_vals = cloud.result(job_id)


            if len(avg_vals) == 1:
                print avg_vals[0]
                continue


            elif basic == True:

                REStimes2.append(REStime)
                avgG, avgN, avgS, avgR, avgQ, avgMaxU, burnIn, lag, avgP, avgMCT, Hurst, pval = avg_vals

                AvgG.append(avgG)
                AvgN.append(avgN/V)
                AvgP.append(avgP)

                AvgS.append(avgS) # mean richness

                AvgR.append(avgR/V)
                AvgQ.append(avgQ)

                AvgMaxU.append(avgMaxU)
                AvgMCT.append(avgMCT)


            elif basic == False:

                REStimes2.append(REStime)
                avgG, avgN, avgR, avgTO, avgQ, avgS, avgCT, avgEvar, avgH, RADinfo, avgAb, avgMaxU, burnIn, lag, avgP, avgMCT, Hurst, pval = avg_vals

                AvgG.append(avgG)
                AvgN.append(avgN/V)
                AvgP.append(avgP)

                AvgAb.append(avgAb)
                AvgR.append(avgR/V)

                AvgTO.append(avgTO)
                AvgQ.append(avgQ)

                AvgS.append(avgS) # mean richness
                AvgEvar.append(avgEvar) # mean evenness

                AvgH.append(avgH) # mean diversity
                AvgCT.append(avgCT) # mean compositional turnover

                AvgMaxU.append(avgMaxU)
                AvgMCT.append(avgMCT)


            print '[Tau:', REStime, ' MCT:',round(avgMCT,2),
            print ' D:', round(REStime**-1,3),
            print ' g(S):', round(avgP,3), ']',

            print ' N:',round(avgN,3),' S:',round(avgS,1),
            print ' R:', round(avgR,2), '[Burn-in:',burnIn,' H:', round(Hurst,3),
            print ' p:', round(pval,3),' lag:',lag,']'

   	    result1 = isinstance(avgN, float)
            result2 = isinstance(avgR, float)

            if result1 == True and result2 == True:
                if avgQ >= Qmax:
                    Qmax = avgQ
                    tauAtQmax = REStime


        tauAtQmaxList.append(tauAtQmax)
        tauAtAvgMuMaxList.append(uptake_mean**-1)
        MaxAvgMaxUptake = max(AvgMaxU)

        index = AvgMaxU.index(MaxAvgMaxUptake)
        ObstauAtAvgMuMaxList.append(REStimes2[index])

        AvgGlist.append(AvgG)
        AvgNlist.append(AvgN)
        AvgPlist.append(AvgP)

        AvgRlist.append(AvgR)
        AvgMCTlist.append(AvgMCT)

        AvgQlist.append(AvgQ)
        TauList.append(REStimes2)

        if basic == False:

            AvgTOlist.append(avgTO)
            AvgSlist.append(avgS) # mean richness
            AvgAblist.append(avgAb)
            AvgEvarlist.append(avgEvar) # mean evenness
            AvgHlist.append(avgH) # mean diversity
            AvgCTlist.append(avgCT) # mean compositional turnover



    if basic == True:
        return [AvgGlist, AvgNlist, AvgPlist, AvgRlist, AvgMCTlist, AvgQlist, TauList, uptake_mean, tauAtQmaxList]

    elif basic == False:
        return [AvgGlist, AvgNlist, AvgPlist, AvgRlist, AvgMCTlist, AvgQlist, AvgTOlist, AvgSlist, AvgAblist, AvgEvarlist, AvgHlist, AvgCTlist, TauList, uptake_mean, tauAtQmaxList]
