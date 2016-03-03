# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import numpy as np
import sys
#from math import modf
#import decimal
import time

limit = 0.2

def GetIndParam(means):
    vals = []

    if isinstance(means, float) or isinstance(means, int):
        std = means/1000.0
        vals = np.random.normal(means, std)
        if vals < 0.00001:
            vals = 0.00001

    else:
        for val in means:
            std = val/1000.0
            i = np.random.normal(val, std)
            if i < 0.00001:
                i = 0.00001
            vals.append(i)

    return vals



def GetRAD(vector):
    RAD = []
    unique = list(set(vector))

    for val in unique:
        RAD.append(vector.count(val)) # the abundance of each Sp_

    return RAD, unique # the rad and the specieslist


def get_color(ID, colorD): # FUNCTION TO ASSIGN COLORS TO Sp_

    r1 = lambda: randint(0,255)
    r2 = lambda: randint(0,255)
    r3 = lambda: randint(0,255)

    color = '#%02X%02X%02X' % (r1(),r2(),r3())
    colorD[ID] = color

    return colorD



def NewTracers(motion, IDs, Xs, Ys, t_In, w, h, u0, ct):

    if ct == 1:
        for i in range(200):
            IDs.append(0)
            t_In.append(0)
            Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
            Xs.append(float(np.random.uniform(0.1*w, 0.101*w)))
    else:                                    
        x = np.random.binomial(1, u0)
        if x == 1:
            IDs.append(0)
            t_In.append(0)
            Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
            Xs.append(float(np.random.uniform(0.1*w, 0.101*w)))

    return [IDs, t_In, Xs, Ys]



def ResIn(motion, Type, Vals, Xs, Ys, ID, IDs, t_In, numr, rmax, nN, nP, nC, w, h, u0):

    for r in range(numr):
        x = np.random.binomial(1, u0/2)

        if x == 1:
            rval = int(np.random.random_integers(1, rmax, 1))
            nr = choice(['N', 'P', 'C'])

            if nr == 'N':
                rtype = int(np.random.random_integers(0, nN-1, 1))
                rtype = 'N'+str(rtype)

            if nr == 'P':
                rtype = int(np.random.random_integers(0, nP-1, 1))
                rtype = 'P'+str(rtype)

            if nr == 'C':
                rtype = int(np.random.random_integers(0, nC-1, 1))
                rtype = 'C'+str(rtype)

            Vals.append(rval)
            IDs.append(ID)
            Type.append(rtype)
            t_In.append(0)
            ID += 1

            if motion == 'random_walk':
                Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
                Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))

            else:
                Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
                Xs.append(float(np.random.uniform(0.1*w, 0.2*w)))


    return [Type, Vals, Xs, Ys, IDs, ID, t_In]



def immigration(d_max, g_max, m_max, motion, seed, ip, Sp, Xs, Ys, w, h, MD,
        EnvD, envGs, GD, DispD, colorD, IDs, ID, t_In, Qs, N_RD, P_RD, C_RD,
        nN, nP, nC, u0, alpha, GList, MList, NList, PList, CList, DList):

    if u0 > 1.0:
        u0 = 1.0

    for m in range(seed):
        x = 0

        if seed > 1:
            x = 1

        else:
            x = np.random.binomial(1, u0*ip)

        if x == 1:

            prop = str(float(np.random.logseries(alpha, 1)))

            Sp.append(prop)

            if motion == 'random_walk':
                Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
                Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))

            else:
                Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
                Xs.append(float(np.random.uniform(0.1*w, 0.2*w)))

            IDs.append(ID)
            t_In.append(0)
            ID += 1
            Qn = float(np.random.uniform(0.01, 0.1))
            Qp = float(np.random.uniform(0.01, 0.1))
            Qc = float(np.random.uniform(0.01, 0.1))

            Qs.append([Qn, Qp, Qc])

            if prop not in colorD:
                # speciescolor
                colorD = get_color(prop, colorD)

                # species growth rate
                GD[prop] = np.random.uniform(g_max/10, g_max)

                # species maintenance
                MD[prop] = np.random.uniform(m_max/10, m_max)

                # species active dispersal rate
                DispD[prop] = np.random.uniform(d_max/10, d_max)

                # species environmental gradient optima
                glist = []
                for g in envGs:
                    x = np.random.uniform(0.0, w)
                    y = np.random.uniform(0.0, h)
                    glist.append([x,y])
                EnvD[prop] = glist

                # species Nitrogen use efficiency
                N_RD[prop] = np.random.uniform(0.01, 1.0, nN)

                # species Phosphorus use efficiency
                P_RD[prop] = np.random.uniform(0.01, 1.0, nP)

                # species Carbon use efficiency
                C_RD[prop] = np.random.uniform(0.01, 1.0, nC)

            means = GD[prop]
            i = GetIndParam(means)
            GList.append(i)

            means = MD[prop]
            i = GetIndParam(means)
            MList.append(i)

            means = N_RD[prop]
            i = GetIndParam(means)
            NList.append(i)

            means = P_RD[prop]
            i = GetIndParam(means)
            PList.append(i)

            means = C_RD[prop]
            i = GetIndParam(means)
            CList.append(i)

            means = DispD[prop]
            i = GetIndParam(means)
            DList.append(i)

    return [Sp, Xs, Ys, MD, EnvD, GD, DispD, colorD, IDs, ID, t_In, Qs, N_RD,
            P_RD, C_RD, GList, MList, NList, PList, CList, DList]



def fluid_movement(TypeOf, List, t_In, xAge, Xs, Ys, ux, uy, w, h, u0):

    Type, IDs, ID, Vals = [], [], int(), []

    if TypeOf == 'resource':
        Type, IDs, ID, Vals = List
    elif TypeOf == 'individual':
        Type, IDs, ID, Vals, DispD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = List
    else:
        IDs = List

    if Xs == []:
        if TypeOf == 'tracer':
            return [IDs, Xs, Ys, xAge, t_In]
        elif TypeOf == 'individual':
            return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals, GrowthList,
                    MaintList, N_RList, P_RList, C_RList, DispList]
        elif TypeOf == 'resource':
            return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals]

    ux = np.reshape(ux, (w*h)) # ux is the macroscopic x velocity
    uy = np.reshape(uy, (w*h)) # uy is the macroscopic y velocity

    # dispersal inside the system

    n = len(Xs)
    for j in range(n):

        i = randint(0, len(Xs)-1)
        X = int(round(Xs[i]))
        Y = int(round(Ys[i]))

        index =  int(round(X + Y * w))

        if index > len(ux) - 1:
            index = len(ux) - 1
        if index > len(uy) - 1:
            index = len(uy) - 1

        k = 0
        if TypeOf == 'individual':
            k = np.random.binomial(1, DispD[Type[i]])

        if k == 0:
            Xs[i] += ux[index]
            Ys[i] += uy[index]

        y = Ys[i]

        if 0.0 > y:
            Ys[i] = 0.0
        elif y >= h:
            Ys[i] = h - 0.0

        t_In[i] += 1
        if Xs[i] <= 0:
            Xs[i] = 0.0

        if Xs[i] >= w - limit:

            xAge.append(t_In[i])
            Xs.pop(i)
            Ys.pop(i)
            t_In.pop(i)
            IDs.pop(i)

            if TypeOf == 'resource' or TypeOf == 'individual':
                Type.pop(i)
                Vals.pop(i)

            if TypeOf == 'individual':
                GrowthList.pop(i)
                MaintList.pop(i)
                N_RList.pop(i)
                P_RList.pop(i)
                C_RList.pop(i)
                DispList.pop(i)

    ux = np.reshape(ux, (h, w))
    uy = np.reshape(uy, (h, w))

    if TypeOf == 'tracer':
        return [IDs, Xs, Ys, xAge, t_In]
    elif TypeOf == 'individual':
        return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals, GrowthList, MaintList,
            N_RList, P_RList, C_RList, DispList]
    elif TypeOf == 'resource':
        return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals]




def nonfluid_movement(TypeOf, motion, List, t_In, xAge, Xs, Ys, ux, uy, w, h, u0):

    Type, IDs, ID, Vals = [], [], int(), []

    if TypeOf == 'resource':
        Type, IDs, ID, Vals = List
    elif TypeOf == 'individual':
        Type, IDs, ID, Vals, DispD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = List
    else:
        IDs = List

    if Xs == []:
        if TypeOf == 'tracer':
            return [IDs, Xs, Ys, xAge, t_In]
        elif TypeOf == 'individual':
            return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals, GrowthList,
                MaintList, N_RList, P_RList, C_RList, DispList]
        elif TypeOf == 'resource':
            return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals]


    limit, distance, direction, pop = 0.1, 0, 0, 'no'
    x, y = 0, 0

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        # get distance
        if TypeOf == 'individual':
            distance = np.random.uniform(0, u0)
        else:
            distance = np.random.uniform(0, u0)

        x, y = Xs[i], Ys[i]

        # go up or down
        #if TypeOf != 'resource':
        if motion == 'unidirectional':
            direction = 1
        if motion == 'random_walk':
            direction = choice([-1, 1])

        y = y + (direction * distance)

        # get distance
        if TypeOf == 'individual':
            distance = np.random.uniform(0, u0)
        else:
            distance = np.random.uniform(0, u0)

        # go forward or backward
        #if TypeOf != 'resource':
        if motion == 'unidirectional':
            direction = 1
        if motion == 'random_walk':
            direction = choice([-1, 1])

        x = x + (direction * distance)

        if x > w - limit or x < limit:
            pop = 'yes'
        elif y > h - limit or y < limit:
            pop = 'yes'


        if pop == 'no':
            Xs[i], Ys[i] = x, y
            t_In[i] = t_In[i]+1

        elif pop == 'yes':
            xAge.append(t_In[i])
            Xs.pop(i)
            Ys.pop(i)
            t_In.pop(i)
            IDs.pop(i)

            if TypeOf == 'resource' or TypeOf == 'individual':
                Type.pop(i)
                Vals.pop(i)

            if TypeOf == 'individual':
                GrowthList.pop(i)
                MaintList.pop(i)
                N_RList.pop(i)
                P_RList.pop(i)
                C_RList.pop(i)
                DispList.pop(i)

    if TypeOf == 'tracer':
        return [IDs, Xs, Ys, xAge, t_In]
    elif TypeOf == 'individual':
        return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals, GrowthList, MaintList,
                N_RList, P_RList, C_RList, DispList]
    elif TypeOf == 'resource':
        return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals]




def predation(P_IDs, P_ID, P_Xs, P_Ys, P_t_In, I_xAge, Sp_IDs, Qs, I_IDs,
        I_ID, I_t_In, I_Xs, I_Ys, w, h):

    """ This function is currently under development """

    I_Boxes, P_Boxes = [], []
    if not len(P_IDs):
        List = [P_IDs, P_ID, P_t_In, P_Xs, P_Ys, Sp_IDs, Qs]
        List += [I_IDs, I_ID, I_t_In, I_Xs, I_Ys]
        return List

    I_Boxes, P_Boxes = [[list([]) for _ in xrange(w*h)]]*2

    index = 0
    for i, val in enumerate(I_IDs):
        rX = int(round(I_Xs[i]))
        rY = int(round(I_Ys[i]))

        index = int(round(rX + (rY * w)))

        if index > len(I_Boxes) - 1:
            index = len(I_Boxes) - 1
        elif index < 0:
            index = 0
        I_Boxes[index].append(val)

    index = 0
    for i, val in enumerate(P_IDs):
        rX = int(round(P_Xs[i]))
        rY = int(round(P_Ys[i]))

        index = int(round(rX + (rY * w)))

        if index > len(P_Boxes) - 1:
            index = len(P_Boxes) - 1
        elif index < 0:
            index = 0
        P_Boxes[index].append(val)


    for i, P_Box in enumerate(P_Boxes):
        if not len(P_Box): continue

        I_Box = I_Boxes[i]

        for P_ in P_Box:
            if not len(I_Box): break

            # The predator
            P_Box.pop(0)

            # The prey
            I_ID = choice(I_Box)
            index = I_Box.index(I_ID)
            I_Box.pop(index)

            j = I_IDs.index(I_ID)
            Qs.pop(j)
            I_xAge.append(I_t_In[j])
            I_t_In.pop(j)
            Sp_IDs.pop(j)
            I_IDs.pop(j)
            I_Xs.pop(j)
            I_Ys.pop(j)
            #GrowthList.pop(i)
            #MaintList.pop(i)
            #N_RList.pop(i)
            #P_RList.pop(i)
            #C_RList.pop(i)
            #DispList.pop(i)

    List = [P_IDs, P_ID, P_t_In, P_Xs, P_Ys, Sp_IDs, Qs, I_IDs]
    List += [I_ID, I_t_In, I_Xs, I_Ys]
    return List





def maintenance(Sp_IDs, Xs, Ys, xAge, colorD, MD, EnvD, IDs, t_In, Qs, GrowthList,
        MaintList, N_RList, P_RList, C_RList, DispList):

    if Sp_IDs == []:
        return [Sp_IDs, Xs, Ys, xAge, IDs, t_In, Qs, GrowthList, MaintList, N_RList,
                P_RList, C_RList, DispList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        val = Qs[i]
        val[0] -= MaintList[i]  # maintanence influenced by species id
        val[1] -= MaintList[i]
        val[2] -= MaintList[i]

        if min(val) <= 0.0001:   # starved

            Qs.pop(i)
            xAge.append(t_In[i])
            t_In.pop(i)
            Sp_IDs.pop(i)
            IDs.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            GrowthList.pop(i)
            MaintList.pop(i)
            N_RList.pop(i)
            P_RList.pop(i)
            C_RList.pop(i)
            DispList.pop(i)

        else: Qs[i] = val

    return [Sp_IDs, Xs, Ys, xAge, IDs, t_In, Qs, GrowthList, MaintList, N_RList,
            P_RList, C_RList, DispList]



def decimate(Sp_IDs, Xs, Ys, xAge, colorD, MD, EnvD, IDs, t_In, Qs, GrowthList,
            MaintList, N_RList, P_RList, C_RList, DispList):

    if Sp_IDs == []:
        return [Sp_IDs, Xs, Ys, xAge, IDs, t_In, Qs, GrowthList, MaintList,
        N_RList, P_RList, C_RList, DispList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        d = np.random.binomial(1, 0.1)

        if d == 1:   # remvoe

            Qs.pop(i)
            xAge.append(t_In[i])
            t_In.pop(i)
            Sp_IDs.pop(i)
            IDs.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            GrowthList.pop(i)
            MaintList.pop(i)
            N_RList.pop(i)
            P_RList.pop(i)
            C_RList.pop(i)
            DispList.pop(i)

    return [Sp_IDs, Xs, Ys, xAge, IDs, t_In, Qs, GrowthList, MaintList, N_RList,
            P_RList, C_RList, DispList]




def consume(R_Types, R_Vals, R_IDs, R_ID, R_Xs, R_Ys, R_t_In, R_xAge, Sp_IDs,
        Qs, I_IDs, I_ID, I_t_In, I_Xs, I_Ys, w, h, GD, N_RD, P_RD, C_RD, DispD,
        GrowthList, MaintList, N_RList, P_RList, C_RList, DispList):

    if not len(R_Types) or not len(Sp_IDs):
        List = [R_Types, R_Vals, R_IDs, R_ID, R_t_In, R_xAge, R_Xs]
        List += [R_Ys, Sp_IDs, Qs, I_IDs, I_ID, I_t_In]
        List += [I_Xs, I_Ys, GrowthList, MaintList, N_RList,
                P_RList, C_RList, DispList]
        return List

    I_Boxes = [list([]) for _ in xrange(w*h)]
    R_Boxes = [list([]) for _ in xrange(w*h)]

    index = 0
    for i, val in enumerate(I_IDs):
        rX = int(round(I_Xs[i]))
        rY = int(round(I_Ys[i]))

        index = int(round(rX + (rY * w)))

        if index > len(I_Boxes) - 1:
            index = len(I_Boxes) - 1
        elif index < 0:
            index = 0

        I_Boxes[index].append(val)

    index = 0
    for i, val in enumerate(R_IDs):

        rX = int(round(R_Xs[i]))
        rY = int(round(R_Ys[i]))
        index = int(round(rX + (rY * w)))

        if index > len(R_Boxes) - 1:
            index = len(R_Boxes) - 1
        elif index < 0:
            index = 0

        R_Boxes[index].append(val)


    for i, box in enumerate(I_Boxes):
        if not len(box): continue

        R_Box = R_Boxes[i]

        for ind in box: # The individuals
            if not len(R_Box): break

            R_ID = choice(R_Box)
            boxI_ex = R_Box.index(R_ID)

            # The food
            j = R_IDs.index(R_ID)
            R_val = R_Vals[j]
            R_type = R_Types[j]

            rtype = list(R_type)
            R = rtype.pop(0)
            rnum = int(''.join(rtype))

            # The Individual
            ID = I_IDs.index(ind)
            # The individual's cell quota

            Q = Qs[ID]
            QN = Q[0]
            QP = Q[1]
            QC = Q[2]

            # the species
            sp = Sp_IDs[ID]
            mu = GD[sp]

            Q = 0.0
            efficiency = 0.0

            if R == 'N':
                efficiency = N_RList[ID][rnum]
                Q = QN

            if R == 'P':
                efficiency = P_RList[ID][rnum]
                Q = QP

            if R == 'C':
                efficiency = C_RList[ID][rnum]
                Q = QC

            mu = mu * efficiency

            if R_val > (mu * Q): # Increase cell quota
                R_val = R_val - (mu * Q)
                Q += (mu * Q)

            else:
                Q += R_val
                R_val = 0.0

            if Q > 1.0:
                R_val = Q - 1.0
                Q = 1.0
                R_Vals[j] = R_val


            if R_val <= 0.0:
                R_Box.pop(boxI_ex)
                R_Vals.pop(j)
                R_xAge.append(R_t_In[j])
                R_t_In.pop(j)
                R_Types.pop(j)
                R_IDs.pop(j)
                R_Xs.pop(j)
                R_Ys.pop(j)


            if Q < 0.0:
                print Q, QN, QP, QC
                sys.exit()

            if R == 'N':
                Qs[ID] = [Q, QP, QC]
            if R == 'P':
                Qs[ID] = [QN, Q, QC]
            if R == 'C':
                Qs[ID] = [QN, QP, Q]


    return [R_Types, R_Vals, R_IDs, R_ID, R_t_In, R_xAge, R_Xs, R_Ys, Sp_IDs,
            Qs, I_IDs, I_ID, I_t_In, I_Xs, I_Ys, GrowthList, MaintList, N_RList,
            P_RList, C_RList, DispList]



def reproduce(repro, spec, Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, w, h, GD, DispD,
        colorD, N_RD, P_RD, C_RD, MD, EnvD, envGs, nN, nP, nC, GList, MList,
        NList, PList, CList, DList):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, GD, DispD, GList, MList,
                NList, PList, CList, DList]

    if repro == 'fission':


        n = len(IDs)
        for j in range(n):

            i = randint(0, len(IDs)-1)

            Q = Qs[i]
            pq = float(np.mean(Q))
            p = np.random.binomial(1, pq)

            if p == 1: # individual is large enough to reproduce

                spID = Sp_IDs[i]
                X = Xs[i]
            	Y = Ys[i]

                pg = []
                sp_opts = EnvD[spID]

                for g, opt in enumerate(sp_opts):

                    x, y = envGs[g]
                    pg.append(1 - (abs(X - x)/max([X,x])))
                    pg.append(1 - (abs(Y - y)/max([Y,y])))


                if np.mean(pg) > 1 or np.mean(pg) < 0:
                    print pg
                    sys.exit()

                p = np.mean(pg)
                p = np.random.binomial(1, p)
                if p == 1: # the environment is suitable for reproduction

                    QN = Q[0]
                    QP = Q[1]
                    QC = Q[2]

                    Qs[i] = [QN/2.0, QP/2.0, QC/2.0]
                    Qs.append([QN/2.0, QP/2.0, QC/2.0])

                    ID += 1
                    IDs.append(ID)
                    t_In.append(t_In[i])

                    p = np.random.binomial(1, spec)
                    if p == 1:

                        # speciate
                        t = time.clock()
                        spID_new = spID +' '+ str(t)

                        # new speciescolor
                        colorD = get_color(spID_new, colorD)

                        # new species growth rate
                        p = np.random.binomial(1, 0.25)
                        if p == 1: GD[spID_new] = np.random.uniform(0.5, 1.0)
                        else: GD[spID_new] = GD[spID]

                        # new speciesmaintenance
                        p = np.random.binomial(1, 0.25)
                        if p == 1: MD[spID_new] = np.random.uniform(0.01, 0.1)
                        else: MD[spID_new] = MD[spID]

                        # species environmental gradient optima
                        glist = []
                        for j, g in enumerate(envGs):
                            p = np.random.binomial(1, 0.25)
                            if p == 1:
                                x = np.random.uniform(0.0, w)
                                y = np.random.uniform(0.0, h)
                            else:
                                x = EnvD[spID][j][0]
                                y = EnvD[spID][j][1]

                            glist.append([x,y])
                            EnvD[spID_new] = glist

                        # new speciesactive dispersal rate
                        p = np.random.binomial(1, 0.25)
                        if p == 1: DispD[spID_new] = np.random.uniform(0.0, 0.1)
                        else: DispD[spID_new] = DispD[spID]

                        # new speciesresource use efficiencies
                        # Nitrogen
                        p = np.random.binomial(1, 0.25)
                        if p == 1: N_RD[spID_new] = np.random.uniform(0.01, 1.0, nN)
                        else: N_RD[spID_new] = N_RD[spID]

                        # Phosphorus
                        p = np.random.binomial(1, 0.25)
                        if p == 1: P_RD[spID_new] = np.random.uniform(0.01, 1.0, nP)
                        else: P_RD[spID_new] = P_RD[spID]

                        # Carbon
                        p = np.random.binomial(1, 0.25)
                        if p == 1: C_RD[spID_new] = np.random.uniform(0.01, 1.0, nC)
                        else: C_RD[spID_new] = C_RD[spID]

                        spID = spID_new

                    means = GD[spID]
                    i = GetIndParam(means)
                    GList.append(i)

                    means = MD[spID]
                    i = GetIndParam(means)
                    MList.append(i)

                    means = N_RD[spID]
                    i = GetIndParam(means)
                    NList.append(i)

                    means = P_RD[spID]
                    i = GetIndParam(means)
                    PList.append(i)

                    means = C_RD[spID]
                    i = GetIndParam(means)
                    CList.append(i)

                    means = DispD[spID]
                    i = GetIndParam(means)
                    DList.append(i)

                    Sp_IDs.append(spID)

                    newX = float(np.random.uniform(X-0.1, X, 1))
                    if limit > newX: newX = 0
                    if newX > w - limit: newX = w - limit
                    Xs.append(newX)

                    newY = float(np.random.uniform(Y-0.1, Y+0.1, 1))
                    if limit > newY: newY = 0
                    elif newY > h: newY = h - limit
                    Ys.append(newY)


    elif repro == 'sexual':

        I_Boxes = []
        spBoxes = []

        I_Boxes = [list([]) for _ in xrange(w*h)]
        spBoxes = [list([]) for _ in xrange(w*h)]

        index = 0
        for i, I_ID in enumerate(IDs):

            spID = Sp_IDs[i]
            rX = int(round(Xs[i]))
            rY = int(round(Ys[i]))

            index = int(round(rX + (rY * w)))

            if index > len(I_Boxes) - 1:
                index = len(I_Boxes) - 1
            elif index < 0:
                index = 0

            I_Boxes[index].append(I_ID)
            spBoxes[index].append(spID)

        for i, box in enumerate(I_Boxes):
            if len(box) < 2: continue

            spbox = spBoxes[i]
            while len(box) > 1:
                # choose an individual at random
                i1 = choice(range(len(box)))
                ind1 = box.pop(i1)
                # remove the speciesID
                sp = spbox.pop(i1)

                index1 = IDs.index(ind1)
                q = Qs[index1]
                q1 = np.mean(q)

                p1 = np.random.binomial(1, q1)
                # individual not large enough to reproduce
                if p1 == 0: continue

                # Find another of the same Sp_
                if spbox.count(sp) > 1:

                    QN = q[0]
                    QP = q[1]
                    QC = q[2]

                    # choose an individual of the same Sp_
                    i2 = spbox.index(sp)
                    # remove the speciesID
                    spbox.pop(i2)
                    # remove the individual
                    box.pop(i2)

                    index2 = IDs.index(i2)
                    q2 = Qs[index2]

                    p2 = np.random.binomial(1, q2)
                    # individual not large enough to reproduce
                    if p2 == 0: continue

                    Qs[i2] = q2/2.0
                    Qs.append(q2/2.0)

                    X = Xs[i2]
                    Y = Ys[i2]

                    newX = float(np.random.uniform(X-0.5, X+0.5, 1))
                    if limit > newX: newX = 0
                    if newX > w - limit: newX = w - limit

                    newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                    if limit > newY: newY = 0
                    elif newY > h: newY = h - limit

                    ID += 1
                    IDs.append(ID)
                    t_In.append(0)
                    Sp_IDs.append(sp)

    return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, GD, DispD, GList, MList,
                NList, PList, CList, DList]




def search(repro, spec, Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys,  w, h, GD, DispD,
        colorD, N_RD, P_RD, C_RD, MD, EnvD, envGs, nN, nP, nC, GList, MList,
        NList, PList, CList, DList):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, GD, DispD, GList,
                    MList, NList, PList, CList, DList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        spID = Sp_IDs[i]
        X = Xs[i]
        Y = Ys[i]

        sp_opts = EnvD[spID]

        for g, opt in enumerate(sp_opts):
            x, y = opt
            dist = 0.5*DispD[spID]

            if x > X:
                X += dist

            elif x < X:
                X -= dist

            if y > Y:
                Y += dist

            elif y < Y:
                Y -= dist

            if X > w: X = w
            elif X < 0: X = 0
            if Y > h: Y = h
            elif Y < 0: Y = 0

            Xs[i] = X
            Ys[i] = Y


    return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, GD, DispD, GList, MList,
                NList, PList, CList, DList]
