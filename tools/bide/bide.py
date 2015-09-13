# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import numpy as np
import sys
#from math import modf
#import decimal
import time

limit = 0.2

def coord(d):
    return float(np.random.uniform(0.1*d, 0.9*d))


def GetIndParam(means):
    vals = []

    if isinstance(means, float) or isinstance(means, int):
        std = means/10.0
        vals = np.random.normal(means, std)
        if vals < 0.00001:
            vals = 0.00001

    else:
        for val in means:
            std = val/10.0
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



def NewTracers(motion, IDs, Xs, Ys, Zs, t_In, w, h, l, u0, D):

    x = np.random.binomial(1, u0)
    if x == 1:
        IDs.append(0)
        t_In.append(0)
        if motion == 'random_walk':
            Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
            Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))
            Zs.append(float(np.random.uniform(0.1*l, 0.9*l)))

        else:
            Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
            Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))
            Zs.append(float(np.random.uniform(0.1*l, 0.9*l)))

    return [IDs, t_In, Xs, Ys, Zs]



def ResIn(motion, Type, Vals, Xs, Ys, Zs, ID, IDs, t_In, numr,
        rmax, nN, nP, nC, w, h, l, u0, D):

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
                Zs.append(float(np.random.uniform(0.1*l, 0.9*l)))

            else:
                Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
                Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))
                Zs.append(float(np.random.uniform(0.1*l, 0.9*l)))


    return [Type, Vals, Xs, Ys, Zs, IDs, ID, t_In]



def immigration(d_max, g_max, m_max, motion, numin, Sp, Xs, Ys, Zs, w, h, l, MD, EnvD, envGs,
        GD, DispD, colorD, IDs, ID, t_In, Qs, N_RD, P_RD, C_RD, nN, nP, nC, u0, alpha, D, GList, MList, NList, PList, CList, DList):

    for m in range(numin):
        x = np.random.binomial(1, u0/2)

        if x == 1:
            prop = str(float(np.random.logseries(0.999, 1)))

            Sp.append(prop)

            if motion == 'random_walk':
                Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
                Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))
                Zs.append(float(np.random.uniform(0.1*l, 0.9*l)))

            else:
                Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
                Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))
                Zs.append(float(np.random.uniform(0.1*l, 0.9*l)))

            IDs.append(ID)
            t_In.append(0)
            ID += 1
            Qn = float(np.random.uniform(0.05, 0.5))
            Qp = float(np.random.uniform(0.05, 0.5))
            Qc = float(np.random.uniform(0.05, 0.5))

            Qs.append([Qn, Qp, Qc])

            if prop not in colorD:
                # speciescolor
                colorD = get_color(prop, colorD)

                # species growth rate
                GD[prop] = np.random.uniform(0.1, g_max)

                # species maintenance
                MD[prop] = np.random.uniform(0.001, m_max)

                # species active dispersal rate
                DispD[prop] = np.random.uniform(0.01, d_max)

                # species environmental gradient optima
                EnvD[prop] = np.random.uniform(0.0, 1.0, len(envGs))

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

    return [Sp, Xs, Ys, Zs, MD, EnvD, GD, DispD, colorD, IDs, ID, t_In, Qs, N_RD, P_RD, C_RD, GList, MList, NList, PList, CList, DList]



def fluid_movement(TypeOf, List, t_In, xAge, Xs, Ys, Zs, ux, uy, w, h, u0):

    Type, IDs, ID, Vals = [], [], int(), []

    if TypeOf == 'resource':
        Type, IDs, ID, Vals = List
    elif TypeOf == 'individual':
        Type, IDs, ID, Vals, DispD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = List
    else:
        IDs = List

    if Xs == []:
        if TypeOf == 'tracer':
            return [IDs, Xs, Ys, Zs, xAge, t_In]
        elif TypeOf == 'individual':
            return [Type, Xs, Ys, Zs, xAge, IDs, ID, t_In, Vals, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]
        elif TypeOf == 'resource':
            return [Type, Xs, Ys, Zs, xAge, IDs, ID, t_In, Vals]

    ux = np.reshape(ux, (w*h)) # ux is the macroscopic x velocity
    uy = np.reshape(uy, (w*h)) # uy is the macroscopic y velocity

    # dispersal inside the system
    for i, val in enumerate(Xs):

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

        if 0 > y:
            Ys[i] = 0
        elif y > h:
            Ys[i] = h

        t_In[i] += 1
        if Xs[i] >= w - limit:

            xAge.append(t_In[i])
            Xs.pop(i)
            Ys.pop(i)
            Zs.pop(i)
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
        return [IDs, Xs, Ys, Zs, xAge, t_In]
    elif TypeOf == 'individual':
        return [Type, Xs, Ys, Zs, xAge, IDs, ID, t_In, Vals, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]
    elif TypeOf == 'resource':
        return [Type, Xs, Ys, Zs, xAge, IDs, ID, t_In, Vals]





def predation(P_IDs, P_ID, P_Xs, P_Ys, P_Zs, P_t_In, I_xAge, Sp_IDs,
        Qs, I_IDs, I_ID, I_t_In, I_Xs, I_Ys, I_Zs, w, h, l, D):

    I_Boxes, P_Boxes = [], []
    if not len(P_IDs):
        List = [P_IDs, P_ID, P_t_In, P_Xs, P_Ys, P_Zs, Sp_IDs, Qs]
        List += [I_IDs, I_ID, I_t_In, I_Xs, I_Ys, I_Zs]
        return List

    if D == 2:
        I_Boxes, P_Boxes = [[list([]) for _ in xrange(w*h)]]*2
    elif D == 3:
        I_Boxes, P_Boxes = [[list([]) for _ in xrange(w*h*l)]]*2

    index = 0
    for i, val in enumerate(I_IDs):
        rX = int(round(I_Xs[i]))
        rY = int(round(I_Ys[i]))
        rZ = int(round(I_Zs[i]))

        if D == 2:
            index = int(round(rX + (rY * w)))
        elif D == 3:
            index = int(round((rY * l * w) + (rX * l) + rZ))

        if index > len(I_Boxes) - 1:
            index = len(I_Boxes) - 1
        elif index < 0:
            index = 0
        I_Boxes[index].append(val)

    index = 0
    for i, val in enumerate(P_IDs):
        rX = int(round(P_Xs[i]))
        rY = int(round(P_Ys[i]))
        rZ = int(round(P_Zs[i]))

        if D == 2:
            index = int(round(rX + (rY * w)))
        elif D == 3:
            index = int(round((rY * l * w) + (rX * l) + rZ))

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
            I_Zs.pop(j)
            #GrowthList.pop(i)
            #MaintList.pop(i)
            #N_RList.pop(i)
            #P_RList.pop(i)
            #C_RList.pop(i)
            #DispList.pop(i)

    List = [P_IDs, P_ID, P_t_In, P_Xs, P_Ys, P_Zs, Sp_IDs, Qs, I_IDs]
    List += [I_ID, I_t_In, I_Xs, I_Ys, I_Zs]
    return List





def maintenance(Sp_IDs, Xs, Ys, Zs, xAge, colorD, MD, EnvD, IDs, t_In, Qs, D,  GrowthList, MaintList, N_RList, P_RList, C_RList, DispList):

    if Sp_IDs == []:
        return [Sp_IDs, Xs, Ys, Zs, xAge, IDs, t_In, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]

    for i, val in enumerate(Qs):

        val[0] -= MD[Sp_IDs[i]]  # maintanence influenced by species id
        val[1] -= MD[Sp_IDs[i]]
        val[2] -= MD[Sp_IDs[i]]

        if min(val) <= 0.01:   # starved

            Qs.pop(i)
            xAge.append(t_In[i])
            t_In.pop(i)
            Sp_IDs.pop(i)
            IDs.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            Zs.pop(i)
            GrowthList.pop(i)
            MaintList.pop(i)
            N_RList.pop(i)
            P_RList.pop(i)
            C_RList.pop(i)
            DispList.pop(i)

        else: Qs[i] = val

    return [Sp_IDs, Xs, Ys, Zs, xAge, IDs, t_In, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]



def decimate(Sp_IDs, Xs, Ys, Zs, xAge, colorD, MD, EnvD, IDs, t_In, Qs, D,  GrowthList, MaintList, N_RList, P_RList, C_RList, DispList):

    if Sp_IDs == []:
        return [Sp_IDs, Xs, Ys, Zs, xAge, IDs, t_In, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]

    for i, val in enumerate(Qs):

        d = np.random.binomial(1, 0.1)

        if d == 1:   # starved

            Qs.pop(i)
            xAge.append(t_In[i])
            t_In.pop(i)
            Sp_IDs.pop(i)
            IDs.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            Zs.pop(i)
            GrowthList.pop(i)
            MaintList.pop(i)
            N_RList.pop(i)
            P_RList.pop(i)
            C_RList.pop(i)
            DispList.pop(i)

        else: Qs[i] = val

    return [Sp_IDs, Xs, Ys, Zs, xAge, IDs, t_In, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]




def consume(R_Types, R_Vals, R_IDs, R_ID, R_Xs, R_Ys, R_Zs, R_t_In,
        R_xAge, Sp_IDs, Qs, I_IDs, I_ID, I_t_In, I_Xs, I_Ys, I_Zs,
        w, h, l, GD, N_RD, P_RD, C_RD, DispD, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList):

    if not len(R_Types) or not len(Sp_IDs):
        List = [R_Types, R_Vals, R_IDs, R_ID, R_t_In, R_xAge, R_Xs]
        List += [R_Ys, R_Zs, Sp_IDs, Qs, I_IDs, I_ID, I_t_In]
        List += [I_Xs, I_Ys, I_Zs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]
        return List

    if D == 2:
        I_Boxes = [list([]) for _ in xrange(w*h)]
        R_Boxes = [list([]) for _ in xrange(w*h)]

    elif D == 3:
        I_Boxes = [list([]) for _ in xrange(w*h*l)]
        R_Boxes = [list([]) for _ in xrange(w*h*l)]

    index = 0
    for i, val in enumerate(I_IDs):
        rX = int(round(I_Xs[i]))
        rY = int(round(I_Ys[i]))
        rZ = int(round(I_Zs[i]))

        if D == 2:
            index = int(round(rX + (rY * w)))
        elif D == 3:
            index = int(round((rY * l * w) + (rX * l) + rZ))

        if index > len(I_Boxes) - 1:
            index = len(I_Boxes) - 1
        elif index < 0:
            index = 0

        I_Boxes[index].append(val)

    index = 0
    for i, val in enumerate(R_IDs):

        rX = int(round(R_Xs[i]))
        rY = int(round(R_Ys[i]))
        rZ = int(round(R_Zs[i]))

        if D == 2:
            index = int(round(rX + (rY * w)))
        elif D == 3:
            index = int(round((rY * l * w) + (rX * l) + rZ))

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
                efficiency = N_RD[sp][rnum]
                Q = QN

            if R == 'P':
                efficiency = P_RD[sp][rnum]
                Q = QP

            if R == 'C':
                efficiency = C_RD[sp][rnum]
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
                R_Zs.pop(j)

            if Q < 0.0:
                print Q, QN, QP, QC
                sys.exit()

            if R == 'N':
                Qs[ID] = [Q, QP, QC]
            if R == 'P':
                Qs[ID] = [QN, Q, QC]
            if R == 'C':
                Qs[ID] = [QN, QP, Q]


    List = [R_Types, R_Vals, R_IDs, R_ID, R_t_In, R_xAge, R_Xs, R_Ys]
    List += [R_Zs,Sp_IDs, Qs, I_IDs, I_ID, I_t_In, I_Xs, I_Ys, I_Zs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]
    return List



def reproduce(repro, spec, Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, Zs, w,
        h, l, GD, DispD, colorD, N_RD, P_RD, C_RD, MD, EnvD, envGs, D, nN, nP, nC, GList, MList, NList, PList, CList, DList):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, Zs, GD, DispD, GList, MList, NList, PList, CList, DList]

    if repro == 'fission':
        for i, Q in enumerate(Qs):

            pq = float(np.min(Q))
            p = np.random.binomial(1, pq)

            if p == 1: # individual is large enough to reproduce

                spID = Sp_IDs[i]
                X = Xs[i]
            	Y = Ys[i]
                Z = Zs[i]

                pg = []
                sp_opts = EnvD[spID]

                for g, opt in enumerate(sp_opts):

                    x, y = envGs[g]
                    pg.append(1 - (abs(X - x)/max([X,x])))
                    pg.append(1 - (abs(Y - y)/max([Y,y])))


                if np.mean(pg) > 1:
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


                    if spec == 'yes':
                        p = np.random.binomial(1, 0.05)
                        if p == 1:

                            # speciate
                            t = time.clock()
                            spID = spID +' '+ str(t)

                            # new speciescolor
                            colorD = get_color(spID, colorD)

                            # new speciesgrowth rate
                            p = np.random.binomial(0.25, 1)
                            GD[spID] = np.random.uniform(0.5, 1.0)

                            # new speciesmaintenance
                            p = np.random.binomial(0.25, 1)
                            MD[spID] = np.random.uniform(0.01, 0.1)

                            # species environmental gradient optima
                            EnvD[spID] = np.random.uniform(0.0, 1.0, len(envGs))

                            # new speciesactive dispersal rate
                            p = np.random.binomial(0.25, 1)
                            DispD[spID] = np.random.uniform(0.0, 0.1)

                            # new speciesresource use efficiencies
                            # Nitrogen
                            N_RD[spID] = np.random.uniform(0.01, 1.0, nN)

                            # Phosphorus
                            P_RD[spID] = np.random.uniform(0.01, 1.0, nP)

                            # Carbon
                            C_RD[spID] = np.random.uniform(0.01, 1.0, nC)

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

                    newX = float(np.random.uniform(X-0.5, X, 1))
                    if limit > newX: newX = 0
                    if newX > w - limit: newX = w - limit
                    Xs.append(newX)

                    newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                    if limit > newY: newY = 0
                    elif newY > h: newY = h - limit
                    Ys.append(newY)

                    newZ = float(np.random.uniform(Z-0.5, Z+0.5, 1))
                    if limit > newZ: newZ = 0
                    elif newZ > l: newZ = l - limit
                    Zs.append(newZ)

    elif repro == 'sexual':

        I_Boxes = []
        spBoxes = []

        if D == 2:
            I_Boxes = [list([]) for _ in xrange(w*h)]
            spBoxes = [list([]) for _ in xrange(w*h)]

        elif D == 3:
            I_Boxes = [list([]) for _ in xrange(w*h*l)]
            spBoxes = [list([]) for _ in xrange(w*h*l)]

        index = 0
        for i, I_ID in enumerate(IDs):

            spID = Sp_IDs[i]
            rX = int(round(Xs[i]))
            rY = int(round(Ys[i]))
            rZ = int(round(Zs[i]))

            if D == 2:
                index = int(round(rX + (rY * w)))
            elif D == 3:
                index = int(round((rY * l * w) + (rX * l) + rZ))

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
                    print 'p2: ',p2
                    # individual not large enough to reproduce
                    if p2 == 0: continue

                    Qs[i2] = q2/2.0
                    Qs.append(q2/2.0)

                    X = Xs[i2]
                    Y = Ys[i2]
                    Z = Zs[i2]

                    newX = float(np.random.uniform(X-0.5, X+0.5, 1))
                    if limit > newX: newX = 0
                    if newX > w - limit: newX = w - limit

                    newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                    if limit > newY: newY = 0
                    elif newY > h: newY = h - limit

                    newZ = float(np.random.uniform(Z-0.5, Z+0.5, 1))
                    if limit > newZ: newZ = 0
                    elif newZ > l: newZ = l - limit

                    ID += 1
                    IDs.append(ID)
                    t_In.append(0)
                    Sp_IDs.append(sp)

    return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, Zs, GD, DispD, GList, MList, NList, PList, CList, DList]



def nonfluid_movement(TypeOf, motion, Lists, xAge, t_In, Xs, Ys, Zs, w, h, l,
            u0, D, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList):

    limit, distance, direction, pop = 0.1, 0, 0, 'no'
    x, y, z = 0, 0, 0
    IDs, Types, Vals = [], [], []
    DispD = {}

    if TypeOf == 'resource':
        Type, IDs, ID, Vals = Lists
    elif TypeOf == 'individual':
        Type, IDs, ID, Vals, DispD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList = Lists
    else:
        IDs = Lists

    if Xs == []:
        if TypeOf == 'tracer':
            return [IDs, Xs, Ys, Zs, xAge, t_In]
        elif TypeOf == 'individual':
            return [Type, Xs, Ys, Zs, xAge, IDs, ID, t_In, Vals, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]
        elif TypeOf == 'resource':
            return [Type, Xs, Ys, Zs, xAge, IDs, ID, t_In, Vals]


    for i, val in enumerate(IDs):

        # get distance
        if TypeOf == 'individual':
            distance = u0
            #np.random.uniform(0, 0.5) #u0 * DispD[Types[i]]
        else:
            distance = u0

        x, y, z = Xs[i], Ys[i], Zs[i]

        # go up or down
        #if TypeOf != 'resource':
        if motion == 'unidirectional':
            direction = np.random.uniform(-1, 1)
        if motion == 'random_walk':
            direction = np.random.uniform(-1, 1)

        y = y + (direction * distance)

        # go forward or backward
        #if TypeOf != 'resource':
        if motion == 'unidirectional':
            direction = np.random.uniform(1, 1)
        if motion == 'random_walk':
            direction = np.random.uniform(-1, 1)

        x = x + (direction * distance)

        if x > w - limit or x < limit:
            pop = 'yes'
        elif y > h - limit or y < limit:
            pop = 'yes'

        # go left or right
        if motion == 'unidirectional':
            direction = np.random.uniform(-0.5, 1.5)
        if motion == 'random_walk':
            direction = np.random.uniform(-1, 1)

        z = z + (direction * distance)

        if z > l - limit or z < limit: pop = 'yes'

        if pop == 'no':
            Xs[i], Ys[i], Zs[i] = x, y, z
            t_In[i] = t_In[i]+1

        elif pop == 'yes':
            xAge.append(t_In[i])
            Xs.pop(i)
            Ys.pop(i)
            Zs.pop(i)
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
        return [IDs, Xs, Ys, Zs, xAge, t_In]
    elif TypeOf == 'individual':
        return [Type, Xs, Ys, Zs, xAge, IDs, ID, t_In, Vals, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList]
    elif TypeOf == 'resource':
        return [Type, Xs, Ys, Zs, xAge, IDs, ID, t_In, Vals]


def search(repro, spec, Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, Zs, w,
        h, l, GD, DispD, colorD, N_RD, P_RD, C_RD, MD, EnvD, envGs, D, nN, nP, nC, GList, MList, NList, PList, CList, DList):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, Zs, GD, DispD, GList, MList, NList, PList, CList, DList]

    for i, spID in enumerate(Sp_IDs):

        X = Xs[i]
        Y = Ys[i]
        Z = Zs[i]

        ex = []
        ey = []
        sp_opts = EnvD[spID]

        for g, opt in enumerate(sp_opts):
            x, y = envGs[g]
            dist = DispD[spID]

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


    return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, Zs, GD, DispD, GList, MList, NList, PList, CList, DList]
