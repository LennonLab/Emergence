from __future__ import division
from random import randint, choice
import numpy as np
#import sys
#from math import modf
#import decimal
import time

limit = 0.2

def coord(d):
    return float(np.random.uniform(0.1*d, 0.9*d))


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
        rmax, nr, w, h, l, u0, D):

    for r in range(numr):
        x = np.random.binomial(1, u0/2)

        if x == 1:
            rval = int(np.random.random_integers(1, rmax, 1))
            rtype = int(np.random.random_integers(0, nr-1, 1))

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



def immigration(motion, numin, Sp, Xs, Ys, Zs, w, h, l, MD,
        GD, DispD, colorD, IDs, ID, t_In, Qs, RD,
        nr, u0, alpha, D):

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
            Q = float(np.random.uniform(0.05, 0.5))
            Qs.append(Q)

            if prop not in colorD:
                # speciescolor
                colorD = get_color(prop, colorD)

                # speciesgrowth rate
                GD[prop] = np.random.uniform(0.01, 0.1)

                # speciesmaintenance
                MD[prop] = np.random.uniform(0.001, 0.01)

                # speciesactive dispersal rate
                DispD[prop] = np.random.uniform(0.01, 0.1)

                # speciesresource use efficiency
                RD[prop] = np.random.uniform(0.01, 1.0, nr)

    return [Sp, Xs, Ys, Zs, MD, GD, DispD, colorD, IDs, ID, t_In, Qs, RD]



def fluid_movement(TypeOf, List, t_In, xAge, Xs, Ys, ux, uy, w, h, u0):

    Type, IDs, ID, Vals = [], [], int(), []

    if TypeOf == 'resource':
        Type, IDs, ID, Vals = List
    elif TypeOf == 'individual':
        Type, IDs, ID, Vals, DispD = List
    else:
        IDs = List

    if Xs == []:
        if TypeOf == 'resource' or TypeOf == 'individual':
            return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals]

        elif TypeOf == 'tracer':
            return [IDs, Xs, Ys, xAge, t_In]

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
        if TypeOf == 'individual': k = np.random.binomial(1, DispD[Type[i]])

        if k == 0:
            Xs[i] += ux[index]
            Ys[i] += uy[index]

        y = Ys[i]
        if 0 > y: Ys[i] = 0
        elif y > h: Ys[i] = h

        t_In[i] += 1
        if Xs[i] >= w - limit:

            xAge.append(t_In[i])
            Xs.pop(i)
            Ys.pop(i)
            t_In.pop(i)
            IDs.pop(i)

            if TypeOf == 'resource' or TypeOf == 'individual':
                Type.pop(i)
                Vals.pop(i)

    ux = np.reshape(ux, (h, w))
    uy = np.reshape(uy, (h, w))

    if TypeOf == 'tracer':
        return [IDs, Xs, Ys, xAge, t_In]
    elif TypeOf == 'resource' or TypeOf == 'individual':
        return [Type, Xs, Ys, xAge, IDs, ID, t_In, Vals]






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

    List = [P_IDs, P_ID, P_t_In, P_Xs, P_Ys, P_Zs, Sp_IDs, Qs, I_IDs]
    List += [I_ID, I_t_In, I_Xs, I_Ys, I_Zs]
    return List





def maintenance(Sp_IDs, Xs, Ys, Zs, xAge, colorD, MD, IDs, t_In, Qs, D):

    if Sp_IDs == []: return [Sp_IDs, Xs, Ys, Zs, xAge, IDs, t_In, Qs]

    for i, val in enumerate(Qs):
        val -= MD[Sp_IDs[i]]  # maintanence influenced by Sp_
        if val <= 0.01:   # starved

            Qs.pop(i)
            xAge.append(t_In[i])
            t_In.pop(i)
            Sp_IDs.pop(i)
            IDs.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            Zs.pop(i)

        else: Qs[i] = val

    return [Sp_IDs, Xs, Ys, Zs, xAge, IDs, t_In, Qs]



def consume(R_Types, R_Vals, R_IDs, R_ID, R_Xs, R_Ys, R_Zs, R_t_In,
        R_xAge, Sp_IDs, Qs, I_IDs, I_ID, I_t_In, I_Xs, I_Ys, I_Zs,
        w, h, l, GD, RD, DispD, D):

    if not len(R_Types) or not len(Sp_IDs):
        List = [R_Types, R_Vals, R_IDs, R_ID, R_t_In, R_xAge, R_Xs]
        List += [R_Ys, R_Zs, Sp_IDs, Qs, I_IDs, I_ID, I_t_In]
        List += [I_Xs, I_Ys, I_Zs]
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
            R_type = R_Types[j]
            R_val = R_Vals[j]

            # The Individual
            ID = I_IDs.index(ind)
            # The individual's cell quota
            Q = Qs[ID]
            # the species
            sp = Sp_IDs[ID]
            mu = GD[sp]
            efficiency = RD[sp][R_type]
            mu = mu * efficiency

            if R_val > (mu * Q): # Increase cell quota
                R_val = R_val - (mu * Q)
                Q = Q + (mu * Q)

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

            Qs[ID] = Q

    List = [R_Types, R_Vals, R_IDs, R_ID, R_t_In, R_xAge, R_Xs, R_Ys]
    List += [R_Zs,Sp_IDs, Qs, I_IDs, I_ID, I_t_In, I_Xs, I_Ys, I_Zs]
    return List



def reproduce(repro, spec, Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, Zs, w,
        h, l, GD, DispD, colorD, RD, MD, D, nr):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, Zs, GD, DispD]

    if repro == 'fission':
        for i, Q in enumerate(Qs):
            p = np.random.binomial(1, Q)
            if p == 1: # individual is large enough to reproduce
                Qs[i] = Q/2.0
                Qs.append(Q/2.0)
                ID += 1
                IDs.append(ID)
                t_In.append(t_In[i])
                spID = Sp_IDs[i]

                X = Xs[i]
                Y = Ys[i]
                Z = Zs[i]

                if spec == 'yes':
                    p = np.random.binomial(1, 0.05)
                    if p == 1:
                        #print 'new Sp_'

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

                        # new speciesactive dispersal rate
                        p = np.random.binomial(0.25, 1)
                        DispD[spID] = np.random.uniform(0.0, 0.1)

                        # new speciesresource use efficiencies
                        p = np.random.binomial(0.25, 1)
                        RD[spID] = np.random.uniform(0.1, 1.0, nr)

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
                q1 = Qs[index1]

                p1 = np.random.binomial(1, q1)
                # individual not large enough to reproduce
                if p1 == 0: continue

                # Find another of the same Sp_
                if spbox.count(sp) > 1:
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


    return [Sp_IDs, Qs, IDs, ID, t_In, Xs, Ys, Zs, GD, DispD]



def nonfluid_movement(TypeOf, motion, Lists, xAge, t_In, Xs, Ys,
        Zs, w, h, l, u0, D):

    limit, distance, direction, pop = 0.1, 0, 0, 'no'
    X, Y, Z = 0, 0, 0
    IDs, Types, Vals = [], [], []
    DispD = {}

    if TypeOf == 'tracer':
        IDs = Lists
    elif TypeOf == 'individual':
        Types, IDs, Vals, DispD = Lists
    elif TypeOf == 'resource':
        Types, IDs, Vals = Lists

    for i, val in enumerate(IDs):

        # get distance
        if TypeOf == 'individual':
            distance = u0
            #np.random.uniform(0, 0.5) #u0 * DispD[Types[i]]
        else:
            distance = u0

        X, Y, Z = Xs[i], Ys[i], Zs[i]

        # go up or down
        #if TypeOf != 'resource':
        if motion == 'unidirectional':
            direction = np.random.uniform(-1, 1)
        if motion == 'random_walk':
            direction = np.random.uniform(-1, 1)

        Y = Y + (direction * distance)

        # go forward or backward
        #if TypeOf != 'resource':
        if motion == 'unidirectional':
            direction = np.random.uniform(1, 1)
        if motion == 'random_walk':
            direction = np.random.uniform(-1, 1)

        X = X + (direction * distance)

        if X > w - limit or X < limit: pop = 'yes'
        elif Y > h - limit or Y < limit: pop = 'yes'

        # go left or right
        if motion == 'unidirectional':
            direction = np.random.uniform(-0.5, 1.5)
        if motion == 'random_walk':
            direction = np.random.uniform(-1, 1)

        Z = Z + (direction * distance)

        if Z > l - limit or Z < limit: pop = 'yes'

        if pop == 'no':
            Xs[i], Ys[i], Zs[i] = X, Y, Z
            t_In[i] = t_In[i]+1

        elif pop == 'yes':
            IDs.pop(i)
            xAge.append(t_In[i])
            t_In.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            Zs.pop(i)

            if TypeOf == 'individual' or TypeOf == 'resource':
                Vals.pop(i)
                Types.pop(i)

    if TypeOf == 'tracer':
        Lists = IDs
    elif TypeOf == 'individual' or TypeOf == 'resource':
        Lists = [Types, IDs, Vals]

    return [Lists, xAge, t_In, Xs, Ys, Zs]
