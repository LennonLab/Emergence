from __future__ import division
from random import randint, choice
import numpy as np
import sys
from math import modf
import decimal
import time

limit = 0.2

def coord(d):
    return float(np.random.uniform(0.1*d, 0.9*d))


def GetRAD(vector):
    RAD = []
    unique = list(set(vector))
    for val in unique:
        RAD.append(vector.count(val)) # the abundance of each species

    return RAD, unique # the rad and the species list


def get_color(ID, color_dict): # FUNCTION TO ASSIGN COLORS TO SPECIES

    r1 = lambda: randint(0,255)
    r2 = lambda: randint(0,255)
    r3 = lambda: randint(0,255)

    color = '#%02X%02X%02X' % (r1(),r2(),r3())
    color_dict[ID] = color

    return color_dict



def NewTracers(IDs, Xcoords, Ycoords, Zcoords, TimeIn, width, height, length, u0, D):

    x = np.random.binomial(1, u0)
    if x == 1:
        TimeIn.append(0)
        Ycoords.append(float(np.random.uniform(0.4*height, 0.6*height)))
        Xcoords.append(float(np.random.uniform(0.1*width, 0.11*width)))
        Zcoords.append(float(np.random.uniform(0.4*length, 0.6*length)))
        IDs.append(0) # used to track the age of the tracer

    return [IDs, TimeIn, Xcoords, Ycoords, Zcoords]



def ResIn(motion, Type, Vals, Xcoords, Ycoords, Zcoords, ID, IDs, TimeIn, numr, rmax, nr, width, height, length, u0, D):

    for r in range(numr):
        x = np.random.binomial(1, u0)

        if x == 1:
            rval = int(np.random.random_integers(1, rmax, 1))
            rtype = int(np.random.random_integers(0, nr-1, 1))

            Vals.append(rval)
            IDs.append(ID)
            Type.append(rtype)
            TimeIn.append(0)
            ID += 1

            if motion == 'random_walk':
                Ycoords.append(float(np.random.uniform(0.4*height, 0.6*height)))
                Xcoords.append(float(np.random.uniform(0.4*width, 0.6*width)))
                Zcoords.append(float(np.random.uniform(0.4*length, 0.6*length)))

            else:
                Ycoords.append(float(np.random.uniform(0.1*height, 0.9*height)))
                Xcoords.append(float(np.random.uniform(0.01*width, 0.6*width)))
                Zcoords.append(float(np.random.uniform(0.2*length, 0.8*length)))


    return [Type, Vals, Xcoords, Ycoords, Zcoords, IDs, ID, TimeIn]



def immigration(numin, Species, Xcoords, Ycoords, Zcoords, width, height, length, MaintDict, GrowthDict, DispDict, color_dict, IDs, ID, TimeIn, Qs, ResUseDict, nr, u0, alpha, D):

    for m in range(numin):
        x = np.random.binomial(1, u0)

        if x == 1:
            prop = str(float(np.random.logseries(0.999, 1)))

            Species.append(prop)

            Ycoords.append(float(np.random.uniform(0.1*height, 0.9*height)))
            Xcoords.append(float(np.random.uniform(0.01*width, 0.6*width)))
            Zcoords.append(float(np.random.uniform(0.2*length, 0.8*length)))

            IDs.append(ID)
            TimeIn.append(0)
            ID += 1
            Q = float(np.random.uniform(0.1, 1.0))
            Qs.append(Q)

            if prop not in color_dict:
                # species color
                color_dict = get_color(prop, color_dict)

                # species growth rate
                GrowthDict[prop] = np.random.uniform(0.1, 0.5)

                # species maintenance
                MaintDict[prop] = np.random.uniform(0.001, 0.01)

                # species active dispersal rate
                DispDict[prop] = np.random.uniform(0, 0.2)

                # species resource use efficiency
                ResUseDict[prop] = np.random.uniform(0.1, 1.0, nr)

    return [Species, Xcoords, Ycoords, Zcoords, MaintDict, GrowthDict, DispDict, color_dict, IDs, ID, TimeIn, Qs, ResUseDict]



def fluid_movement(TypeOf, List, TimeIn, ExitAge, Xcoords, Ycoords, ux, uy, width, height, u0):

    Type, IDs, ID, Vals = [], [], int(), []

    if TypeOf == 'resource': Type, IDs, ID, Vals = List
    elif TypeOf == 'individual': Type, IDs, ID, Vals, DispDict = List
    else: IDs = List

    if Xcoords == []:
        if TypeOf == 'resource' or TypeOf == 'individual':
            return [Type, Xcoords, Ycoords, ExitAge, IDs, ID, TimeIn, Vals]

        elif TypeOf == 'tracer':
            return [IDs, Xcoords, Ycoords, ExitAge, TimeIn]

    ux = np.reshape(ux, (width*height)) # ux is the macroscopic x velocity
    uy = np.reshape(uy, (width*height)) # uy is the macroscopic y velocity

    # dispersal inside the system
    for i, val in enumerate(Xcoords):

        X = int(round(Xcoords[i]))
        Y = int(round(Ycoords[i]))

        index =  int(round(X + Y * width))

        if index > len(ux) - 1: index = len(ux) - 1
        if index > len(uy) - 1: index = len(uy) - 1

        k = 0
        if TypeOf == 'individual': k = np.random.binomial(1, DispDict[Type[i]])

        if k == 0:
            Xcoords[i] += ux[index]
            Ycoords[i] += uy[index]

        y = Ycoords[i]
        if 0 > y: Ycoords[i] = 0
        elif y > height: Ycoords[i] = height

        TimeIn[i] += 1
        if Xcoords[i] >= width - limit:

            ExitAge.append(TimeIn[i])
            Xcoords.pop(i)
            Ycoords.pop(i)
            TimeIn.pop(i)
            IDs.pop(i)

            if TypeOf == 'resource' or TypeOf == 'individual':
                Type.pop(i)
                Vals.pop(i)

    ux = np.reshape(ux, (height, width))
    uy = np.reshape(uy, (height, width))

    if TypeOf == 'tracer':
        return [IDs, Xcoords, Ycoords, ExitAge, TimeIn]
    elif TypeOf == 'resource' or TypeOf == 'individual':
        return [Type, Xcoords, Ycoords, ExitAge, IDs, ID, TimeIn, Vals]



def maintenance(SpeciesIDs, Xcoords, Ycoords, Zcoords, ExitAge, color_dict, MaintDict, IDs, TimeIn, Qs, D):

    if SpeciesIDs == []: return [SpeciesIDs, Xcoords, Ycoords, Zcoords, ExitAge, IDs, TimeIn, Qs]

    for i, val in enumerate(Qs):
        val -= MaintDict[SpeciesIDs[i]]  # maintanence influenced by species
        if val <= 0.01:   # starved

            Qs.pop(i)
            ExitAge.append(TimeIn[i])
            TimeIn.pop(i)
            SpeciesIDs.pop(i)
            IDs.pop(i)
            Xcoords.pop(i)
            Ycoords.pop(i)
            Zcoords.pop(i)

        else: Qs[i] = val

    return [SpeciesIDs, Xcoords, Ycoords, Zcoords, ExitAge, IDs, TimeIn, Qs]



def predation(PredIDs, PredID, PredXcoords, PredYcoords, PredZcoords, PredTimeIn, IndExitAge, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords, width, height, length, D):

    IndBoxes, PredBoxes = [], []
    if not len(PredIDs): return [PredIDs, PredID, PredTimeIn, PredXcoords, PredYcoords, PredZcoords, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords]

    if D == 2: IndBoxes, PredBoxes = [[list([]) for _ in xrange(width*height)]]*2
    elif D == 3: IndBoxes, PredBoxes = [[list([]) for _ in xrange(width*height*length)]]*2

    index = 0
    for i, val in enumerate(IndIDs):
        roundedX = int(round(IndXcoords[i]))
        roundedY = int(round(IndYcoords[i]))
        roundedZ = int(round(IndZcoords[i]))

        if D == 2: index = int(round(roundedX + (roundedY * width)))
        elif D == 3: index = int(round((roundedY * length * width) + (roundedX * length) + roundedZ))

        if index > len(IndBoxes) - 1: index = len(IndBoxes) - 1
        elif index < 0: index = 0
        IndBoxes[index].append(val)

    index = 0
    for i, val in enumerate(PredIDs):
        roundedX = int(round(PredXcoords[i]))
        roundedY = int(round(PredYcoords[i]))
        roundedZ = int(round(PredZcoords[i]))

        if D == 2: index = int(round(roundedX + (roundedY * width)))
        elif D == 3: index = int(round((roundedY * length * width) + (roundedX * length) + roundedZ))

        if index > len(PredBoxes) - 1: index = len(PredBoxes) - 1
        elif index < 0: index = 0
        PredBoxes[index].append(val)


    for i, PredBox in enumerate(PredBoxes):
        if not len(PredBox): continue

        IndBox = IndBoxes[i]

        for pred in PredBox:
            if not len(IndBox): break

            # The predator
            PredBox.pop(0)

            # The prey
            IndID = choice(IndBox)
            index = IndBox.index(IndID)
            IndBox.pop(index)

            j = IndIDs.index(IndID)
            Qs.pop(j)
            IndExitAge.append(IndTimeIn[j])
            IndTimeIn.pop(j)
            SpeciesIDs.pop(j)
            IndIDs.pop(j)
            IndXcoords.pop(j)
            IndYcoords.pop(j)
            IndZcoords.pop(j)

    return [PredIDs, PredID, PredTimeIn, PredXcoords, PredYcoords, PredZcoords, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords]



def consume(ResTypes, ResVals, ResIDs, ResID, ResXcoords, ResYcoords, ResZcoords, ResTimeIn, ResExitAge, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords, width, height, length, GrowthDict, ResUseDict, DispDict, D):

    if not len(ResTypes) or not len(SpeciesIDs): return [ResTypes, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords]

    if D == 2:
        IndBoxes = [list([]) for _ in xrange(width*height)]
        ResBoxes = [list([]) for _ in xrange(width*height)]

    elif D == 3:
        IndBoxes = [list([]) for _ in xrange(width*height*length)]
        ResBoxes = [list([]) for _ in xrange(width*height*length)]

    index = 0
    for i, val in enumerate(IndIDs):
        roundedX = int(round(IndXcoords[i]))
        roundedY = int(round(IndYcoords[i]))
        roundedZ = int(round(IndZcoords[i]))

        if D == 2: index = int(round(roundedX + (roundedY * width)))
        elif D == 3: index = int(round((roundedY * length * width) + (roundedX * length) + roundedZ))

        if index > len(IndBoxes) - 1: index = len(IndBoxes) - 1
        elif index < 0: index = 0

        IndBoxes[index].append(val)

    index = 0
    for i, val in enumerate(ResIDs):

        roundedX = int(round(ResXcoords[i]))
        roundedY = int(round(ResYcoords[i]))
        roundedZ = int(round(ResZcoords[i]))

        if D == 2: index = int(round(roundedX + (roundedY * width)))
        elif D == 3: index = int(round((roundedY * length * width) + (roundedX * length) + roundedZ))

        if index > len(ResBoxes) - 1: index = len(ResBoxes) - 1
        elif index < 0: index = 0

        ResBoxes[index].append(val)


    for i, box in enumerate(IndBoxes):
        if not len(box): continue

        ResBox = ResBoxes[i]

        for ind in box: # The individuals
            if not len(ResBox): break

            resID = choice(ResBox)
            boxIndex = ResBox.index(resID)

            # The food
            j = ResIDs.index(resID)
            restype = ResTypes[j]
            resval = ResVals[j]

            # The Individual
            ID = IndIDs.index(ind)
            Q = Qs[ID]
            species = SpeciesIDs[ID]
            mu = GrowthDict[species]
            efficiency = ResUseDict[species][restype]
            mu = mu * efficiency

            if resval > (mu * Q): # Increase cell quota
                resval = resval - (mu * Q)
                Q = Q + (mu * Q)

            else:
                Q += resval
                resval = 0.0

            if Q > 1.0:
                resval = Q - 1.0
                Q = 1.0
                ResVals[j] = resval


            if resval <= 0.0:
                ResBox.pop(boxIndex)
                ResVals.pop(j)
                ResExitAge.append(ResTimeIn[j])
                ResTimeIn.pop(j)
                ResTypes.pop(j)
                ResIDs.pop(j)
                ResXcoords.pop(j)
                ResYcoords.pop(j)
                ResZcoords.pop(j)

            Qs[ID] = Q

    return [ResTypes, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords]



def reproduce(reproduction, speciation, SpeciesIDs, Qs, IDs, ID, TimeIn, Xcoords, Ycoords, Zcoords, width, height, length, GrowthDict, DispDict, color_dict, ResUseDict, MaintDict, D, nr):

    if SpeciesIDs == []:
        return [SpeciesIDs, Qs, IDs, ID, TimeIn, Xcoords, Ycoords, Zcoords, GrowthDict, DispDict]

    if reproduction == 'fission':
        for i, Q in enumerate(Qs):
            p = np.random.binomial(1, Q)
            if p == 1: # individual is large enough to reproduce
                Qs[i] = Q/2.0
                Qs.append(Q/2.0)
                ID += 1
                IDs.append(ID)
                TimeIn.append(TimeIn[i])
                spID = SpeciesIDs[i]

                X = Xcoords[i]
                Y = Ycoords[i]
                Z = Zcoords[i]

                if speciation == 'yes':
                    p = np.random.binomial(1, 0.01)
                    if p == 1:
                        #print 'new species'

                        # speciate
                        t = time.clock()
                        spID = spID +' '+ str(t)

                        # new species color
                        color_dict = get_color(spID, color_dict)

                        # new species growth rate
                        p = np.random.binomial(0.25, 1)
                        GrowthDict[spID] = np.random.uniform(0.5, 1.0)

                        # new species maintenance
                        p = np.random.binomial(0.25, 1)
                        MaintDict[spID] = np.random.uniform(0.01, 0.1)

                        # new species active dispersal rate
                        p = np.random.binomial(0.25, 1)
                        DispDict[spID] = np.random.uniform(0.0, 0.1)

                        # new species resource use efficiencies
                        p = np.random.binomial(0.25, 1)
                        ResUseDict[spID] = np.random.uniform(0.1, 1.0, nr)

                SpeciesIDs.append(spID)

                newX = float(np.random.uniform(X-0.5, X, 1))
                if limit > newX: newX = 0
                if newX > width - limit: newX = width - limit
                Xcoords.append(newX)

                newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                if limit > newY: newY = 0
                elif newY > height: newY = height - limit
                Ycoords.append(newY)

                newZ = float(np.random.uniform(Z-0.5, Z+0.5, 1))
                if limit > newZ: newZ = 0
                elif newZ > length: newZ = length - limit
                Zcoords.append(newZ)

    elif reproduction == 'sexual':

        indBoxes = []
        spBoxes = []

        if D == 2:
            indBoxes = [list([]) for _ in xrange(width*height)]
            spBoxes = [list([]) for _ in xrange(width*height)]

        elif D == 3:
            indBoxes = [list([]) for _ in xrange(width*height*length)]
            spBoxes = [list([]) for _ in xrange(width*height*length)]

        index = 0
        for i, indID in enumerate(IDs):

            spID = SpeciesIDs[i]
            roundedX = int(round(Xcoords[i]))
            roundedY = int(round(Ycoords[i]))
            roundedZ = int(round(Zcoords[i]))

            if D == 2: index = int(round(roundedX + (roundedY * width)))
            elif D == 3: index = int(round((roundedY * length * width) + (roundedX * length) + roundedZ))

            if index > len(indBoxes) - 1: index = len(indBoxes) - 1
            elif index < 0: index = 0

            indBoxes[index].append(indID)
            spBoxes[index].append(spID)

        for i, box in enumerate(indBoxes):
            if len(box) < 2: continue

            spbox = spBoxes[i]
            while len(box) > 1:
                i1 = choice(range(len(box))) # choose an individual at random
                ind1 = box.pop(i1)
                sp = spbox.pop(i1) # remove the species ID

                index1 = IDs.index(ind1)
                q1 = Qs[index1]

                p1 = np.random.binomial(1, q1)
                if p1 == 0: continue # individual not large enough to reproduce

                # Find another of the same species
                if spbox.count(sp) > 1:
                    i2 = spbox.index(sp) # choose an individual of the same species
                    spbox.pop(i2) # remove the species ID
                    box.pop(i2) # remove the individual

                    index2 = IDs.index(i2)
                    q2 = Qs[index2]

                    p2 = np.random.binomial(1, q2)
                    if p2 == 0: continue # individual not large enough to reproduce

                    Qs[i2] = q2/2.0
                    Qs.append(q2/2.0)

                    X = Xcoords[i2]
                    Y = Ycoords[i2]
                    Z = Zcoords[i2]

                    newX = float(np.random.uniform(X-0.5, X+0.5, 1))
                    if limit > newX: newX = 0
                    if newX > width - limit: newX = width - limit

                    newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                    if limit > newY: newY = 0
                    elif newY > height: newY = height - limit

                    newZ = float(np.random.uniform(Z-0.5, Z+0.5, 1))
                    if limit > newZ: newZ = 0
                    elif newZ > length: newZ = length - limit

                    ID += 1
                    IDs.append(ID)
                    TimeIn.append(0)
                    SpeciesIDs.append(sp)


    return [SpeciesIDs, Qs, IDs, ID, TimeIn, Xcoords, Ycoords, Zcoords, GrowthDict, DispDict]



def nonfluid_movement(TypeOf, motion, Lists, ExitAge, TimeIn, Xcoords, Ycoords, Zcoords, width, height, length, u0, D):

    limit, distance, direction, pop = 0.1, 0, 0, 'no'
    X, Y, Z = 0, 0, 0
    IDs, Types, Vals = [], [], []
    DispDict = {}

    if TypeOf == 'tracer': IDs = Lists
    elif TypeOf == 'individual': Types, IDs, Vals, DispDict = Lists
    elif TypeOf == 'resource': Types, IDs, Vals = Lists

    for i, val in enumerate(IDs):

        # get distance
        if TypeOf == 'individual': distance = u0 #np.random.uniform(0, 0.5) #u0 * DispDict[Types[i]]
        else: distance = u0

        X, Y, Z = Xcoords[i], Ycoords[i], Zcoords[i]

        # go up or down
        #if TypeOf != 'resource':
        if motion == 'unidirectional': direction = np.random.uniform(-1, 1)
        if motion == 'random_walk': direction = np.random.uniform(-1, 1)
        Y = Y + (direction * distance)

        # go forward or backward
        #if TypeOf != 'resource':
        if motion == 'unidirectional': direction = np.random.uniform(1, 1)
        if motion == 'random_walk': direction = np.random.uniform(-1, 1)
        X = X + (direction * distance)

        if X > width - limit or X < limit: pop = 'yes'
        elif Y > height - limit or Y < limit: pop = 'yes'

        # go left or right
        if motion == 'unidirectional': direction = np.random.uniform(-0.5, 1.5)
        if motion == 'random_walk': direction = np.random.uniform(-1, 1)
        Z = Z + (direction * distance)
        if Z > length - limit or Z < limit: pop = 'yes'


        if pop == 'no': Xcoords[i], Ycoords[i], Zcoords[i] = X, Y, Z
        elif pop == 'yes':
            IDs.pop(i)
            ExitAge.append(TimeIn[i])
            TimeIn.pop(i)
            Xcoords.pop(i)
            Ycoords.pop(i)
            Zcoords.pop(i)

            if TypeOf == 'individual' or TypeOf == 'resource':
                Vals.pop(i)
                Types.pop(i)

    if TypeOf == 'tracer': Lists = IDs
    elif TypeOf == 'individual' or TypeOf == 'resource': Lists = [Types, IDs, Vals]

    return [Lists, ExitAge, TimeIn, Xcoords, Ycoords, Zcoords]
