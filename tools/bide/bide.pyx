from __future__ import division
from random import randint, choice
import numpy as np
import sys
from scipy import stats
import os

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "tools/metrics")
import metrics

limit = 0.5


def coord(d):
    return float(np.random.uniform(0.2*d, 0.8*d))


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



def NewTracers(IDs, coords, width, height, length, u0, D):

    if D == 2:
        Xcoords, Ycoords = coords
    elif D == 3:
        Xcoords, Ycoords, Zcoords = coords

    x = np.random.binomial(1, u0)

    if x == 1:
        y = coord(height)
        Ycoords.append(y)

        x = coord(width)
        Xcoords.append(x)

        IDs.append(0) # used to track the age of the tracer

    if D == 3:
        z = coord(length)
        Zcoords.append(z)

    if D == 2:
        coords = [Xcoords, Ycoords]
    elif D == 3:
        coords = [Xcoords, Ycoords, Zcoords]

    return [TracerIDs, coords]



def ResIn(Type, Vals, coords, ID, IDs, ExitAge, TimeIn, numr, rmax, nr, width, height, length, u0, D):

    if D == 2:
        Xcoords, Ycoords = coords
    elif D == 3:
        Xcoords, Ycoords, Zcoords = coords

    for r in numr:
        x = np.random.binomial(1, u0)

        if x == 1:
            rval = np.random.random_integers(1, rmax, 1)
            rtype = np.random.random_integers(0, nr-1, 1)

            Vals.append(val)
            IDs.append(ID)
            Type.append(rtype)
            ID += 1

            y = coord(height)
            Ycoords.append(y)

            x = coord(width)
            Xcoords.append(x)

            if D == 3:
                z = coord(length)
                Zcoords.append(z)

    if D == 2:
        coords = [Xcoords, Ycoords]
    elif D == 3:
        coords = [Xcoords, Ycoords, Zcoords]

    return [Type, Vals, coords, ExitAge, IDs, ID, TimeIn]



def immigration(numin, Species, coords, width, height, length, MaintDict, GrowthDict, DispParamDict, color_dict, IDs, ID, TimeIn, Qs, ResUseDict, nr, u0, alpha, D):

    if D == 2:
        Xcoords, Ycoords = coords
    elif D == 3:
        Xcoords, Ycoords, Zcoords = coords

    N = len(COM)

    for m in numin:
        x = np.random.binomial(1, u0)

        if x == 1:
            prop = np.random.logseries(alpha, 1)
            #props = list(np.random.randint(1, 1000, 1))

            if prop > 1000: continue

            Species.append(prop)

            y = coord(height)
            Ycoords.append(y)

            x = coord(width)
            Xcoords.append(x)

            if D == 3:
                z = coord(length)
                Zcoords.append(z)

            IDs.append(IndID)
            TimeIn.append(0)
            ID += 1
            Q = float(np.random.uniform(50, 100))
            Qs.append(Q)

            if prop not in color_dict:
                # species color
                color_dict = get_color(prop, color_dict)

                # species growth rate
                GrowthDict[prop] = np.random.uniform(.5, 1)

                # species maintenance
                MaintDict[prop] = np.random.uniform(5, 15)

                # species active dispersal rate
                DispParamDict[prop] = np.random.uniform(.2, 2)

                # species resource use efficiency
                ResUseDict[prop] = np.random.uniform(0, 1, nr)

    if D == 2:
        coords = [Xcoords, Ycoords]
    elif D == 3:
        coords = [Xcoords, Ycoords, Zcoords]

    return [Species, coords, width, height, MaintDict, GrowthDict, DispParamDict, color_dict, IDs, ID, TimeIn, Qs, ResUseDict]




def fluid_movement(Type, IDlist, ID, TimeIn, Vals, ExitAge, ux, uy, Xcoords, Ycoords, width, height, u0):

    ux = np.reshape(ux, (width*height)) # ux is the macroscopic x velocity
    uy = np.reshape(uy, (width*height)) # uy is the macroscopic y velocity

    # dispersal inside the system
    for i, val in enumerate(Xcoords):

        X = int(round(Xcoords[i]))
        Y = int(round(Ycoords[i]))

        index =  int(round(X + Y * width))

        if index > len(ux) - 1: index = len(ux) - 1
        if index > len(uy) - 1: index = len(uy) - 1

        k = np.random.binomial(1, 1)
        if k == 1:
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
            IDs.pop(i)
            Vals.pop(i)
            TimeIn.pop(i)
            Type.pop(i)

    ux = np.reshape(ux, (height, width))
    uy = np.reshape(uy, (height, width))

    return [Type, Xcoords, Ycoords, ExitAge, IDs, ID, TimeIn, Vals]



def maintenance(Species, coords, ExitAge, color_dict, MaintDict, IDs, TimeIn, Qs, D):

    if D == 2:
        Xcoords, Ycoords = coords
    elif D == 3:
        Xcoords, Ycoords, Zcoords = coords

    for i, val in enumerate(Qs):

        val -= MaintDict[SpeciesTypes[i]]  # maintanence influenced by species

        if val < 1:   # starved

            Qs.pop(i)
            ExitAge.append(TimeIn[i])
            TimeIn.pop(i)
            Species.pop(i)
            IDs.pop(i)

            Xcoords.pop(i)
            Ycoords.pop(i)

            if D == 3:
                Zcoords.pop(i)

        else: Qs[i] = val

    if D == 2:
        coords = [Xcoords, Ycoords]
    elif D == 3:
        coords = [Xcoords, Ycoords, Zcoords]

    return [Species, coords, ExitAge, IDs, TimeIn, Qs]




def ConsumeAndReproduce(ResType, ResVals, ResIDs, ResID, ResCoords, ResTimeIn, ResExitAge, IndType, IndQs, IndIDs, IndID, IndTimeIn, IndCoords, width, height, length, GrowthDict, ResUseDict):

    if D == 2:
        ResXcoords, ResYcoords = ResCoords
        IndXcoords, IndYcoords = IndCoords

    elif D == 3:
        ResXcoords, ResYcoords, ResZcoords = ResCoords
        IndXcoords, IndYcoords, IndZcoords = IndCoords

    if D == 2:
        IndBoxes = [list([]) for _ in xrange(width*height)]
        ResBoxes = [list([]) for _ in xrange(width*height)]

    elif D == 3:
        IndBoxes = [list([]) for _ in xrange(width*height*length)]
        ResBoxes = [list([]) for _ in xrange(width*height*length)]

    for i, val in enumerate(IndIDs):

        roundedX = int(round(IndXcoords[i]))
        roundedY = int(round(IndYcoords[i]))

        if D == 2:
            index = int(round(roundedX + (roundedY * width)))

        if D == 3:
            roundedZ = int(round(IndZcoords[i]))
            index = int(round((roundedY * length * width) + (roundedX * length) + roundedZ))

        if index > len(IndBoxes) - 1:
            index = len(IndBoxes) - 1
        elif index < 0:
            index = 0

        IndBoxes[index].append(val)

    for i, val in enumerate(ResIDs):

        roundedX = int(round(ResXcoords[i]))
        roundedY = int(round(ResYcoords[i]))

        if D == 2:
            index = int(round(roundedX + (roundedY * width)))

        elif D == 3:
            roundedZ = int(round(ResZcoords[i]))
            index = int(round((roundedY * length * width) + (roundedX * length) + roundedZ))

        if index > len(ResBoxes) - 1:
            index = len(ResBoxes) - 1
        elif index < 0:
            index = 0

        ResBoxes[index].append(val)


    for i, box in enumerate(IndBoxes):
        ResBox = ResBoxes[i]

        while box: # The resource
            if len(ResBox): pass
            else: break

            resID = choice(ResBox)

            # The food
            j = ResIDs.index(resID)
            restype = ResType[j]
            resval = ResVal[j]

            # The Individual
            spID = choice(box)
            box.remove(spID)
            index = IndIDs.index(spID)

            Q = IndQs[index]
            species = IndType[index]
            mu = GrowthDict[species]
            efficiency = ResUseDict[species][restype]
            mu = mu * efficiency

            if resval > mu * Q: # Increase cell quota
                Q += mu * Q
                resval -= mu * Q
                ResBoxes[i][j] = resval
                ResVal[j] = resval
                ResTimeIn[j] += 1

            else:
                Q += resval
                ResBoxes[i].pop(j)

                ResXcoords.pop(j)
                ResYcoords.pop(j)
                ResZcoords.pop(j)

                ResVal.pop(j)
                ResTimeIn.pop(j)
                ResType.pop(j)
                ResIDs.pop(j)
                ResExitAge.append(ResTimeIn[j])

            IndQs[index] = Q

            if Q > 500: # reproduction

                spID = IndType[index]
                X =IndXcoords[index]
                Y =IndYcoords[index]
                Z =IndZcoords[index]

                IndQs[index] = Q/2.0

                newX = float(np.random.uniform(X-0.5, X+0.5, 1))
                if newX > width - limit:
                    newX = width - limit

                newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                if 0 > newY: newY = 0
                elif newY > height: newY = height

                IndXcoords.append(newX)
                IndYcoords.append(newY)

                if D == 3:
                    newZ = float(np.random.uniform(Z-0.5, Z+0.5, 1))
                    if 0 > newZ: newZ = 0
                    elif newZ > height: newZ = height
                    IndZcoords.append(newZ)

                IndQs.append(Q/2.0)
                IndType.append(spID)

                IndIDs.append(IndID)
                IndTimeIn.append(0)
                IndID += 1

    IndBoxes, ResBoxes = [], []

    ResLists = ResType, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords
    IndLists = IndType, IndQs,   IndIDs, IndID, IndTimeIn,             IndXcoords, IndYcoords, IndZcoords

    return [ResLists, IndLists]



def nonfluid_movement(Type, Lists, ExitAge, TimeIn, coords, width, height, length, u0, D):

    if D == 2:
        Xcoords, Ycoords = coords
    elif D == 3:
        Xcoords, Ycoords, Zcoords = coords

    if Type == 'tracer':
        IDs = Lists

    elif Type == 'individual':
        Types, IDs, Vals, DispParamDict = Lists

    elif Type == 'resource':
        Types, IDs, Vals = Lists

    for i, val in Types:

        distance, direction, pop = [0, 0, 'no']

        # get distance
        if Type == 'individual':
            distance = u0 * DispParamsDict[val]
        else: distance = u0

        X, Y = Xcoords[i], Ycoords[i]
        # go up or down
        direction = choice([-1, 1])
        Y += direction * distance

        # go forward or backward
        direction = choice([-1, 1])
        X += direction * distance

        if X > width - limit or X < limit:
            pop = 'yes'

        elif Y > height - limit or Y < limit:
            pop = 'yes'

        if D == 3 and pop == 'no':
            Z = Zcoords[i]
            # go left or right
            direction = choice([-1, 1])
            Z += direction * distance

            if Z > length - limit or Z < limit: pop = 'yes'

        if pop == 'no':
            Xcoords[i] = X
            Ycoords[i] = Y

            if D == 3: Zcoords[i] = Z

        elif pop == 'yes':

            IDs.pop(i)
            ExitAge.append(TimeIn[i])
            TimeIn.pop(i)

            Xcoords.pop(i)
            Ycoords.pop(i)

            if D == 3:
                Zcoords.pop(i)

            if Type == 'individual' or Type == 'resource':
                Vals.pop(i)
                Types.pop(i)


    if D == 2:
        coords = [Xcoords, Ycoords]
    elif D == 3:
        coords = [Xcoords, Ycoords, Zcoords]

    if Type == 'tracer':
        Lists = [IDs]

    if Type == 'individual' or Type == 'resource':
        Lists = [Types, IDs, Vals]

    return [Lists, ExitAge, TimeIn, coords]
