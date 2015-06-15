from __future__ import division
from random import randint, choice
import numpy as np
import sys

limit = 0.5

def coord(d):
    return float(np.random.uniform(0.05*d, 0.8*d))


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



def NewTracers(IDs, coords, TimeIn, width, height, length, u0, D):

    Xcoords, Ycoords, Zcoords = [], [], []
    if D == 2:
        Xcoords, Ycoords = coords
    elif D == 3:
        Xcoords, Ycoords, Zcoords = coords

    listlengths = [len(Xcoords), len(Ycoords), len(IDs)]
    if max(listlengths) != min(listlengths):
        print 'Tracer list lengths are different sizes:', listlengths
        sys.exit()

    x = np.random.binomial(1, u0)

    if x == 1:
        TimeIn.append(0)

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

    listlengths = [len(Xcoords), len(Ycoords), len(IDs)]
    if max(listlengths) != min(listlengths):
        print 'Resulting tracer list lengths are different sizes:', listlengths
        sys.exit()

    return [IDs, TimeIn, coords]



def ResIn(Type, Vals, coords, ID, IDs, TimeIn, numr, rmax, nr, width, height, length, u0, D):

    Xcoords, Ycoords, Zcoords = [], [], []
    if D == 2: Xcoords, Ycoords = coords

    elif D == 3:
        Xcoords, Ycoords, Zcoords = coords

        listlengths = [len(Type), len(IDs), len(Vals), len(Xcoords), len(Ycoords), len(Zcoords)]
        if max(listlengths) != min(listlengths):
            print 'Resource list lengths are different sizes:', listlengths
            sys.exit()


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

            y = coord(height)
            Ycoords.append(y)

            x = coord(width)
            Xcoords.append(x)

            if D == 3:
                z = coord(length)
                Zcoords.append(z)

    if D == 2:
        coords = [Xcoords, Ycoords]
        listlengths = [len(Type), len(IDs), len(Vals), len(Xcoords), len(Ycoords)]

    elif D == 3:
        coords = [Xcoords, Ycoords, Zcoords]
        listlengths = [len(Type), len(IDs), len(Vals), len(Xcoords), len(Ycoords), len(Zcoords)]

    if max(listlengths) != min(listlengths):
        print 'ResIn, Resulting resource list lengths are different sizes:', listlengths
        sys.exit()

    return [Type, Vals, coords, IDs, ID, TimeIn]



def immigration(numin, Species, coords, width, height, length, MaintDict, GrowthDict, DispParamDict, color_dict, IDs, ID, TimeIn, Qs, ResUseDict, nr, u0, alpha, D):

    Xcoords, Ycoords, Zcoords = [], [], []
    if D == 2:
        Xcoords, Ycoords = coords
    elif D == 3:
        Xcoords, Ycoords, Zcoords = coords

    for m in range(numin):
        x = np.random.binomial(1, u0)

        if x == 1:
            prop = int(np.random.logseries(0.999, 1))

            if prop > 1000: continue

            Species.append(prop)

            y = coord(height)
            Ycoords.append(y)

            x = coord(width)
            Xcoords.append(x)

            if D == 3:
                z = coord(length)
                Zcoords.append(z)

            IDs.append(ID)
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
                MaintDict[prop] = np.random.uniform(10, 15)

                # species active dispersal rate
                DispParamDict[prop] = np.random.uniform(.02, 0.2)

                # species resource use efficiency
                ResUseDict[prop] = np.random.uniform(0.2, 1.0, nr)

    if D == 2:
        coords = [Xcoords, Ycoords]
    elif D == 3:
        coords = [Xcoords, Ycoords, Zcoords]

    return [Species, coords, MaintDict, GrowthDict, DispParamDict, color_dict, IDs, ID, TimeIn, Qs, ResUseDict]


def fluid_movement(TypeOf, List, TimeIn, ExitAge, Xcoords, Ycoords, ux, uy, width, height, u0):

    Type, IDs, ID, Vals = [], [], int(), []

    if TypeOf == 'resource' or TypeOf == 'individual':
        Type, IDs, ID, Vals = List
        listlengths = [len(Type), len(IDs), len(Vals)]
        if max(listlengths) != min(listlengths):
            print TypeOf, ', List lengths of different sizes in fluid_movement of bide_test.py:', listlengths
            sys.exit()

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
            TimeIn.pop(i)
            IDs.pop(i)

            if TypeOf == 'resource' or TypeOf == 'individual':
                Type.pop(i)
                Vals.pop(i)

    ux = np.reshape(ux, (height, width))
    uy = np.reshape(uy, (height, width))

    listlengths = [len(Xcoords), len(Ycoords), len(IDs)]

    if max(listlengths) != min(listlengths):
        print TypeOf, ', In function fluid_movement, resulting list lengths are different sizes:', listlengths
        sys.exit()

    if TypeOf == 'tracer':
        return [IDs, Xcoords, Ycoords, ExitAge, TimeIn]

    elif TypeOf == 'resource' or TypeOf == 'individual':
        return [Type, Xcoords, Ycoords, ExitAge, IDs, ID, TimeIn, Vals]



def maintenance(Species, coords, ExitAge, color_dict, MaintDict, IDs, TimeIn, Qs, D):

    if Species == []:
        return [Species, coords, ExitAge, IDs, TimeIn, Qs]

    Xcoords, Ycoords, Zcoords = [], [], []
    if D == 2:
        Xcoords, Ycoords = coords
    elif D == 3:
        Xcoords, Ycoords, Zcoords = coords

    for i, val in enumerate(Qs):

        val -= MaintDict[Species[i]]  # maintanence influenced by species

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




def Consume(ResTypes, ResVals, ResIDs, ResID, ResCoords, ResTimeIn, ResExitAge, SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndCoords, width, height, length, GrowthDict, ResUseDict, DispParamsDict, D):

    ResXcoords, ResYcoords, ResZcoords, IndXcoords, IndYcoords, IndZcoords, IndBoxes, ResBoxes = [[], [], [], [], [], [], [], []]

    if D == 2:
        ResXcoords, ResYcoords = ResCoords
        IndXcoords, IndYcoords = IndCoords

    elif D == 3:
        ResXcoords, ResYcoords, ResZcoords = ResCoords
        IndXcoords, IndYcoords, IndZcoords = IndCoords

    if not len(ResTypes):
        ResLists = [ResTypes, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords]
        IndLists = [SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords]
        return [ResLists, IndLists]

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

        if D == 2: index = int(round(roundedX + (roundedY * width)))

        elif D == 3:
            roundedZ = int(round(IndZcoords[i]))
            index = int(round((roundedY * length * width) + (roundedX * length) + roundedZ))

        if index > len(IndBoxes) - 1: index = len(IndBoxes) - 1
        elif index < 0: index = 0

        IndBoxes[index].append(val)

    index = 0
    for i, val in enumerate(ResIDs):
        roundedX = int(round(ResXcoords[i]))
        roundedY = int(round(ResYcoords[i]))

        if D == 2: index = int(round(roundedX + (roundedY * width)))

        elif D == 3:
            roundedZ = int(round(ResZcoords[i]))
            index = int(round((roundedY * length * width) + (roundedX * length) + roundedZ))

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
                ResVals[j] = resval
                #ResBox.pop(boxIndex)

            else:
                Q += resval
                ResBox.pop(boxIndex)
                ResVals.pop(j)
                ResExitAge.append(ResTimeIn[j])
                ResTimeIn.pop(j)
                ResTypes.pop(j)
                ResIDs.pop(j)
                ResXcoords.pop(j)
                ResYcoords.pop(j)
                if D == 3: ResZcoords.pop(j)

            Qs[ID] = Q

    ResLists = ResTypes, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords
    IndLists = SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords

    return [ResLists, IndLists]



def reproduce(reproduction, SpeciesIDs, Qs, IDs, ID, TimeIn, coords, width, height, length, GrowthDict, DispParamsDict, D, qmin):

    if SpeciesIDs == []:
        return [SpeciesIDs, coords, ExitAge, IDs, TimeIn, Qs]

    Xcoords, Ycoords, Zcoords = [], [], []
    if D == 2: Xcoords, Ycoords = coords
    elif D == 3: Xcoords, Ycoords, Zcoords = coords

    for i, Q in enumerate(Qs):
        if reproduction == 'clonal':

            if Q >= qmin: # reproduction
                Qs[i] = Q/2.0
                Qs.append(Q/2.0)
                IndID += 1
                IDs.append(IDs[i])
                TimeIn.append(TimeIn[i])


                SpeciesIDs.append(SpeciesIDs[ID])
                direction = choice([-1,1])
                IndXcoords.append(IndXcoords[ID] + direction*DispParamsDict[SpeciesIDs[ID]])
                direction = choice([-1,1])
                IndYcoords.append(IndYcoords[ID] + direction*DispParamsDict[SpeciesIDs[ID]])

        if D == 3: IndZcoords.append(IndZcoords[ID])






        elif reproduction == 'sexual':






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



    ResLists = ResTypes, ResVals, ResIDs, ResID, ResTimeIn, ResExitAge, ResXcoords, ResYcoords, ResZcoords
    IndLists = SpeciesIDs, Qs, IndIDs, IndID, IndTimeIn, IndXcoords, IndYcoords, IndZcoords

    return [ResLists, IndLists]


def nonfluid_movement(TypeOf, motion, Lists, ExitAge, TimeIn, coords, width, height, length, u0, D):

    limit = 0.5
    Xcoords, Ycoords, Zcoords = [], [], []
    distance, direction, pop = 0, 0, 'no'
    X, Y, Z = 0, 0, 0
    IDs, Types, Vals = [], [], []
    DispParamsDict = {}

    if D == 2: Xcoords, Ycoords = coords
    elif D == 3: Xcoords, Ycoords, Zcoords = coords

    if TypeOf == 'tracer': IDs = Lists
    elif TypeOf == 'individual': Types, IDs, Vals, DispParamsDict = Lists
    elif TypeOf == 'resource': Types, IDs, Vals = Lists

    for i, val in enumerate(IDs):

        # get distance
        if TypeOf == 'individual':
            distance = u0 * DispParamsDict[Types[i]]
        else:
            distance = u0

        X, Y = Xcoords[i], Ycoords[i]

        if D == 3: Z = Zcoords[i]

        # go up or down
        #if TypeOf != 'resource':
        if motion == 'unidirectional': direction = np.random.uniform(-1, 1)
        if motion == 'random_walk': direction = np.random.uniform(-height, height)
        Y = Y + (direction * distance)

        # go forward or backward
        #if TypeOf != 'resource':
        if motion == 'unidirectional': direction = np.random.uniform(1, 1)
        if motion == 'random_walk': direction = np.random.uniform(-width, width)
        X = X + (direction * distance)

        if X > width - limit or X < limit: pop = 'yes'
        elif Y > height - limit or Y < limit: pop = 'yes'

        if D == 3:
            # go left or right
            if motion == 'unidirectional': direction = np.random.uniform(-0.5, 1.5)
            if motion == 'random_walk': direction = np.random.uniform(-length, length)
            Z = Z + (direction * distance)

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

            if D == 3: Zcoords.pop(i)

            if TypeOf == 'individual' or TypeOf == 'resource':
                Vals.pop(i)
                Types.pop(i)

    if D == 2: coords = [Xcoords, Ycoords]
    elif D == 3: coords = [Xcoords, Ycoords, Zcoords]

    if TypeOf == 'tracer': Lists = IDs
    if TypeOf == 'individual' or TypeOf == 'resource': Lists = [Types, IDs, Vals]

    return [Lists, ExitAge, TimeIn, coords]
