from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from random import randint, choice
from scipy import stats
import numpy as np
import sys
import os
import psutil

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "tools/metrics")
import metrics

sys.path.append(mydir + "/GitHub/hydrobide/tools/bide/bideMovie")
import bide


def threeD(args):

    """ local environment as a 3-dimensional space with the maximum coordinate values X, Y, and Z."""

    Xmax = 1000
    Ymax = 1000
    Zmax = 1000

    """ Declare lists to hold spatial coordinates. The reason for using multiple lists is a matter of computing speed
        and efficiency, and will be made clear below.

        Every individual will have a coordinate in one of three lists, depending on whether it is active or dormant.

        Example:
        XcoordsActive[0]  will hold all active individuals with an x coordinate of 0
        YcoordsActive[32] will hold all active individuals with a y coordinate of 32
        ZcoordsActive[81] will hold all active individuals with a z coordinate of 81 """

    XcoordsActive = [list([]) for _ in xrange(X+1)] # creates a list of X+1 empty lists: [[], [], [], ... ]
    YcoordsActive = [list([]) for _ in xrange(Y+1)] # creates a list of Y+1 empty lists
    ZcoordsActive = [list([]) for _ in xrange(Z+1)] # ...and so forth

    XcoordsDormant = [list([]) for _ in xrange(X+1)]
    YcoordsDormant = [list([]) for _ in xrange(Y+1)]
    ZcoordsDormant = [list([]) for _ in xrange(Z+1)]


    """ Declare some stuff """
    J = {}     # The local community will start as an empty dictionary. In Python dictionaries are often
            # more efficient data structures than lists (or lists of lists to be more specific)
    m = 100    # number of expected immigrants
    ID = 0  # a place holder used to add a unique identifier to each individual
    t = 0   # a place holder to track the number of generations that passed


    # Simulate over some number of generations
    while t <= 30:
        print 'generation',t,' J =',len(J)
        """ Immigration of dormant propagules:

            This occurs by randomly choosing individuals to migrate from a regional pool.
            The distribution of abundance for this regional pool (metacommunity) is
            assumed to be a log-series distribution. That is, a canonical hollow-curve
            reflecting that few species are highly abundant and that most are rare.
            This log-series assumption is common, e.g. Hubbell 2001.

            A propagule can land in a spot that is already occupied, but cannot become
            active until the spot is vacated. In this way, dormant propagules can always
            be dumped into the local community, but cannot be competitively dominant in
            the presence of an already active cell. This assumption is just a starting point."""

        lgp = 0.90 # log-series parameter, generally > 0.9 and < 1.0
        propagules = np.random.logseries(lgp, m)  # Using the class 'random' in Numpy to get a
        # random sample of the logseries distribution. 'propagules' will be a list of log-series
        # distributed integers, where each integer represents a different species.

        for i, sp_label in enumerate(propagules):

            prop_info = [sp_label, 0] # 0 denotes that the propagule is dormant

            # We need to assign each individual to a randomly chosen 3D spatial coordinate. This assumes that propagules
            # have equal chances of landing anywhere in the local environment. This does not have to be the case, but
            # it is the simplest scenario)

            x = random.randint(0,X)
            y = random.randint(0,Y)
            z = random.randint(0,Z)

            XcoordsDormant[x].append(ID) #
            YcoordsDormant[y].append(ID) #
            ZcoordsDormant[z].append(ID) #

            prop_info.append([x,y,z]) # append the individual's spatial coordinates
            J[ID] = prop_info # adding the individual's information to the local community so that, for example:
            # J[455] = {455: [the_species_label, Dormant_or_Active, xyz_coordinates]}
            # Then, instead of sorting through a list for individual 455, we can go directly to it.

            ID += 1 # ID will always increase so that no two individuals can ever have the same identifier.

        """ Individual replacement and transitions to and from dormancy:

            A simple scenario where death and emmigration are effectively the same thing, and
            where individual replacement occurs by effectively outcompeting adjacent individuals
            for the opportunity to reproduce. Reproduction occurs via binary fission. One daughter
            cell occupies the coordinate of the parent cell, the other daughter cell occupies a
            randomly chosen adjacent coordinate (replacing any other active individual that might
            already be there). This could all be done differently and deserves some thought and
            exploration.

            I HAVE NOT ADDED THE TRANSITION FROM ACTIVE TO DORMANT; IT WILL TAKE SOME THOUGHT.
            ACTUALLY, THIS MIGHT BE THE PLACE WHERE KAYLA JUMPS IN.

            The transition to activity is induced when a dormant individual is the sole occupant
            of a location (sole competitor for resources at that location).

            Note: Here is where a computational nightmare can occur. This is because we need to
            know whether an individual is about to be replaced, i.e. whether the randomly chosen
            coordinates for a daughter cell are the same as another active member of J. Below are
            slow and fast ways to do this, and the justication for storing X,Y,Z coordinates in
            separate lists:

            THE SLOW MO-FO WAY:
            Get the coordinates for the daughter cell, then loop through the list of community
            members (J) asking whether the coordinates for each individual match those for the
            daughter cell. There can only be one match. If the community is big, then the
            amount of looping will be ridiculous. This will eat up time and greatly limit the
            size of the community you can play with.

            THE FAST 'NEVER-THOUGHT-ABOUT-DOING-IT-THIS-WAY-BEFORE' WAY:
            Store the XYZ coordinates with the individual and also in separate X,Y,Z lists according
            to whether the individual is active or dormant. That way, when we choose an individual we
            can immediately know it's xyz coordinates, which will allow us to quickly choose adjacent
            coordinates for it's daughter cell. Then, we need see what active individual the daughter
            cell is going to replace, if any. Instead of looping through the entire community, we can
            look at the indices of the X Y Z lists that match the coordinates of the daughter cell and
            see whether an active individual already occurs at daughter cell's coordinates.

                HERE'S AN EXAMPLE:
                daughter cell coordinates = [0,0,0]

                XcoordsActive[0] = [2, 4, 7]
                YcoordsActive[0] = [2, 3, 89]
                ZcoordsActive[0] = [2, 4, 3]

                So, individual 2 is active and must have the same coordinates as picked for our daughter cell.
                So we replace individual 2 with the daughter cell, meaning that the parent cell outcompeted
                individual 2 for the chance to reproduce, sensu good ol' neutral theory.

            """

        random.seed() # use current system time to initiate a random seed (ensures that autocorrelation
                    # doesn't crop-up, as is possible when using pseudorandom number generators)
        Jsize = len(J)
        t2 = 0
        while t2 < Jsize: # Induce a generation's worth of replacement events or transitions to/from dormancy
            #print t2, Jsize
            ind = random.choice(J.keys()) # choose an individual at random
            ind_info = J[ind] # get the individual's information
            # Remember that ind is dictionary key pointing to a list.
            # Example  J[455] = {455: [the_species_label, Dormant_or_Active, xyz_coordinates]}

            spLabel = ind_info[0] # the species label
            DorA = ind_info[1] # dormant/active status
            x = ind_info[2][0]
            y = ind_info[2][1]
            z = ind_info[2][2]

            if DorA == 0: # if the inidividual is dormant
                # Find the intersection for the XcoordsActive[x], YcoordsActive[y], and ZcoordsActive[z] lists
                # If there is an intersection, it means an active individual occurs at the xyz coordinates
                # 'intersection' is a list of integers representing individual ID's
                intersection = list(set(XcoordsActive[x]) & set(YcoordsActive[y]) & set(ZcoordsActive[z]))
                if len(intersection) == 0: # if there are no active individuals at the location
                    J[ind][1] = 1 # make the individual active


            elif DorA == 1: # if the individual is active, it will reproduce. One cell will occupy the parent cell's
                                # location and the other will occupy an adjacent location, replacing whatever cell
                                # already occurs there.

                x2 = random.choice([x-1,x+1])
                y2 = random.choice([y-1,y+1])
                z2 = random.choice([z-1,z+1])

                if x2 == -1: x2 = 1 # if one of the coordinates of the parent cell are zero, then the parent cell
                if y2 == -1: y2 = 1 # is already at an edge of the environment. It must reproduce into an
                if z2 == -1: z2 = 1 # available location. '-1' is not a location in the environment

                if x2 == X+1: x2 = int(X) # if one of the coordinates of the parent cell are greater than X,Y, or Z,
                if y2 == Y+1: y2 = int(Y) # then the parent cell is already at an edge of the environment.
                if z2 == Z+1: z2 = int(Z) # It must reproduce into an available location.

                # Find the intersection for the XcoordsActive[x], YcoordsActive[y], and ZcoordsActive[z] lists
                # If there is an intersection, it means an active individual occurs at the xyz coordinates.
                # 'intersection' is a list of integers representing individual ID's
                intersection = list(set(XcoordsActive[x2]) & set(YcoordsActive[y2]) & set(ZcoordsActive[z2]))
                if len(intersection) > 1:
                    print "Error: More than one active individual occupying the exact same spot.\nThe universe is now going to implode. Nice job."
                    sys.exit() # kill the script

                elif len(intersection) == 1:
                    prop_info = [ID, 'A', [x2,y2,z2]] # append the daughter cell's spatial coordinates
                    J[intersection] = prop_info # effectively replacing one individual with another

            t2+=1
            """ Why not add Brownian movement through the local environment? """
        t+=1


    """ The community is now sufficiently simulated. We can go ahead and plot it.
        We'll generate two 3D subplots. One for individuals color-coded by species.
        The other for individuals color-coded by active/dormant. """

    fig = plt.figure()
    ax = fig.add_subplot(221, projection='3d') # 121 means there's one row, two columns, and this is the first plot

    sp_labels = []
    for i in J:
        sp_labels.append(J[i][0])

    sp_labels = list(set(sp_labels))

    colors = {}
    for label in sp_labels:
        r = lambda: random.randint(0,255)
        colors[label] = '#%02X%02X%02X' % (r(),r(),r())

    for ind in J:
        x = J[ind][2][0]
        y = J[ind][2][1]
        z = J[ind][2][2]
        sp_label = J[ind][0]
        ax.scatter(x, y, z, c=colors[sp_label], s=3, marker='o', lw=0)

    plt.tick_params(axis='both', which='major', labelsize=7)
    plt.setp(ax, xticks=[], yticks=[], zticks = [])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax = fig.add_subplot(222, projection='3d') # 121 means there's one row, two columns, and this is the second plot
    colors = ['b','r']
    for ind in J:
        x = J[ind][2][0]
        y = J[ind][2][1]
        z = J[ind][2][2]
        DorA = J[ind][1]
        ax.scatter(x, y, z, c=colors[DorA], s=3, marker='o', lw=0)

    plt.tick_params(axis='both', which='major', labelsize=7)
    plt.setp(ax, xticks=[], yticks=[], zticks = [])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    return [stuff]
