from __future__ import division
import sys

def get_SADs(path, name, closedref=True):

    SADdict = {}
    DATA = path + name + '-data.txt'

    with open(DATA) as f:

        for d in f:

            if d.strip():
                d = d.split()
                length = len(d)

                if name == 'GENTRY':
                    site = d[0]
                    #species = d[1] # Dataset name plus species identifier
                    abundance = float(d[-1])

                else:
                    site = d[0]
                    #year = d[1]

                    if closedref == True:
                        for i in d:
                            if 'unclassified' in i:
                                print 'unclassified'
                                continue
                            elif 'unidentified' in i:
                                print 'unidentified'
                                continue

                    abundance = float(d[-1])


                if abundance > 0:
                    if site in SADdict:
                        SADdict[site].append(abundance)
                    else:
                        SADdict[site] = [abundance]

    SADs = SADdict.values()
    filteredSADs = []
    for sad in SADs:
        if len(sad) >= 10:
            filteredSADs.append(sad)


    return filteredSADs




def EMP_SADs(path, name):

    minS = 10

    IN = path + '/' + name + '-SSADdata.txt'
    n = sum(1 for line in open(IN))

    SiteDict = {}

    with open(IN) as f:

        for d in f:

            n -= 1

            if d.strip():

                d = d.split()
                #species = d[0]
                sample = d[1]
                abundance = float(d[2])

                if abundance > 0:
                    if sample not in SiteDict:

                        SiteDict[sample] = [abundance]

                    else:
                        SiteDict[sample].append(abundance)


    SADs = SiteDict.values()
    filteredSADs = []
    for sad in SADs:
        if len(sad) >= minS:
            filteredSADs.append(sad)

    return filteredSADs
