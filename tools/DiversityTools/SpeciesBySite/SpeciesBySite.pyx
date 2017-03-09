import os

mydir = os.path.expanduser("~/data")

OrC = 'closed'

def GetSADsFromBiom_labeled(path, dataset):

    minS = 2

    IN = path + '/' + dataset + '-SSADdata.txt'
    n = sum(1 for line in open(IN))

    SpeciesDict = {}

    print 'Starting build of SpeciesDict'
    with open(IN) as f:

        for d in f:

            #print 'Reading in SSAD data. Lines left:', n
            #n -= 1

            if d.strip():

                d = d.split()
                species = d[0]
                sample = d[1]
                abundance = float(d[2])

                if abundance > 0:
                    if sample not in SpeciesDict:

                        SpeciesDict[sample] = [species]

                    else:
                        SpeciesDict[sample].append(species)

    print 'Finished building SpeciesDict'


    OUT = open(path + '/' + dataset + '-SbyS.txt','w+')

    SpeciesLists = SpeciesDict.items()
    n = len(SpeciesLists)

    for i, site in enumerate(SpeciesLists):

        if len(site) >= minS:
            print n - i
            print >> OUT, site

    OUT.close()

    print 'Finished generating SbyS file'
    return

GetSADsFromBiom_labeled(mydir +'/micro/EMP'+OrC, 'EMP'+OrC)
