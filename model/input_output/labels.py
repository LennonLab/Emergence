from os.path import expanduser


def processes():

    ps = []
    ps.append('resource_inflow')
    ps.append('resource_flow')
    ps.append('immigration')
    ps.append('passive_dispersal')
    ps.append('active_dispersal')
    ps.append('search')
    ps.append('consume')
    ps.append('growth')
    ps.append('transition')
    ps.append('maintenance')
    ps.append('reproduction')
    ps.append('disturb')

    return ps


def headings():

    headings = 'sim,'
    headings += 'ct,'
    headings += 'immigration.rate,'
    headings += 'res.inflow,'
    headings += 'N.types,'
    headings += 'max.res.val,'
    headings += 'max.growth.rate,'
    headings += 'max.met.maint,'
    headings += 'max.active.dispersal,'
    headings += 'starting.seed,'
    headings += 'flow.rate,'
    headings += 'height,'
    headings += 'length,'
    headings += 'width,'
    headings += 'total.abundance,'
    headings += 'ind.production,'
    headings += 'biomass.prod.N,'
    headings += 'resource.particles,'
    headings += 'resource.concentration,'
    headings += 'species.richness,'
    headings += 'simpson.e,'
    headings += 'e.var,'
    headings += 'avg.pop.size,'
    headings += 'pop.var,'
    headings += 'N.max,'
    headings += 'logmod.skew,'
    headings += 'Whittakers.turnover,'
    headings += 'amplitude,'
    headings += 'frequency,'
    headings += 'phase,'
    headings += 'all.biomass,'
    headings += 'active.biomass,'
    headings += 'dormant.biomass,'
    headings += 'all.avg.per.capita.growth,'
    headings += 'active.avg.per.capita.growth,'
    headings += 'dormant.avg.per.capita.growth,'
    headings += 'all.avg.per.capita.maint,'
    headings += 'active.avg.per.capita.maint,'
    headings += 'dormant.avg.per.capita.maint,'
    headings += 'all.avg.per.capita.efficiency,'
    headings += 'active.avg.per.capita.efficiency,'
    headings += 'dormant.avg.per.capita.efficiency,'
    headings += 'all.avg.per.capita.active.dispersal,'
    headings += 'active.avg.per.capita.active.dispersal,'
    headings += 'dormant.avg.per.capita.active.dispersal,'
    headings += 'all.avg.per.capita.rpf,'
    headings += 'active.avg.per.capita.rpf,'
    headings += 'dormant.avg.per.capita.rpf,'
    headings += 'all.avg.per.capita.mf,'
    headings += 'active.avg.per.capita.mf,'
    headings += 'dormant.avg.per.capita.mf,'
    headings += 'all.size,'
    headings += 'active.size,'
    headings += 'dormant.size,'
    headings += 'N.active,'
    headings += 'S.active,'
    headings += 'N.dormant,'
    headings += 'S.Dormant,'
    headings += 'Percent.Dormant,'
    headings += 'dorm.limit'

    return headings



def clear():
    mydir = expanduser("~/")
    GenPath = mydir + 'GitHub/simplex/results/simulated_data/'

    OUT = open(GenPath + 'RAD-Data.csv', 'w+').close()
    OUT = open(GenPath + 'SAR-Data.csv', 'w+').close()
    OUT = open(GenPath + 'SimData.csv','w+')
    h = headings()
    print>>OUT, h
    OUT.close()

    return
