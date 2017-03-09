from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
from scipy import stats


mydir = os.path.expanduser('~/GitHub/residence-time')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')

df2 = pd.DataFrame({'width' : df['width']})
df2['flow'] = df['flow.rate']
df2['tau'] = (df2['width']**3)/df2['flow']
df2['dil'] = np.log10(1/df2['tau'])

df2['N'] = df['total.abundance']
df2['S'] = df['species.richness']
df2['Prod'] = df['ind.production']
df2['E'] = np.log10(df['simpson.e'])
df2['W'] = np.log10(df['Whittakers.turnover'])
df2['Dorm'] = df['Percent.Dormant']

df2['AvgG'] = df['avg.per.capita.growth']
df2['AvgDisp'] = df['avg.per.capita.active.dispersal']
df2['AvgRPF'] = df['avg.per.capita.RPF']
df2['AvgE'] = df['avg.per.capita.N.efficiency']
df2['AvgMaint'] = df['avg.per.capita.maint']
df2['MF'] = df['avg.per.capita.MF']/np.mean(df['avg.per.capita.MF'])

E = 0.04
'''
df2['P'] = df2['AvgMaint'] * df2['AvgRPF'] * df2['MF']
df2['G'] = df2['AvgG'] * df2['AvgDisp'] * df2['AvgE']
df2['phi'] = np.log10(df2['P'] / (df2['G'] + E))
'''

df2['P'] = (1/df2['AvgMaint']) * (1-df2['AvgRPF']) * df2['MF']
df2['G'] = (df2['AvgG'] * df2['AvgDisp'])
df2['phi'] = np.log10(df2['G'] / df2['P'])

df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

#exp_phi = np.mean(df2['phi']) # average of logs
#exp_dil  = np.mean(df2['dil']) # log of the average

print '\nMean and SE for phi:', np.mean(df2['phi']), stats.sem(df2['phi'])
print 'Mean and SE for dil:', np.mean(df2['dil']), stats.sem(df2['dil'])
print (np.abs(-np.mean(df2['phi']) - -np.mean(df2['dil']))/np.mean([-np.mean(df2['phi']), -np.mean(df2['dil'])])) * 100

'''
mP = np.mean(df2['P'])
gP = np.mean(np.log10(df2['P']))
lP = np.log10(np.mean(df2['P']))

mG = np.mean(df2['G'])
gG = np.mean(np.log10(df2['G']))
lG = np.log10(np.mean(df2['G']))

print '1. G:', gG
print '2. P:', gP
print '3. G-P:', gG - gP
print '4. G/P:', gG/gP
print '5. G/P:', gG/gP
print '6. P/G:', gP/gG
print '7. P/G:', gP/gG
print '\n'

print '8. log(G) - E:', gG - E
print '9. log(P) - E:', gP - E
print '10. log(G-P)  - E:', gG - gP - E
print '11. log(G/P)  - E:', gG/gP - E
print '12. log(G)/log(P)  - E:', gG/gP - E
print '13. log(P/G)  - E:', gP/gG - E
print '14. log(P)/log(G)  - E:', gP/gG - E
print '\n'

print '15. G:', lG
print '16. P:', lP
print '17. G-P:', lG - lP
print '18. G/P:', lG/lP
print '19. G/P:', lG/lP
print '20. P/G:', lP/lG
print '21. P/G:', lP/lG
print '\n'

print '22. log(G) - E:', lG - E
print '23. log(P) - E:', lP - E
print '24. log(G-P)  - E:', lG - lP - E
print '25. log(G/P)  - E:', lG/lP - E
print '26. log(G)/log(P)  - E:', lG/lP - E
print '27. log(P/G)  - E:', lP/lG - E
print '28. log(P)/log(G)  - E:', lP/lG - E
print '\n'

print '29. log(G):', np.log10(mG)
print '30. log(P):', np.log10(mP)
print '31. log(G-P):', np.log10(mG - mP)
print '32. log(G/P):', np.log10(mG/mP)
print '33. log(G)/log(P):', np.log10(G)/np.log10(P)
print '34. log(P/G):', np.log10(P/G)
print '35. log(P)/log(G):', np.log10(P)/np.log10(G)
print '\n'

print '36. log(G) - E:', np.log10(G) - E
print '37. log(P) - E:', np.log10(P) - E
print '38. log(G-P)  - E:', np.log10(G - P) - E
print '39. log(G/P)  - E:', np.log10(G/P) - E
print '40. log(G)/log(P)  - E:', np.log10(G)/np.log10(P) - E
print '41. log(P)/log(G)  - E:', np.log10(P)/np.log10(G) - E
print '42. log(P/(G-E)):', np.mean(np.log10((df2['G'] - E)/(df2['P'])))
print '43. log(P/(G-E)):', np.mean(np.log10(df2['P']/(df2['G'] - E)))
'''
