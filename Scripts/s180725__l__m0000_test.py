"""
180725
Testing m0000

"""

import M as x
import Models.m0000 as m0000

x.datadirec = '../../ModelaData'

########## Single simulaton <- good

model = m0000.Model(m0000.p0)
x.alg_singlesim(model, 2, 0, 0, compression=0)
x.sliderplot(2, 0, 0)
#
# model = m.Model(m.p0)
# model.params.pA = 0
# alg_singlesim(model, 2, 0, 1, compression=0)
# sliderplot(2, 0, 1)


########## Data compression <- good

# model = m.Model(m.p0)
# alg_singlesim(model, 2, 0, 0, compression=1)
# sliderplot(2, 0, 0)
#
# model = m.Model(m.p0)
# model.params.pA = 0
# alg_singlesim(model, 2, 0, 1, compression=1)
# sliderplot(2, 0, 1)



########## Parallel simulation <- good

# params = ['eAneg', 'ePneg']
# vals = [[1, 2], [1, 2]]
# alg_parsim(m.Model(m.p0), params=params, vals=vals, jobid=2, subjobid=2)
#
# sliderplot(2, 2, 0)
# sliderplot(2, 2, 1)
# sliderplot(2, 2, 2)
# sliderplot(2, 2, 3)


########### Parallel, random <- good

# params = ['pA', 'pP']
# ranges = [[0, 2], [0, 2]]
# alg_parsim_rand(m.Model(m.p0), params=params, ranges=ranges, nsims=4, jobid=2, subjobid=2)
#
# sliderplot(2, 2, 0)
# sliderplot(2, 2, 1)
# sliderplot(2, 2, 2)
# sliderplot(2, 2, 3)


########### Parallel, random, with *args, compression=2  <- good

params = ['pA', 'pP']
ranges = [[0, 2], [0, 2]]
alg_parsim_rand(m0000.Model(m0000.p0), params=params, ranges=ranges, nsims=4, jobid=2, subjobid=2, compression=2, funcs=[mse, asi_a, asi_p])
print(loaddata(2, 2, 0).params.pA)
print(loaddata(2, 2, 0).params.pP)
print(loaddata(2, 2, 0).scores)


########## Genetic algorithm <- sort of works, but has a tendancy to bootleneck

# params = ['pA', 'pP']
# ranges = [[0, 2], [0, 2]]
# gen_alg(m=m.Model(m.p0), func=func_genalg_0000, params=params, ranges=ranges, pop=16, gens=10, jobid=3)
# genalg_plot1(3)


########## Genetic algorithm 2

# alg_singlesim(m.Model(m.p0), 9999, 0, 0)
#
# params = ['pA', 'pP']
# ranges = [[0, 2], [0, 2]]
# gen_alg(m=m.Model(m.p0), func=func_genalg_0001, params=params, ranges=ranges, pop=16, gens=5, jobid=4)
# genalg_plot1(3, [9999, 0, 0])

from pprint import pprint
a = loaddata(2, 2, 0)
pprint(vars(a))


