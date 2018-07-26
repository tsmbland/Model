"""
180725
Testing m0008

"""

from M import *
import Models.m0008 as m

######### Single simulation <- good

# model = m.Model(m.p0)
# alg_singlesim(model, 2, 0, 0, compression=0)
# sliderplot(2, 0, 0)


########## Data compression <- good

# model = m.Model(m.p0)
# alg_singlesim(model, 2, 0, 0, compression=1)
# sliderplot(2, 0, 0)



########## Parallel simulation <- good

# params = ['kon_a', 'kon_p']
# vals = [[1, 2], [1, 2]]
# alg_parsim(m.Model(m.p0), params=params, vals=vals, jobid=2, subjobid=2)
#
# sliderplot(2, 2, 0)
# sliderplot(2, 2, 1)
# sliderplot(2, 2, 2)
# sliderplot(2, 2, 3)


########### Parallel, random <- good

# params = ['kon_a', 'kon_p']
# ranges = [[0, 2], [0, 2]]
# alg_parsim_rand(m.Model(m.p0), params=params, ranges=ranges, nsims=4, jobid=2, subjobid=2)
#
# sliderplot(2, 2, 0)
# sliderplot(2, 2, 1)
# sliderplot(2, 2, 2)
# sliderplot(2, 2, 3)


########## Genetic algorithm <- sort of works, but has a tendancy to bootleneck

# params = ['kon_a', 'kon_p']
# ranges = [[0, 2], [0, 2]]
# gen_alg(m=m.Model(m.p0), func=func_genalg_0000, params=params, ranges=ranges, pop=8, gens=2, jobid=3)
#
# sliderplot(3, 0, 0)

