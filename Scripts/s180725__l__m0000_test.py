"""
180725
Testing m0000

"""

import M as x
import Models.m0000 as m0000

x.datadirec = '../../ModelData'


########## Single simulaton <- good

p0 = m0000.Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, ePneg=1, eAneg=2,
            pA=1.56, pP=1, L=67.3, xsteps=500, psi=0.174, Tmax=5000, deltat=0.1, Aeqmin=0, Aeqmax=0.5, Peqmin=0.5,
            Peqmax=1)

model = m0000.Model(m0000.p0)
x.alg_singlesim(model, 2, 0, 0, funcs=x.all_analysis)
# x.plot_singlesim(2, 0, 0)
x.parplot_norm(2, 0, 0)
# x.sliderplot(2, 0, 0)
# res = x.loaddata(2, 0, 0)
# print(res.scores)


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

# params = ['pA', 'pP']
# ranges = [[0, 2], [0, 2]]
# alg_parsim_rand(m0000.Model(m0000.p0), params=params, ranges=ranges, nsims=4, jobid=2, subjobid=2, compression=2, funcs=[mse, asi_a, asi_p])
# print(loaddata(2, 2, 0).params.pA)
# print(loaddata(2, 2, 0).params.pP)
# print(loaddata(2, 2, 0).scores)


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

# from pprint import pprint
# a = loaddata(2, 2, 0)
# pprint(vars(a))


