"""
180725
Testing m0008

"""

import M as x
import Models.m0008 as m0008

x.datadirec = '../../ModelData'

######### Single simulation <- good

p = m0008.Params(Da=1, kon_a=1, koff_a=1, ra=1, Dp=1, kon_p=1, koff_p=1, kon_p_2=5, kd_f=2, kd_b=1, rp=1, L=50, xsteps=500,
            psi=0.3, Tmax=100, deltat=0.01, starts=[0, 1, 0, 0, 0, 1, 0])

p.Tmax = 100

# p.kon_a = 0
# p.koff_a = 0
#
# p.kon_p = 0
# p.kon_p_2 = 0
# p.koff_p = 0
# p.kd_b = 0
# p.kd_f = 0

model = m0008.Model(p)
x.alg_singlesim(model, 5, 0, 0, compression=0)
x.sliderplot(5, 0, 0)
res = x.loaddata(5, 0, 0)
x.plt.plot(res.atot)
x.plt.plot(res.ptot)
x.plt.show()



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
