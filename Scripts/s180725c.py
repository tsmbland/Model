"""
180725
Testing m0008
Using genetic algorithm to find parameter set that polarises

"""

from M import *
import Models.m0008 as m0008
import Models.m0000 as m0000

alg_singlesim(m0000.Model(m0000.p0), 9999, 0, 0)

p0 = m0008.Params(Da=1, kon_a=1, koff_a=1, ra=1, Dp=1, kon_p=1, koff_p=1, kon_p_2=5, kd_f=2, kd_b=1, rp=1, L=50,
                  xsteps=500, psi=0.3, Tmax=100, deltat=0.01, starts=[0, 1, 0, 0, 0, 1, 0])

params = ['kon_a', 'koff_a', 'ra', 'kon_p', 'koff_p', 'kon_p_2', 'kd_f', 'kd_b', 'rp']
ranges = [[0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10]]

gen_alg(m=m0008.Model(m0008.p0), func=func_genalg_0001, params=params, ranges=ranges, pop=8, gens=1, jobid=6)

genalg_plot1(6, [9999, 0, 0])

