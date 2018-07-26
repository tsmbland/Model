"""
180725
m0008
Testing random parameter sets on cluster

"""

from Simulation import *
from Analysis import *
import Models.m0008 as m

p0 = m.Params(Da=0.28, kon_a=1, koff_a=1, ra=1, Dp=0.15, kon_p=1, koff_p=1, kon_p_2=1, kd_f=1, kd_b=1, rp=1, L=67.3,
              xsteps=500, psi=0.174, Tmax=1000, deltat=0.01, starts=[0, 1.56, 0, 0, 0, 1, 0])

params = ['kon_a', 'koff_a', 'ra', 'kon_p', 'koff_p', 'kon_p_2', 'kd_f', 'kd_b', 'rp']
ranges = [[0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10]]

alg_parsim_rand_clust(m=m.Model(p0), params=params, ranges=ranges, nsims=64, jobid=3, subjobid=0, cores=32,
                      node=int(sys.argv[1]))
