"""
180725
m0009
Testing random parameter sets on cluster

"""

import M as x
import Models.m0009 as m0009
import sys

x.datadirec = '../working/Tom/ModelData'

p0 = m0009.Params(pA=1.56, Da=0.28, konA=1, koffA=1, kAP=1, eAP=1, pP=1, Dp=0.15, konP=1, koffP=1, kPA=1, ePA=2, pS=1,
                  Ds=0, konS=1, koffS=1, kSA=1, eSA=2, L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.01)

params = ['konA', 'koffA', 'kAP', 'koffP', 'kPA', 'pS', 'konS', 'koffS', 'kSA']
ranges = [[0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10], [0.1, 10]]

x.alg_parsim_rand_clust(m=m0009.Model(p0), params=params, ranges=ranges, nsims=960, jobid=5, subjobid=0, cores=32,
                        compression=2, node=int(sys.argv[1]), funcs=[x.mse, x.asi_a, x.asi_p])

x.save_scores_batch(5, 0)
