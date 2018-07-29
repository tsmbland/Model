"""
180725
m0000
Testing genetic algorithm on cluster
Not run yet

"""

import M as x
import Models.m0000 as m0000
import sys

x.datadirec = '../working/Tom/ModelData'

params = ['pA', 'pP']
ranges = [[0, 2], [0, 2]]
x.gen_alg_clust(m=m0000.Model(m0000.p0), func=x.func_genalg_0000, params=params, ranges=ranges, jobid=3, cores=32, innergens=3,
                node=int(sys.argv[1]), nodes=int(sys.argv[2]))
