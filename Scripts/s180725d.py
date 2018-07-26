"""
180725
m0000
Testing genetic algorithm on cluster
Not run yet

"""

from Simulation import *
from Analysis import *
import Models.m0000 as m

params = ['pA', 'pP']
ranges = [[0, 2], [0, 2]]
gen_alg_clust(m=m.Model(m.p0), func=func_genalg_0000, params=params, ranges=ranges, jobid=3, cores=32, innergens=3,
              node=int(sys.argv[1]), nodes=int(sys.argv[2]))
