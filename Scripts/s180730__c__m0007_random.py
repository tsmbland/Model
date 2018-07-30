"""
180730
m0007
Testing random parameter sets on cluster

"""

import M as x
import Models.m0007 as m0007
import sys

x.datadirec = '../working/Tom/ModelData'

p0 = m0007.Params(Da=0.28, konan=1, koffan=1, konam=1, koffam=1, konap=1, koffap=1, kphosan=1, kphosam=1, kdephosap=1,
                  kdephosam=1, Dp=0.15,
                  konpn=1, koffpn=1, konpm=1, koffpm=1, konpp=1, koffpp=1, kphospn=1, kphospm=1, kdephospp=1,
                  kdephospm=1,
                  L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.01, starts=[0, 0, 0, 1.56, 0, 0, 0, 0, 0, 1, 0, 0])

params = ['konan', 'koffan', 'konam', 'koffam', 'konap', 'koffap', 'kphosan', 'kphosam', 'kdephosap', 'kdephosam',
          'konpn', 'koffpn', 'konpm', 'koffpm', 'konpp', 'koffpp', 'kphospn', 'kphospm', 'kdephospp', 'kdephospm']
ranges = [[0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1],
          [0.001, 1],
          [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1],
          [0.001, 1]]

x.alg_parsim_rand_clust(m=m0007.Model(p0), params=params, ranges=ranges, nsims=960, jobid=6, subjobid=0, cores=32,
                        compression=2, node=int(sys.argv[1]), funcs=[x.mse, x.asi_a, x.asi_p])

x.save_scores_batch(6, 0)
