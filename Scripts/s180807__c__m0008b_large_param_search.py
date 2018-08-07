"""

"""

import M as x
import Models.m0008b as m0008b
import sys

x.datadirec = '../working/Tom/ModelData'

# Base parameter set
p0 = m0008b.Params(Da=0.25, kon_a=0.0085, koff_a=0.0054, ra=1, Dp=0.15, kon_p=0.0474, koff_p=0.0073, kon_p_2=0.0474,
                   kd_f=1, kd_b=1, rp=1, L=67.3,
                   xsteps=500, psi=0.174, Tmax=1000, deltat=0.01, starts=[0, 1.56, 0, 0, 0, 1, 0])

# Free parameters
params = ['ra', 'kon_p', 'kon_p_2', 'kd_f', 'kd_b', 'rp']
# ranges = [[0.001, 10], [0.001, 10], [0.001, 10], [0.001, 10], [0.001, 10], [0.001, 10]]
ranges = [[0.001, 0.1], [0.0001, 0.01], [0.1, 10], [0.01, 1], [0.01, 1], [0.0001, 0.01]]

# # Standard simulation + analysis
# x.alg_parsim_rand_clust(m=m0008b.Model(p0), params=params, ranges=ranges, nsims=960, jobid=11, subjobid=0, cores=32,
#                         compression=0, node=int(sys.argv[1]))
# x.batch_analysis(11, 0, funcs=x.all_analysis)
# x.save_scores_batch(11, 0)

# # Retesting in PAR-1 RNAi mode + analysis
# x.alg_parsim_clust_duplicate(m=m0008b.Model(p0), changes={'ra': 0}, oldjobid=11, oldsubjobid=0, newjobid=11,
#                              newsubjobid=1,
#                              cores=32, node=int(sys.argv[1]), compression=0)
# x.batch_analysis(11, 1, funcs=x.all_analysis)
# x.save_scores_batch(11, 1)

# Retesting in PKC RNAi mode + analysis
x.alg_parsim_clust_duplicate(m=m0008b.Model(p0), changes={'ac_0': 0}, oldjobid=11, oldsubjobid=0, newjobid=11,
                             newsubjobid=2,
                             cores=32, node=int(sys.argv[1]), compression=0)
x.batch_analysis(11, 2, funcs=x.all_analysis)
x.save_scores_batch(11, 2)


