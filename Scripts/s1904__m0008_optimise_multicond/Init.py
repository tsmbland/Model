import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import InVivo.s190318__wtdata as a
import M as x
import pandas as pd
import numpy as np

"""
m0008 random params

"""

direc = x.ddirec + 's1904__m0008_optimise_multicond/'
params_base = {'Da': 0.28, 'kon_a': 0, 'koff_a': 0.0054, 'ra': 0, 'Dp': 0.15, 'kon_p': 0, 'koff_p': 0.0073,
               'kon_p_2': 0, 'kd_f': 0, 'kd_b': 0, 'rp': 0, 'xsteps': 500, 'psi': a.svr, 'Tmax': 1000, 'deltat': 0.01,
               'deltal': a.spres, 'radii': a.r_local, 'am_0': a.p6_mem, 'ac_0': a.p6_cyt, 'pm1_0': a.p2_mem,
               'pm2s_0': 0 * a.p2_mem, 'pm2d_0': 0 * a.p2_mem, 'pc1_0': a.p2_cyt, 'pc2_0': 0 * a.p2_cyt}
params = ['kon_a', 'ra', 'kon_p', 'kon_p_2', 'kd_f', 'kd_b', 'rp']
ranges = [[0.001, 0.1], [0.0001, 1], [0.0001, 0.1], [0.01, 10], [0.01, 10], [0.01, 10], [0.0001, 1]]

# x.init_rand(p=params_base, direc=direc, params=params, ranges=ranges, seeds=None,
#             nsims=100000)

valsarray = x.params_array_rand_lhs(params, ranges, 100000)
d = dict(zip(params, valsarray.T))
pd.DataFrame(d).to_csv(direc + 'db.csv')
