import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import InVivo.s1905__par2_ng as a
import M as x
import pandas as pd
import pickle

"""
m0000 random params

"""

direc = x.ddirec + 's1905__m0001_optimise'

params_base = {'Da': 0.28, 'Dp': 0.15, 'konA': 0.0, 'koffA': 0.0054, 'konP': 0.0951952, 'koffP': 0.0073,
               'kposP': 0.06296296,
               'kAP': 0.0, 'kPA': 0.0, 'ePneg': 1, 'eAneg': 2, 'xsteps': 500, 'psi': a.svr, 'Tmax': 10000,
               'deltat': 0.01, 'deltal': a.spres, 'radii': a.r_local, 'am_0': a.p6_mem, 'ac_0': a.p6_cyt,
               'pm_0': a.p2_mem, 'pc_0': a.p2_cyt}

params = ['konA', 'kAP', 'kPA']
ranges = [[0.001, 0.1], [0.0001, 0.1], [0.0001, 0.1]]
popsize = 10000

os.mkdir(direc)

valsarray = x.params_array_rand_lhs(params, ranges, popsize)
d = dict(zip(params, valsarray.T))
pd.DataFrame(d).to_csv('db.csv')
pickle.dump(params_base, open('%s/_base.pkl' % direc, 'wb'))
