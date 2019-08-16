import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import InVivo.s190318__wtdata as a
import M as x
import numpy as np

"""
INPUT

"""

direc = x.ddirec + 's1904__m0000_robustness/'

params_base = {'Da': 0.28, 'Dp': 0.15, 'konA': 0.01377691, 'koffA': 0.0054, 'konP': 0.04790167, 'koffP': 0.0073,
               'kAP': 0.0001,
               'kPA': 0.0000001, 'ePneg': 1, 'eAneg': 1, 'xsteps': 500, 'psi': a.svr, 'Tmax': 1000, 'deltat': 0.01,
               'deltal': a.spres, 'radii': a.r_local, 'am_0': a.p6_mem, 'ac_0': a.p6_cyt, 'pm_0': a.p2_mem,
               'pc_0': a.p2_cyt}
params = ['eAneg', 'kPA', 'kAP']
ranges = [[1, 2, 3], np.linspace(0, 0.005, 100), np.linspace(0, 0.001, 100)]

x.init_multi(p=params_base, direc=direc, params=params, vals=ranges)
