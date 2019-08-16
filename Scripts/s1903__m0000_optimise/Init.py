import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import InVivo.s190318__wtdata as a
import M as x

"""
m0000 random params

"""

direc = x.ddirec + 's1903__m0000_optimise'

params_base = {'Da': 0.28, 'Dp': 0.15, 'konA': 0.0001, 'koffA': 0.0054, 'konP': 0.0001, 'koffP': 0.0073, 'kAP': 0.0001,
               'kPA': 0.0000001, 'ePneg': 1, 'eAneg': 2, 'xsteps': 500, 'psi': a.svr, 'Tmax': 10000, 'deltat': 0.01,
               'deltal': a.spres, 'radii': a.r_local, 'am_0': a.p6_mem, 'ac_0': a.p6_cyt, 'pm_0': a.p2_mem,
               'pc_0': a.p2_cyt}

os.mkdir(direc)

x.init_rand(p=params_base, direc=direc, params=['konA', 'konP', 'kAP', 'kPA'],
            ranges=[[0.0000001, 1], [0.0000001, 1], [0.0000001, 1], [0.0000001, 1]], seeds=None,
            nsims=10000)
