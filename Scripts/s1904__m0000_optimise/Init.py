import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import InVivo.s190318__wtdata as a
import M as x

"""
INPUT

"""

direc = x.ddirec + 's1904__m0000_optimise_B/'
params_base = {'Da': 0.28, 'Dp': 0.15, 'konA': 0.0001, 'koffA': 0.0054, 'konP': 0.0001, 'koffP': 0.0073, 'kAP': 0.0001,
               'kPA': 0.0000001, 'ePneg': 1, 'eAneg': 2, 'xsteps': 500, 'psi': a.svr, 'Tmax': 1000, 'deltat': 0.01,
               'deltal': a.spres, 'radii': a.r_local, 'am_0': a.p6_mem, 'ac_0': a.p6_cyt, 'pm_0': a.p2_mem,
               'pc_0': a.p2_cyt}
params = ['konA', 'konP', 'kAP', 'kPA']
ranges = [[0.001, 0.1], [0.001, 0.1], [0.0001, 0.1], [0.0001, 0.1]]
popsize = 1000
gen = len(x.direcslist(direc))

"""
GENERATIONS

"""

os.mkdir(direc + '/g%s' % '{:03d}'.format(int(gen)))
x.init_rand(p=params_base, direc=direc + '/g%s' % '{:03d}'.format(int(gen)), params=params, ranges=ranges, seeds=None,
            nsims=popsize)
