import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import InVivo.s190318__wtdata as a
import M as x

"""
m0008 random params

"""

direc = x.ddirec + 's1903__m0000_optimise'

params_base = {'ac_on': 0, 'am_off': 0.0054, 'sc_on': 0, 'sm_off': 0, 'pc_on': 0, 'pm_off': 0, 'k_dim_f': 0,
               'k_dim_b': 0, 'pc_on2': 0, 'sc_on2': 0, 'pm_off2': 0, 'sm_off2': 0, 'r_am_p': 0, 'r_sm_a': 0,
               'r_pm_a': 0, 'r_spm10_a': 0, 'r_spm11_a': 0, 'r_spm01_a': 0, 'd_am': 0.28, 'd_sm': 0.15, 'd_pm': 0.15,
               'd_spm10': 0.15, 'd_spm11': 0.15, 'd_spm01': 0.15, 'xsteps': 500, 'psi': a.svr, 'Tmax': 100,
               'deltat': 0.01,
               'deltal': a.spres, 'radii': a.r_local, 'am_0': 0, 'ac_0': 0, 'pm_0': 0, 'pc_0': 0, 'sm_0': 0, 'sc_0': 0,
               'spm10_0': 0, 'spm11_0': 0, 'spm01_0': 0, 'spc_0': 0}

os.mkdir(direc)

# x.init_rand(p=params_base, direc=direc, params=['kon_a', 'ra', 'kon_p', 'kon_p_2', 'kd_f', 'kd_b', 'r_p'],
#             ranges=[[0.0000001, 1], [0.0000001, 1], [0.0000001, 1], [0.0000001, 1], [0.0000001, 1], [0.0000001, 1],
#                     [0.0000001, 1]], seeds=None, nsims=50)
