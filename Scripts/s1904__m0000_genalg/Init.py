import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import InVivo.s190318__wtdata as a
import M as x
import numpy as np
import pickle
import copy
import random

"""
INPUT

"""

direc = x.ddirec + 's1904__m0000_genalg/'
params_base = {'Da': 0.28, 'Dp': 0.15, 'konA': 0.0001, 'koffA': 0.0054, 'konP': 0.0001, 'koffP': 0.0073, 'kAP': 0.0001,
               'kPA': 0.0000001, 'ePneg': 1, 'eAneg': 2, 'xsteps': 500, 'psi': a.svr, 'Tmax': 5000, 'deltat': 0.01,
               'deltal': a.spres, 'radii': a.r_local, 'am_0': a.p6_mem, 'ac_0': a.p6_cyt, 'pm_0': a.p2_mem,
               'pc_0': a.p2_cyt}
params = ['konA', 'konP', 'kAP', 'kPA']
ranges = [[0.0000001, 1], [0.0000001, 1], [0.0000001, 1], [0.0000001, 1]]
popsize = 1000
gen = len(x.direcslist(direc))
mut_rate = 0.2
survival_rate = 0.2
mix_rate = 1  # don't currently have option for this

"""
FIRST GENERATION

"""

if int(gen) == 0:
    os.mkdir(direc + '/g000/')
    x.init_rand(p=params_base, direc=direc + '/g000/', params=params, ranges=ranges, seeds=None, nsims=popsize)

"""
SUBSEQUENT GENERATIONS

"""

if int(gen) != 0:
    os.mkdir(direc + '/g%s' % '{:03d}'.format(int(gen)))

    # Import scores for previous gen
    d = direc + '/g%s' % '{:03d}'.format(int(gen) - 1)
    scores = np.zeros([popsize])
    for p in range(popsize):
        scores[p] = np.loadtxt(d + '/%s/mse.txt' % int(p))

    # Import parameters from previous gen
    valsarray = np.zeros([popsize, len(params)])
    for member in range(popsize):
        p = pickle.load(open(d + '/%s/_params.pkl' % int(member), 'rb'))
        for i, param in enumerate(params):
            valsarray[member, i] = p[param]

    # Create new parameters array
    valsarraytop = valsarray[np.array(scores).argsort()[:int(len(valsarray[:, 0]) * survival_rate)], :]
    valsarraynew = np.zeros([popsize, len(params)])

    for param in range(len(valsarray[0, :])):

        # Mixing
        valsarraynew[:, param] = np.random.choice(valsarraytop[:, param], len(valsarraynew[:, 0]))

        # Mutation
        for simid in range(len(valsarray[:, 0])):
            if random.uniform(0, 1) < mut_rate:
                valsarraynew[simid, param] = random.uniform(ranges[param][0], ranges[param][1])

    # Save
    for j, s in enumerate(valsarraynew):

        # Create copy of original parameter dict
        pcopy = copy.deepcopy(params_base)

        # Make changes
        for i, param in enumerate(params):
            pcopy[param] = s[i]

        # Save parameters file
        x.init_single(pcopy, direc + '/g%s' % '{:03d}'.format(int(gen)), j)
