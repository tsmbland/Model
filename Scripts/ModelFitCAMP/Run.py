import matplotlib

matplotlib.use('Agg')

import copy
import multiprocessing
import itertools
import random
import sys
import os
import numpy as np

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')

from Models.SimplePosFeedback import Model

"""
Funcs

"""


class ModelFit:
    def __init__(self, base_model, kAP_vals, kPA_vals, res, cores, n_sims, direc):
        self.base_model = base_model
        self.kAP_vals = kAP_vals
        self.kPA_vals = kPA_vals
        self.cores = cores
        self.n_sims = n_sims
        self.direc = direc
        self.res = res

    def sim(self, kAP, kPA):
        name = '%s_%s' % (kAP, kPA)
        print(name)

        # Set up model
        model = copy.deepcopy(self.base_model)

        # Add antagonism
        model.kAP = 10 ** kAP
        model.kPA = 10 ** kPA

        # Simulation
        for t in range(int(model.Tmax / model.deltat)):
            model.react()
            model.time = (t + 1) * model.deltat

        # Save error
        error = np.mean((model.am - self.base_model.am_0) ** 2 + (model.pm - self.base_model.pm_0) ** 2)
        np.savetxt(self.direc + '/' + name + '.txt', [error])

    def run(self):
        combinations = np.array(list(itertools.product(self.kAP_vals, self.kPA_vals)))
        filtered_combinations = combinations[self.res.flatten() == 0, :]
        random.seed(1)
        sims_array = filtered_combinations[random.sample(range(len(filtered_combinations)), self.n_sims), :]

        # Run simulations
        pool = multiprocessing.Pool(self.cores)
        pool.starmap(self.sim, iter(sims_array))


"""
Import data

"""

kAP_vals = np.loadtxt(home_direc + '/kAP_vals.txt')
kPA_vals = np.loadtxt(home_direc + '/kPA_vals.txt')
res = np.loadtxt(home_direc + '/Res.txt')

a_mem = np.loadtxt(home_direc + '/a_mem.txt')
a_cyt = 0.17147142857142855
p_mem = np.loadtxt(home_direc + '/p_mem.txt')
p_cyt = 0.06860000000000001

"""
Run

"""

koffA = 0.0092
koffP = 0.0073
psi = 0.10318684114244771
dosP = 0.294005475663175
dosA = 1.05143336288 * dosP

base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
                   koffP=koffP,
                   kposP=12.7357711156 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01,
                   deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

a = ModelFit(base_model=base_model, cores=32, direc=home_direc + '/Res', kAP_vals=kAP_vals, kPA_vals=kPA_vals, res=res,
             n_sims=1000)
a.run()
