import matplotlib

# matplotlib.use('Agg')

import numpy as np
import copy
from Funcs import Bifurcation2D
from Models.SimplePosFeedback import Model

"""
Parameter run

"""


def func(kAP, kPA):
    # Set up model
    model = copy.deepcopy(base_model)

    # Remove antagonism
    model.kAP = 0
    model.kPA = 0

    # Initial equilibration (no antagonism)
    for t in range(10000):
        model.react()

    # Polarise
    model.am *= 2 * np.r_[np.ones([model.xsteps // 2]), np.zeros([model.xsteps // 2])]
    model.pm *= 2 * np.r_[np.zeros([model.xsteps // 2]), np.ones([model.xsteps // 2])]

    # Add antagonism
    model.kAP = kAP
    model.kPA = kPA

    # Simulation
    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat

    # Return state
    if sum(model.am > model.pm) == len(model.am):
        return 1
    elif sum(model.am > model.pm) == 0:
        return -1
    else:
        return 0


koffA = 0.0092
koffP = 0.0073
psi = 0.10318684114244771
dosP = 0.294005475663175
dosA = 1.05143336288 * dosP

base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP,
                   koffP=koffP, kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01,
                   deltax=0.5,
                   psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

a = Bifurcation2D(func, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=8,
                  resolution0=5, resolution_step=2, n_iterations=10, direc='_temp')
a.run()
