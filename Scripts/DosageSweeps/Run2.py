import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/DosageSweeps2/'

from Models.Simple import Model as M1
from Models.WavePinning import Model as M2

import numpy as np
import itertools
import multiprocessing
import copy

print(sys.argv[1])

"""
Dosage sweep

"""


class Sweep:
    def __init__(self, base_model, param, vals, cores, direc):
        self.base_model = base_model
        self.param = param
        self.vals = vals
        self.cores = cores
        self.direc = direc

    def run(self, val):
        name = '%s' % float('%.4g' % val)

        # Test if simulation has already been performed
        if not os.path.exists(self.direc + '/' + name):
            print(self.direc + '/' + name)

            # Set up model
            model = copy.deepcopy(self.base_model)
            setattr(model, self.param, val)

            # Remove antagonism
            kAP, kPA = model.kAP, model.kPA
            model.kAP = 0
            model.kPA = 0

            # Initial equilibration (no antagonism)
            for t in range(10000):
                model.react()

            # Add antagonism
            model.kAP = kAP
            model.kPA = kPA

            # Polarise
            model.am *= 2 * np.r_[np.ones([250]), np.zeros([250])]
            model.pm *= 2 * np.r_[np.zeros([250]), np.ones([250])]

            # Simulation
            for t in range(int(model.Tmax / model.deltat)):
                model.react()
                model.time = (t + 1) * model.deltat

            # Save results
            os.mkdir(self.direc + '/' + name)
            model.save(self.direc + '/' + name + '/')

    def sweep(self):

        # Run
        pool = multiprocessing.Pool(self.cores)
        pool.map(self.run, self.vals)


"""
1: Generic PAR (Goehring 2011)

"""

if int(sys.argv[1]) == 1:
    base_model = M1(Da=1, Dp=1, konA=1, koffA=0.3, konP=1, koffP=0.3, kAP=1.5, kPA=1.5, ePneg=2,
                    eAneg=2, xsteps=500, psi=0.3, Tmax=1000, deltat=0.001, deltax=0.1, am_0=np.zeros([500]),
                    ac_0=1, pm_0=np.zeros([500]), pc_0=1)
    s = Sweep(base_model=base_model, param='ac', vals=np.linspace(0, 2, 100), cores=32, direc=save_direc + '/1')
    s.sweep()

if int(sys.argv[1]) == 2:
    base_model = M1(Da=1, Dp=1, konA=1, koffA=0.3, konP=1, koffP=0.3, kAP=1.5, kPA=1.5, ePneg=2,
                    eAneg=2, xsteps=500, psi=0.3, Tmax=1000, deltat=0.001, deltax=0.1, am_0=np.zeros([500]),
                    ac_0=1, pm_0=np.zeros([500]), pc_0=1)
    s = Sweep(base_model=base_model, param='pc', vals=np.linspace(0, 2, 100), cores=32, direc=save_direc + '/2')
    s.sweep()

"""
2: Generic PAR (Hubatsch 2019)

"""

if int(sys.argv[1]) == 3:
    base_model = M1(Da=0.1, Dp=0.1, konA=0.006, koffA=0.005, konP=0.006, koffP=0.005, kAP=1, kPA=1, ePneg=2,
                    eAneg=2, xsteps=500, psi=0.127, Tmax=1000, deltat=0.01, deltax=0.06, am_0=np.zeros([500]),
                    ac_0=1, pm_0=np.zeros([500]), pc_0=1)
    s = Sweep(base_model=base_model, param='ac', vals=np.linspace(0, 2, 100), cores=32, direc=save_direc + '/3')
    s.sweep()

if int(sys.argv[1]) == 4:
    base_model = M1(Da=0.1, Dp=0.1, konA=0.006, koffA=0.005, konP=0.006, koffP=0.005, kAP=1, kPA=1, ePneg=2,
                    eAneg=2, xsteps=500, psi=0.127, Tmax=1000, deltat=0.01, deltax=0.06, am_0=np.zeros([500]),
                    ac_0=1, pm_0=np.zeros([500]), pc_0=1)
    s = Sweep(base_model=base_model, param='pc', vals=np.linspace(0, 2, 100), cores=32, direc=save_direc + '/4')
    s.sweep()

"""
3: Generic PAR (new parameters)

"""

if int(sys.argv[1]) == 5:
    base_model = M1(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, konP=0.1, koffP=0.01, kAP=0.1, kPA=0.1, ePneg=2,
                    eAneg=2, xsteps=500, psi=0.1, Tmax=1000, deltat=0.01, deltax=0.1, am_0=np.zeros([500]),
                    ac_0=1, pm_0=np.zeros([500]), pc_0=1)
    s = Sweep(base_model=base_model, param='ac', vals=np.linspace(0, 2, 100), cores=32, direc=save_direc + '/5')
    s.sweep()

if int(sys.argv[1]) == 6:
    base_model = M1(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, konP=0.1, koffP=0.01, kAP=0.1, kPA=0.1, ePneg=2,
                    eAneg=2, xsteps=500, psi=0.1, Tmax=1000, deltat=0.01, deltax=0.1, am_0=np.zeros([500]),
                    ac_0=1, pm_0=np.zeros([500]), pc_0=1)
    s = Sweep(base_model=base_model, param='pc', vals=np.linspace(0, 2, 100), cores=32, direc=save_direc + '/6')
    s.sweep()

"""
4: Goehring 2011

"""

if int(sys.argv[1]) == 7:
    base_model = M1(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, ePneg=1,
                    eAneg=2, xsteps=500, psi=0.127, Tmax=1000, deltat=0.01, deltax=0.135, am_0=np.zeros([500]),
                    ac_0=1.56, pm_0=np.zeros([500]), pc_0=1)
    s = Sweep(base_model=base_model, param='ac', vals=np.linspace(0, 3, 100), cores=32, direc=save_direc + '/7')
    s.sweep()

if int(sys.argv[1]) == 8:
    base_model = M1(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, ePneg=1,
                    eAneg=2, xsteps=500, psi=0.127, Tmax=1000, deltat=0.01, deltax=0.135, am_0=np.zeros([500]),
                    ac_0=1.56, pm_0=np.zeros([500]), pc_0=1)
    s = Sweep(base_model=base_model, param='pc', vals=np.linspace(0, 2, 100), cores=32, direc=save_direc + '/8')
    s.sweep()
