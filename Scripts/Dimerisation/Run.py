import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')

from Models.pPAR_Dimerisation_NonSpatial import Model
import numpy as np
import copy
import itertools
import multiprocessing

print(sys.argv[1])

save_direc = home_direc + '/../../../../ModelData/Dimerisation/'

"""
Functions

"""


class Sweep:
    def __init__(self, base_model, p1, p1_vals, p2, p2_vals, cores, direc):
        self.base_model = base_model
        self.p1 = p1
        self.p2 = p2
        self.p1_vals = p1_vals
        self.p2_vals = p2_vals
        self.cores = cores
        self.direc = direc

    def run(self, p1_val, p2_val):
        name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))

        # Test if simulation has already been performed
        if not os.path.exists(self.direc + '/' + name):
            print(self.direc + '/' + name)

            # Set up model
            model = copy.deepcopy(self.base_model)
            setattr(model, self.p1, 10 ** p1_val)
            setattr(model, self.p2, 10 ** p2_val)
            dosages = np.linspace(0, 1, 100)
            model.pc1 = dosages
            model.pc2 = 0 * dosages
            model.pm1 = 0 * dosages
            model.pm2s = 0 * dosages
            model.pm2d = 0 * dosages

            # Simulation
            for t in range(int(model.Tmax / model.deltat)):
                model.react()
                model.time = (t + 1) * model.deltat
            model.pool()

            # Save results
            os.mkdir(self.direc + '/' + name)
            model.save(self.direc + '/' + name + '/')

    def sweep(self):

        # Run
        sims_array = np.array(list(itertools.product(self.p1_vals, self.p2_vals)))
        pool = multiprocessing.Pool(self.cores)
        pool.starmap(self.run, iter(sims_array))


class Sweep2:
    """
    For sweeping kd_f and kd_b together


    """

    def __init__(self, base_model, vals, cores, direc):
        self.base_model = base_model
        self.vals = vals
        self.cores = cores
        self.direc = direc

    def run(self, val):
        name = float('%.4g' % val)

        # Test if simulation has already been performed
        if not os.path.exists(self.direc + '/' + name):
            print(self.direc + '/' + name)

            # Set up model
            model = copy.deepcopy(self.base_model)
            setattr(model, 'kd_f', 10 ** val)
            setattr(model, 'kd_b', 10 ** val)
            dosages = np.linspace(0, 1, 100)
            model.pc1 = dosages
            model.pc2 = 0 * dosages
            model.pm1 = 0 * dosages
            model.pm2s = 0 * dosages
            model.pm2d = 0 * dosages

            # Simulation
            for t in range(int(model.Tmax / model.deltat)):
                model.react()
                model.time = (t + 1) * model.deltat
            model.pool()

            # Save results
            os.mkdir(self.direc + '/' + name)
            model.save(self.direc + '/' + name + '/')

    def sweep(self):

        # Run
        pool = multiprocessing.Pool(self.cores)
        pool.starmap(self.run, iter(self.vals))


"""
Run

"""

if int(sys.argv[1]) == 0:
    base_model = Model(kon_p=1, kon_p2=1, koff_p=1, kd_f=1, kd_b=1, psi=0.1, Tmax=1000, deltat=0.001,
                       pm1_0=0, pm2s_0=0, pm2d_0=0, pc1_0=0, pc2_0=0)

    a = Sweep(base_model=base_model, p1='kd_f', p2='kd_b', p1_vals=np.linspace(-2, 2, 20),
              p2_vals=np.linspace(-2, 2, 20), cores=32, direc=save_direc + '/0')
    a.sweep()

if int(sys.argv[1]) == 1:
    base_model = Model(kon_p=1, kon_p2=10, koff_p=1, kd_f=1, kd_b=1, psi=0.1, Tmax=1000, deltat=0.001,
                       pm1_0=0, pm2s_0=0, pm2d_0=0, pc1_0=0, pc2_0=0)

    a = Sweep(base_model=base_model, p1='kd_f', p2='kd_b', p1_vals=np.linspace(-2, 2, 20),
              p2_vals=np.linspace(-2, 2, 20), cores=32, direc=save_direc + '/1')
    a.sweep()

if int(sys.argv[1]) == 2:
    base_model = Model(kon_p=1, kon_p2=100, koff_p=1, kd_f=1, kd_b=1, psi=0.1, Tmax=1000, deltat=0.001,
                       pm1_0=0, pm2s_0=0, pm2d_0=0, pc1_0=0, pc2_0=0)

    a = Sweep(base_model=base_model, p1='kd_f', p2='kd_b', p1_vals=np.linspace(-2, 2, 20),
              p2_vals=np.linspace(-2, 2, 20), cores=32, direc=save_direc + '/2')
    a.sweep()

if int(sys.argv[1]) == 3:
    base_model = Model(kon_p=1, kon_p2=10, koff_p=0.1, kd_f=0.1, kd_b=0.1, psi=0.1, Tmax=1000, deltat=0.001,
                       pm1_0=0, pm2s_0=0, pm2d_0=0, pc1_0=0, pc2_0=0)

    a = Sweep(base_model=base_model, p1='kon_p', p2='kon_p2', p1_vals=np.linspace(-2, 2, 20),
              p2_vals=np.linspace(-2, 2, 20), cores=32, direc=save_direc + '/3')
    a.sweep()

if int(sys.argv[1]) == 4:
    base_model = Model(kon_p=1, kon_p2=10, koff_p=1, kd_f=1, kd_b=1, psi=0.1, Tmax=1000, deltat=0.001,
                       pm1_0=0, pm2s_0=0, pm2d_0=0, pc1_0=0, pc2_0=0)

    a = Sweep(base_model=base_model, p1='kon_p', p2='kon_p2', p1_vals=np.linspace(-2, 2, 20),
              p2_vals=np.linspace(-2, 2, 20), cores=32, direc=save_direc + '/4')
    a.sweep()

if int(sys.argv[1]) == 5:
    base_model = Model(kon_p=1, kon_p2=10, koff_p=10, kd_f=10, kd_b=10, psi=0.1, Tmax=1000, deltat=0.001,
                       pm1_0=0, pm2s_0=0, pm2d_0=0, pc1_0=0, pc2_0=0)

    a = Sweep(base_model=base_model, p1='kon_p', p2='kon_p2', p1_vals=np.linspace(-2, 2, 20),
              p2_vals=np.linspace(-2, 2, 20), cores=32, direc=save_direc + '/5')
    a.sweep()
