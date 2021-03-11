import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc + '/../..')

import numpy as np
from polaritymodel import pdeRK, diffusion


class WP:
    def __init__(self, D=0.1, k0=0.067, gamma=1, K=1, delta=1, p0=2.3, xsteps=100, Tmax=1000, deltat=0.01, deltax=0.1):

        # Species
        self.U = np.zeros([int(xsteps)])
        self.time = 0

        # Diffusion
        self.D = D

        # Dosage
        self.p0 = p0

        # Membrane exchange
        self.k0 = k0
        self.gamma = gamma
        self.K = K
        self.delta = delta

        # Misc
        self.xsteps = int(xsteps)
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s
        self.deltax = deltax  # um

    def dxdt(self, X):
        U = X[0]
        V = self.p0 - np.mean(U)
        dU = V * (self.k0 + (self.gamma * (U ** 2)) / ((self.K ** 2) + (U ** 2))) - (self.delta * U) + (
                self.D * diffusion(U, self.deltax))
        return [dU]

    def initiate(self):

        # Initial equilibration
        Tmax = self.Tmax / 10
        soln, time, solns, times = pdeRK(dxdt=self.dxdt, X0=[self.U], Tmax=Tmax, deltat=self.deltat,
                                         t_eval=np.arange(0, Tmax + 0.0001, Tmax))
        self.U = soln[0]

        # Polarise
        self.U *= 2 * np.r_[np.zeros([self.xsteps // 2]), np.ones([self.xsteps // 2])]

    def run(self, save_direc=None, save_gap=None):
        if save_gap is None:
            save_gap = self.Tmax

        # Run
        soln, time, solns, times = pdeRK(dxdt=self.dxdt, X0=[self.U], Tmax=self.Tmax, deltat=self.deltat,
                                         t_eval=np.arange(0, self.Tmax + 0.0001, save_gap))
        self.U = soln[0]

        # Save
        if save_direc is not None:
            np.savetxt(save_direc + '/U.txt', solns[0])
            np.savetxt(save_direc + '/times.txt', times)
