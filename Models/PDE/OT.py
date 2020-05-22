import numpy as np
from Funcs import pdeRK, diffusion


class OT:
    def __init__(self, D, a1, a2, s, p0, xsteps=100, Tmax=1000, deltat=0.01, deltax=0.1):

        # Species
        self.U = np.zeros([int(xsteps)])
        self.time = 0

        # Diffusion
        self.D = D

        # Dosage
        self.p0 = p0

        # Membrane exchange
        self.a1 = a1
        self.a2 = a2
        self.s = s

        # Misc
        self.xsteps = int(xsteps)
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s
        self.deltax = deltax  # um

    def dxdt(self, X):
        U = X[0]
        V = self.p0 - np.mean(U)
        dUdt = self.a1 * (V - (U + V) / ((self.a2 * self.s * (U + V) + 1) ** 2)) + (self.D * diffusion(U, self.deltax))
        return [dUdt]

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

