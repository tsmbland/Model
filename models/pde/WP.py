import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc + '/../..')

import numpy as np
from parmodel import pdeRK, diffusion


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

# from Funcs import Bifurcation2D
# import time
# import matplotlib.pyplot as plt
#
#
# def evaluate1(k0, p0):
#     m = WP(k0=k0, p0=p0)
#     m.initiate()
#     m.run()
#     return m.state()
#
#
# a = Bifurcation2D(evaluate1, p1_range=(0.0001, 0.12), p2_range=(0.0001, 4), log=False, cores=4, resolution0=5,
#                   resolution_step=2, n_iterations=2, direc='_test', parallel=True, crange=[1, 2])
#
# t = time.time()
# a.run()
# print(time.time() - t)


# m = WP()
# m.initiate()
# m.run(save_direc='_test', save_gap=10)
#
# from matplotlib.widgets import Slider
#
#
# def animate(direc):
#     """
#     direc: directory to results
#
#     """
#
#     times = np.loadtxt(direc + '/times.txt')
#     U = np.loadtxt(direc + '/U.txt')
#
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     plt.subplots_adjust(bottom=0.25, wspace=0.5)
#     axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
#     sframe = Slider(axframe, 'Time (s)', 0, times[-1], valinit=0, valfmt='%d')
#
#     def update(i):
#         ax.clear()
#         tpoint = np.argmin(abs(times - int(i)))
#         u = U[tpoint, :]
#         ax.plot(u, c='tab:red')
#         # ax.set_ylim(bottom=0)
#         ax.set_ylabel('Cortical concentration (a.u.)')
#         ax.set_xlabel('Position (μm)')
#
#     sframe.on_changed(update)
#     plt.show()
#
#
# animate('_test')
