import numpy as np
from Funcs import pdeRK, diffusion


class PAR:
    def __init__(self, Da, Dp, konA, koffA, kposA, konP, koffP, kposP, kPA, kAP, eAneg, ePneg, xsteps, psi, Tmax,
                 deltat, L, pA, pP):
        # Species
        self.A = np.zeros([int(xsteps)])
        self.P = np.zeros([int(xsteps)])
        self.time = 0

        # Dosages
        self.pA = pA
        self.pP = pP

        # Diffusion
        self.Da = Da  # input is um2 s-1
        self.Dp = Dp  # um2 s-1

        # Membrane exchange
        self.konA = konA  # um s-1
        self.koffA = koffA  # s-1
        self.konP = konP  # um s-1
        self.koffP = koffP  # s-1

        # Positive feedback
        self.kposA = kposA
        self.kposP = kposP

        # Antagonism
        self.kPA = kPA  # um4 s-1
        self.kAP = kAP  # um2 s-1
        self.eAneg = eAneg
        self.ePneg = ePneg

        # Misc
        self.L = L
        self.xsteps = int(xsteps)
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s
        self.deltax = self.L / xsteps  # um
        self.psi = psi  # um-1

    def dxdt(self, X):
        A = X[0]
        P = X[1]
        ac = self.pA - self.psi * np.mean(A)
        pc = self.pP - self.psi * np.mean(P)
        dA = ((self.konA * ac) - (self.koffA * A) - (self.kAP * (P ** self.ePneg) * A) + (self.kposA * A * ac) + (
            self.Da * diffusion(A, self.deltax)))
        dP = ((self.konP * pc) - (self.koffP * P) - (self.kPA * (A ** self.eAneg) * P) + (self.kposP * P * pc) + (
            self.Dp * diffusion(P, self.deltax)))
        return [dA, dP]

    def initiate(self):
        # Remove antagonism
        kAP, kPA = self.kAP, self.kPA
        self.kAP = 0
        self.kPA = 0

        # Initial equilibration (no antagonism)
        Tmax = self.Tmax / 10
        soln, times = pdeRK(dxdt=self.dxdt, X0=[self.A, self.P], Tmax=Tmax, deltat=self.deltat,
                            t_eval=np.arange(0, Tmax + 0.0001, Tmax))
        self.A = soln[0][-1, :]
        self.P = soln[1][-1, :]

        # Polarise
        self.A *= 2 * np.r_[np.ones([self.xsteps // 2]), np.zeros([self.xsteps // 2])]
        self.P *= 2 * np.r_[np.zeros([self.xsteps // 2]), np.ones([self.xsteps // 2])]

        # Add back antagonism
        self.kAP = kAP
        self.kPA = kPA

    def run(self, save_direc=None, save_gap=None):
        if save_gap is None:
            save_gap = self.Tmax

        # Run
        soln, times = pdeRK(dxdt=self.dxdt, X0=[self.A, self.P], Tmax=self.Tmax, deltat=self.deltat,
                            t_eval=np.arange(0, self.Tmax + 0.0001, save_gap))
        self.A = soln[0][-1, :]
        self.P = soln[1][-1, :]

        # Save
        if save_direc is not None:
            np.savetxt(save_direc + '/A.txt', soln[0])
            np.savetxt(save_direc + '/P.txt', soln[1])
            np.savetxt(save_direc + '/times.txt', times)

    def polarised(self):
        if sum(self.A > self.P) == len(self.A):
            return 3
        elif sum(self.A > self.P) == 0:
            return 1
        else:
            return 2


# from Funcs import Bifurcation2D
# import time
# import matplotlib.pyplot as plt
#
#
# def evaluate1(pA, pP):
#     m = PAR(Da=1, Dp=1, konA=1, koffA=0.3, kposA=0, konP=1,
#             koffP=0.3, kposP=0, kAP=1, kPA=1, ePneg=2, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01,
#             L=50, psi=0.3, pA=pA, pP=pP)
#     m.initiate()
#     m.run()
#     return m.polarised()
#
#
# def evaluate2(kAP, kPA):
#     m = PAR(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
#             koffP=0.01, kposP=0, kAP=kAP, kPA=kPA, ePneg=1, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01,
#             L=50, psi=0.1, pA=1, pP=1)
#     m.initiate()
#     m.run()
#     return m.polarised()
#
#
# a = Bifurcation2D(evaluate1, p1_range=(0, 2), p2_range=(0, 2), log=False, cores=4, resolution0=5,
#                   resolution_step=2, n_iterations=1, direc='_test', parallel=True, crange=[1, 3])
#
# t = time.time()
# a.run()
# print(time.time() - t)
