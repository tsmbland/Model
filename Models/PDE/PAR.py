import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc + '/../..')

import numpy as np
from Funcs import pdeRK, diffusion
from Models.ODE.PAR import PAR as ODE
from scipy.integrate import odeint


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

        # Solve ODE, no antagonism
        o = ODE(konA=self.konA, koffA=self.koffA, kposA=self.kposA, konP=self.konP, koffP=self.koffP, kposP=self.kposP,
                ePneg=self.ePneg, eAneg=self.eAneg, psi=self.psi, pA=self.pA, pP=self.pP, kAP=0, kPA=0)
        soln = odeint(o.dxdt, (0, 0), t=np.linspace(0, 10000, 100000))[-1]

        self.A = soln[0]
        self.P = soln[1]

        # Polarise
        self.A *= 2 * np.r_[np.ones([self.xsteps // 2]), np.zeros([self.xsteps // 2])]
        self.P *= 2 * np.r_[np.zeros([self.xsteps // 2]), np.ones([self.xsteps // 2])]

    def initiate2(self):

        # Solve ODE
        o = ODE(konA=self.konA, koffA=self.koffA, kposA=self.kposA, konP=self.konP, koffP=self.koffP, kposP=self.kposP,
                ePneg=self.ePneg, eAneg=self.eAneg, psi=self.psi, pA=self.pA, pP=self.pP, kAP=self.kAP, kPA=self.kPA)
        soln = odeint(o.dxdt, (o.pA / o.psi, 0), t=np.linspace(0, 10000, 100000))[-1]

        # Set concentrations
        self.A[:] = soln[0]
        self.P[:] = soln[1]

        # Polarise
        self.A *= np.linspace(1.01, 0.99, self.xsteps)
        self.P *= np.linspace(0.99, 1.01, self.xsteps)

    def initiate3(self, asi):

        # Solve ODE
        o = ODE(konA=self.konA, koffA=self.koffA, kposA=self.kposA, konP=self.konP, koffP=self.koffP, kposP=self.kposP,
                ePneg=self.ePneg, eAneg=self.eAneg, psi=self.psi, pA=self.pA, pP=self.pP, kAP=self.kAP, kPA=self.kPA)
        sol = odeint(o.dxdt, (o.pA / o.psi, 0), t=np.linspace(0, 10000, 100000))[-1]

        # Calculate asymmetry
        x = sol[0]
        y = (asi + 0.5) / (0.5 - asi)
        a = x * y / (1 + y)
        p = x - a

        # Set concentrations
        self.A[:] = np.r_[a * np.ones([self.xsteps // 2]), p * np.ones([self.xsteps // 2])]
        self.P[:] = sol[1]

    def run(self, save_direc=None, save_gap=None, kill_uni=False, kill_stab=False):
        """

        :param save_direc: if given, will save A and P distributions over time according to save_gap
        :param save_gap: gap in model time between save points
        :param kill_uni: terminate once polarity is lost. Generally can assume models never regain polarity once lost
        :param kill_stab: terminate when patterns are stable
        :return:
        """
        if save_gap is None:
            save_gap = self.Tmax

        # Kill when uniform
        if kill_uni:
            def killfunc(X):
                if sum(X[0] > X[1]) == len(X[0]) or sum(X[0] > X[1]) == 0:
                    return True
                return False
        else:
            killfunc = None

        # Run
        soln, time, solns, times = pdeRK(dxdt=self.dxdt, X0=[self.A, self.P], Tmax=self.Tmax, deltat=self.deltat,
                                         t_eval=np.arange(0, self.Tmax + 0.0001, save_gap), killfunc=killfunc,
                                         stabilitycheck=kill_stab)
        self.A = soln[0]
        self.P = soln[1]

        # Save
        if save_direc is not None:
            np.savetxt(save_direc + '/A.txt', solns[0])
            np.savetxt(save_direc + '/P.txt', solns[1])
            np.savetxt(save_direc + '/times.txt', times)

# import matplotlib.pyplot as plt
# from Funcs import animatePAR
#
# kon0 = -1.75
# x = 0.85
#
# m = PAR(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.0101, kposP=0, kAP=0.5, kPA=0.5,
#         ePneg=1, eAneg=1, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)
#
# kon0 = 10 ** kon0
# m.konP = kon0 * (1 - x)
# m.konA = kon0 * (1 - x)
# m.kposP = x * (m.psi * kon0 + m.koffP) / 1
# m.kposA = x * (m.psi * kon0 + m.koffA) / 1
#
# m.initiate()
# m.run(save_direc='_test', save_gap=10)
# animatePAR('_test')
#
# # plt.plot(m.A)
# # plt.plot(m.P)
# # plt.show()
