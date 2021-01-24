import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc + '/../..')

import numpy as np
from parmodel import pdeRK, diffusion
from models.ode.Goehring2011 import Goehring2011 as ODE
from scipy.integrate import odeint


class Goehring2011:
    def __init__(self, Da=0.28, Dp=0.15, konA=0.00858, koffA=0.0054, konP=0.0474, koffP=0.0073, kPA=2, kAP=0.19,
                 alpha=1, beta=2, xsteps=100, psi=0.174, Tmax=1000, deltat=0.01, L=134.6, pA=1.56, pP=1):
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

        # Antagonism
        self.kPA = kPA  # um4 s-1
        self.kAP = kAP  # um2 s-1
        self.alpha = alpha
        self.beta = beta

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
        dA = ((self.konA * ac) - (self.koffA * A) - (self.kAP * (P ** self.alpha) * A) + (
                self.Da * diffusion(A, self.deltax)))
        dP = ((self.konP * pc) - (self.koffP * P) - (self.kPA * (A ** self.beta) * P) + (
                self.Dp * diffusion(P, self.deltax)))
        return [dA, dP]

    def initiate(self):
        """
        Initiating the system polarised

        """

        # Solve ode, no antagonism
        o = ODE(konA=self.konA, koffA=self.koffA, konP=self.konP, koffP=self.koffP, alpha=self.alpha, beta=self.beta,
                psi=self.psi, pA=self.pA, pP=self.pP, kAP=0, kPA=0)
        soln = odeint(o.dxdt, (0, 0), t=np.linspace(0, 10000, 100000))[-1]

        self.A = soln[0]
        self.P = soln[1]

        # Polarise
        self.A *= 2 * np.r_[np.ones([self.xsteps // 2]), np.zeros([self.xsteps // 2])]
        self.P *= 2 * np.r_[np.zeros([self.xsteps // 2]), np.ones([self.xsteps // 2])]

    def initiate2(self):
        """
        Initiating the system (near) uniform, A dominant

        """

        # Solve ode
        o = ODE(konA=self.konA, koffA=self.koffA, konP=self.konP, koffP=self.koffP, alpha=self.alpha, beta=self.beta,
                psi=self.psi, pA=self.pA, pP=self.pP, kAP=0, kPA=0)
        soln = odeint(o.dxdt, (o.pA / o.psi, 0), t=np.linspace(0, 10000, 100000))[-1]

        # Set concentrations
        self.A[:] = soln[0]
        self.P[:] = soln[1]

        # Polarise
        self.A *= np.linspace(1.01, 0.99, self.xsteps)
        self.P *= np.linspace(0.99, 1.01, self.xsteps)

    def run(self, save_direc=None, save_gap=None, kill_uni=False, kill_stab=False):
        """

        :param save_direc: if given, will save A and P distributions over time according to save_gap
        :param save_gap: gap in model time between save points
        :param kill_uni: terminate once polarity is lost. Generally can assume models never regain polarity once lost
        :param kill_stab: terminate when patterns are stable. I'd advise against for phase-space diagrams, can get
            fuzzy boundaries
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
