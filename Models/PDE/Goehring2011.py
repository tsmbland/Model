import numpy as np
from Funcs import pdeRK, diffusion


class Goehring2011:
    def __init__(self, Da, Dp, konA, koffA, konP, koffP, kPA, kAP, alpha, beta, xsteps, psi, Tmax, deltat, L, pA, pP):
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
        dA = ((self.konA * ac) - (self.koffA * A) - (self.kAP * (P ** self.beta) * A) + (
            self.Da * diffusion(A, self.deltax)))
        dP = ((self.konP * pc) - (self.koffP * P) - (self.kPA * (A ** self.alpha) * P) + (
            self.Dp * diffusion(P, self.deltax)))
        return [dA, dP]

    def initiate(self):
        # Remove antagonism
        kAP, kPA = self.kAP, self.kPA
        self.kAP = 0
        self.kPA = 0

        # Initial equilibration (no antagonism)
        soln, time, solns, times = pdeRK(dxdt=self.dxdt, X0=[self.A, self.P], Tmax=self.Tmax, deltat=self.deltat,
                                         t_eval=np.arange(0, self.Tmax + 0.0001, self.Tmax), stabilitycheck=True)
        self.A = soln[0]
        self.P = soln[1]

        # Polarise
        self.A *= 2 * np.r_[np.ones([self.xsteps // 2]), np.zeros([self.xsteps // 2])]
        self.P *= 2 * np.r_[np.zeros([self.xsteps // 2]), np.ones([self.xsteps // 2])]

        # Add back antagonism
        self.kAP = kAP
        self.kPA = kPA

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
