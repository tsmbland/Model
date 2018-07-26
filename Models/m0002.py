import numpy as np
import pickle

"""
Model with Hill functions for antagonism and positive feedback

"""


class Params:
    def __init__(self, Da, Dp, konA1, koffA, konP1, koffP, kAP, kPA, kAneg, kPneg, ePneg, eAneg, konA2, konP2,
                 kApos, kPpos, eApos, ePpos, pA, pP, L, xsteps, psi, Tmax, deltat, Aeqmin, Aeqmax, Peqmin, Peqmax):
        # Diffusion
        self.Da = Da  # um2 s-1
        self.Dp = Dp  # um2 s-1

        # Membrane exchange
        self.konA1 = konA1  # um s-1
        self.koffA = koffA  # s-1
        self.konP1 = konP1  # um s-1
        self.koffP = koffP  # s-1

        # Antagonism
        self.kAP = kAP  # s-1
        self.kPA = kPA  # s-1
        self.kAneg = kAneg  # um-2
        self.kPneg = kPneg  # um-2
        self.ePneg = ePneg  #
        self.eAneg = eAneg  #

        # Positive feedback
        self.konA2 = konA2  # um2 s-1
        self.konP2 = konP2  # um2 s-1
        self.kApos = kApos  # um-2
        self.kPpos = kPpos  # um-2
        self.eApos = eApos
        self.ePpos = ePpos

        # Pools
        self.pA = pA  # um-3
        self.pP = pP  # um-3

        # Misc
        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s

        # Equilibration
        self.Aeqmin = Aeqmin
        self.Aeqmax = Aeqmax
        self.Peqmin = Peqmin
        self.Peqmax = Peqmax


class Model:
    def __init__(self, p):
        self.aco = np.zeros([p.xsteps])
        self.pco = np.zeros([p.xsteps])
        self.res = self.Res(p)
        self.params = p

    def diffusion(self, concs, coeff):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_aco(self, p):
        diff = self.diffusion(self.aco, p.Da)
        off = (p.koffA * self.aco)
        on = (p.konA1 * (p.pA - p.psi * np.mean(self.aco))) + (
            p.konA2 * (p.pA - p.psi * np.mean(self.aco)) * (self.aco ** p.eApos) / (
                (p.kApos ** p.eApos) + (self.aco ** p.eApos)))
        ant = p.kAP * self.aco * (self.pco ** p.ePneg) / (
            (p.kPneg ** p.ePneg) + (self.pco ** p.ePneg))
        self.aco += ((diff + on - off - ant) * p.deltat)

    def update_pco(self, p):
        diff = self.diffusion(self.pco, p.Dp)
        off = (p.koffP * self.pco)
        on = (p.konP1 * (p.pP - p.psi * np.mean(self.pco))) + (
            p.konP2 * (p.pP - p.psi * np.mean(self.pco)) * (self.pco ** p.ePpos) / (
                (p.kPpos ** p.ePpos) + (self.pco ** p.ePpos)))
        ant = p.kPA * self.pco * (self.aco ** p.eAneg) / (
            (p.kAneg ** p.eAneg) + (self.aco ** p.eAneg))
        self.pco += (diff + on - off - ant) * p.deltat

    def equilibrate_aco(self, p):
        for t in range(5000):
            self.update_aco(p)
            self.aco[:int(p.xsteps * p.Aeqmin)] = 0
            self.aco[int(p.xsteps * p.Aeqmax):] = 0

    def equilibrate_pco(self, p):
        for t in range(5000):
            self.update_pco(p)
            self.pco[:int(p.xsteps * p.Peqmin)] = 0
            self.pco[int(p.xsteps * p.Peqmax):] = 0

    def get_all(self):
        return [self.aco, self.pco]

    def run(self):

        # Equilibrate
        self.equilibrate_aco(self.params)
        self.equilibrate_pco(self.params)
        self.res.update(-1, self.get_all())

        # Run model
        for t in range(int(self.params.Tmax / self.params.deltat)):
            self.update_aco(p)
            self.update_pco(p)
            res.update(t, self.get_all())

        return self.res

    class Res:
        def __init__(self, p):
            self.params = p
            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

        def update(self, t, c):
            self.aco[t + 1] = c[0]
            self.pco[t + 1, :] = c[1]

        def compress(self):
            self.aco = np.asarray([self.aco[-1, :], ])
            self.pco = np.asarray([self.pco[-1, :], ])


###################################################################################


p1 = Params(Da=0.28, Dp=0.15, konA1=0.0085, koffA=0.0054, konP1=0.0474, koffP=0.0073, kAP=0, kPA=0.15,
            kAneg=0.3, kPneg=0.3, ePneg=1, eAneg=10, konA2=0, konP2=1, kApos=0, kPpos=2, eApos=0,
            ePpos=10, pA=1.56, pP=1, L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.1, Aeqmin=0, Aeqmax=0.5,
            Peqmin=0.5, Peqmax=1)

p2 = Params(Da=0.28, Dp=0.15, konA1=0.0085, koffA=0.0054, konP1=0.02539162, koffP=0.0073, kAP=0,
            kPA=0.15, kAneg=0.3, kPneg=0.3, ePneg=1, eAneg=10, konA2=0, konP2=0.00257667, kApos=0,
            kPpos=1, eApos=0, ePpos=10, pA=1.56, pP=1, L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.1, Aeqmin=0,
            Aeqmax=0.5, Peqmin=0.5, Peqmax=1)
