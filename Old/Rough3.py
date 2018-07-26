import numpy as np


class Model0:
    """Originial model with mass action antagonism and no positive feedback"""

    class Params:
        def __init__(self, Da, Dp, konA, koffA, konP, koffP, kAP, kPA, ePneg, eAneg, pA, pP, pgen):
            # Diffusion
            self.Da = Da  # um2 s-1
            self.Dp = Dp  # um2 s-1

            # Membrane exchange
            self.konA = konA  # um s-1
            self.koffA = koffA  # s-1
            self.konP = konP  # um s-1
            self.koffP = koffP  # s-1

            # Antagonism
            self.kAP = kAP  # um2 s-1
            self.kPA = kPA  # um4 s-1
            self.ePneg = ePneg
            self.eAneg = eAneg

            # Pools
            self.pA = pA  # um-3
            self.pP = pP  # um-3

            # Misc
            self.L = pgen.L  # um
            self.xsteps = pgen.xsteps
            self.psi = pgen.psi  # um-1
            self.Tmax = pgen.Tmax  # s
            self.deltat = pgen.deltat  # s

            # Equilibration
            self.Aeqmin = pgen.Aeqmin
            self.Aeqmax = pgen.Aeqmax
            self.Peqmin = pgen.Peqmin
            self.Peqmax = pgen.Peqmax

    class Molecule:
        def __init__(self, cort0, cyt0):
            self.cort = cort0
            self.cyt = cyt0

    def __init__(self, p):
        self.a = Molecule()
        self.p1 = Molecule()
        self.p2 = Molecule()
        self.p3 = Molecule()

    def update_a(self, p):
        diff = diffusion(self.a.cort, p.Da, p)
        off = (p.koffA * self.a.cort)
        on = (p.konA * (p.pA - p.psi * np.mean(self.a.cort)))
        ant = p.kAP * (self.pco ** p.ePneg) * self.a.cort
        self.aco += ((diff + on - off - ant) * p.deltat)

    def update_p1(self, p):
        diff = diffusion(self.pco, p.Dp, p)
        off = (p.koffP * self.pco)
        on = (p.konP * (p.pP - p.psi * np.mean(self.pco)))
        ant = p.kPA * (self.aco ** p.eAneg) * self.pco
        self.pco += ((diff + on - off - ant) * p.deltat)

    def update_p2(self, p):
        diff = diffusion(self.pco, p.Dp, p)
        off = (p.koffP * self.pco)
        on = (p.konP * (p.pP - p.psi * np.mean(self.pco)))
        ant = p.kPA * (self.aco ** p.eAneg) * self.pco
        self.pco += ((diff + on - off - ant) * p.deltat)

    def update_p3(self, p):
        diff = diffusion(self.pco, p.Dp, p)
        off = (p.koffP * self.pco)
        on = (p.konP * (p.pP - p.psi * np.mean(self.pco)))
        ant = p.kPA * (self.aco ** p.eAneg) * self.pco
        self.pco += ((diff + on - off - ant) * p.deltat)


class k:
    pona = 0.01
    ponb = 0.1
    poffa = 0.01
    poffb = 1
    kphos = 1
    kdephos = 1
    aon = 0.01
    aoff = 0.01
    kAP = 0
    Da = 1
    Dp = 1


def diffusion(concs, coeff, p):
    if hasattr(concs, '__len__'):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (p.L / p.xsteps)
    else:
        diff = 0
    return diff


class MiscParams:
    """
    General parameters shared between all models

    """

    def __init__(self, L, xsteps, psi, Tmax, deltat):
        # Misc
        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s


p = MiscParams(L=50, xsteps=500, psi=0.3, Tmax=100, deltat=0.01)
