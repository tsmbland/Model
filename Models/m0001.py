import numpy as np

"""
Model with Hill functions for antagonism

"""


class Params:
    def __init__(self, Da, Dp, konA, koffA, konP, koffP, kAP, kPA, kAneg, kPneg, ePneg, eAneg, pA, pP, L, xsteps, psi,
                 Tmax, deltat, Aeqmin, Aeqmax, Peqmin, Peqmax):
        # Diffusion
        self.Da = Da  # um2 s-1
        self.Dp = Dp  # um2 s-1

        # Membrane exchange
        self.konA = konA  # um s-1
        self.koffA = koffA  # s-1
        self.konP = konP  # um s-1
        self.koffP = koffP  # s-1

        # Antagonism
        self.kAP = kAP  # s-1
        self.kPA = kPA  # s-1
        self.kAneg = kAneg  # um-2
        self.kPneg = kPneg  # um-2
        self.ePneg = ePneg  #
        self.eAneg = eAneg  #

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
        self.params = p
        self.res = self.Res(p)

    def diffusion(self, concs, coeff):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_aco(self, p):
        diff = self.diffusion(self.aco, p.Da)
        off = (p.koffA * self.aco)
        on = (p.konA * (p.pA - p.psi * np.mean(self.aco)))
        ant = p.kAP * self.aco * (self.pco ** p.ePneg) / (
            (p.kPneg ** p.ePneg) + (self.pco ** p.ePneg))
        self.aco += ((diff + on - off - ant) * p.deltat)

    def update_pco(self, p):
        diff = self.diffusion(self.pco, p.Dp)
        off = (p.koffP * self.pco)
        on = (p.konP * (p.pP - p.psi * np.mean(self.pco)))
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
        self.equilibrate_aco(self.params.p)
        self.equilibrate_pco(self.params.p)
        self.res.update(-1, self.get_all())

        # Run model
        for t in range(int(self.params.p.Tmax / self.params.p.deltat)):
            self.update_aco(self.params.p)
            self.update_pco(self.params.p)
            self.res.update(t, self.get_all())

        return self.res

    class Res:
        def __init__(self, p):
            self.params = p
            self.scores = {}
            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

        def update(self, t, c):
            self.aco[t + 1] = c[0]
            self.pco[t + 1, :] = c[1]

        def compress(self):
            self.aco = np.asarray([self.aco[-1, :], ])
            self.pco = np.asarray([self.pco[-1, :], ])
