import numpy as np

"""
Model with antagonism on positive feedback

"""


class Params:
    def __init__(self, Da, Dp, konA1, koffA, konP1, koffP, kAneg, kPneg, ePneg, eAneg, konA2, konP2, kApos, kPpos,
                 eApos, ePpos, pA, pP, L, xsteps, psi, Tmax, deltat, Aeqmin, Aeqmax, Peqmin, Peqmax):
        # Diffusion
        self.Da = Da  # um2 s-1
        self.Dp = Dp  # um2 s-1

        # Membrane exchange
        self.konA1 = konA1  # um s-1
        self.koffA = koffA  # s-1
        self.konP1 = konP1  # um s-1
        self.koffP = koffP  # s-1

        # Positive feedback
        self.konA2 = konA2  # um2 s-1
        self.konP2 = konP2  # um2 s-1
        self.kApos = kApos  # um2
        self.kPpos = kPpos  # um2
        self.eApos = eApos
        self.ePpos = ePpos

        # Antagonism
        self.kPneg = kPneg  # um2
        self.kAneg = kAneg  # um2
        self.ePneg = ePneg
        self.eAneg = eAneg

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
        int_on = (p.konA1 * (p.pA - p.psi * np.mean(self.aco)))
        pf_on = (p.konA2 * (p.pA - p.psi * np.mean(self.aco)) * (self.aco ** p.eApos) / (
            (p.kApos ** p.eApos) + (self.aco ** p.eApos)))
        ant = (p.kPneg ** p.ePneg) / ((p.kPneg ** p.ePneg) + (self.pco ** p.ePneg))
        self.aco += ((diff + int_on - off + pf_on * ant) * p.deltat)

    def update_pco(self, p):
        diff = self.diffusion(self.pco, p.Dp)
        off = (p.koffP * self.pco)
        int_on = (p.konP1 * (p.pP - p.psi * np.mean(self.pco)))
        pf_on = (p.konP2 * (p.pP - p.psi * np.mean(self.pco)) * (self.pco ** p.ePpos) / (
            (p.kPpos ** p.ePpos) + (self.pco ** p.ePpos)))
        ant = (p.kAneg ** p.eAneg) / ((p.kAneg ** p.eAneg) + (self.aco ** p.eAneg))
        self.pco += ((diff + int_on - off + pf_on * ant) * p.deltat)

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
        for t in range(int(self.params.Tmax / self.params.deltat)):
            self.update_aco(self.params)
            self.update_pco(self.params)
            self.res.update(t, self.get_all())

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
