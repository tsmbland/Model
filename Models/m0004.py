import numpy as np

"""
Model with pPAR state switch

"""


class Params:
    def __init__(self, Da, Dp, konA, koffA, konP, koffP, kAP, kPA, kPneg1, kPneg2, kAneg1, kAneg2, ePneg, eAneg,
                 kApos, kPpos, eApos, ePpos, pA, pP, L, xsteps, psi, Tmax, deltat, Aeqmin, Aeqmax, Peqmin, Peqmax):
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
        self.kPneg1 = kPneg1  # um2
        self.kPneg2 = kPneg2  # um2
        self.kAneg1 = kAneg1  # um2
        self.kAneg2 = kAneg2  # um2
        self.ePneg = ePneg
        self.eAneg = eAneg

        # Positive feedback
        self.kApos = kApos  # um2
        self.kPpos = kPpos  # um2
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
        self.params = p
        self.aco = np.zeros([p.xsteps])
        self.pco = np.zeros([p.xsteps])
        self.res = self.Res(p)

    def diffusion(self, concs, coeff):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_aco(self, p):
        diff = self.diffusion(self.aco, p.Da)
        off = (p.koffA * self.aco)
        on = (p.konA * (p.pA - p.psi * np.mean(self.aco)))
        k = p.kPneg1 + p.kPneg2 * ((self.aco ** p.eApos) / ((p.kApos ** p.eApos) + (self.aco ** p.eApos)))
        ant = p.kAP * self.aco * ((self.pco ** p.ePneg) / ((k ** p.ePneg) + (self.pco ** p.ePneg)))
        self.aco += ((diff + on - off - ant) * p.deltat)

    def update_pco(self, p):
        diff = self.diffusion(self.pco, p.Dp)
        off = (p.koffP * self.pco)
        on = (p.konP * (p.pP - p.psi * np.mean(self.pco)))
        k = p.kAneg1 + p.kAneg2 * (self.pco ** p.ePpos) / ((p.kPpos ** p.ePpos) + (self.pco ** p.ePpos))
        ant = p.kPA * self.pco * (self.aco ** p.eAneg) / ((k ** p.eAneg) + (self.aco ** p.eAneg))
        self.pco += ((diff + on - off - ant) * p.deltat)

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


###################################################################################


p1 = Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0, kPA=0.8,
            kPneg1=1.5, kPneg2=0, kAneg1=0.5, kAneg2=3.5, ePneg=10, eAneg=10, kApos=1, kPpos=0.5,
            eApos=1, ePpos=10, pA=1.56, pP=1, L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.1, Aeqmin=0,
            Aeqmax=0.5, Peqmin=0.5, Peqmax=1)

p0 = Params(Da=1, Dp=1, konA=1, koffA=0.3, konP=1, koffP=0.3, kAP=3, kPA=3,
            kPneg1=1, kPneg2=1, kAneg1=1, kAneg2=1, ePneg=20, eAneg=20, kApos=1, kPpos=1, eApos=20,
            ePpos=20, pA=1, pP=1, L=50, xsteps=500, psi=0.3, Tmax=1000, deltat=0.01, Aeqmin=0, Aeqmax=0.5, Peqmin=0.5,
            Peqmax=1)
