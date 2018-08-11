"""
Original model, but coded differently to explicitly model cytoplasmic species


"""

import numpy as np

"""
Originial model with mass action antagonism and no positive feedback

"""


class Params:
    def __init__(self, Da, Dp, konA, koffA, konP, koffP, kAP, kPA, ePneg, eAneg, L, xsteps, psi, Tmax, deltat,
                 starts):
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

        # Misc
        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s

        # Starts
        self.am_0 = starts[0]
        self.ac_0 = starts[1]
        self.pm_0 = starts[2]
        self.pc_0 = starts[3]


class Model:
    def __init__(self, p):
        """

        :param p:
        """

        self.params = p
        self.am = np.zeros([p.xsteps]) * self.params.am_0
        self.pm = np.zeros([p.xsteps]) * self.params.am_0
        self.ac = self.params.ac_0
        self.pc = self.params.pc_0
        self.res = self.Res(p)

    def diffusion(self, concs, coeff):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def reactions(self, lb, ub):
        """
        r0: a on
        r1: a off
        r2: p to a antagnism
        r3: a diffusion

        r4: p on
        r5: p off
        r6: a to p antagonism
        r7: p diffusion

        """

        r = [None] * 8
        r[0] = self.params.konA * self.ac
        r[1] = self.params.koffA * self.am[lb:ub]
        r[2] = self.params.kAP * (self.pm[lb:ub] ** self.params.ePneg) * self.am[lb:ub]
        r[3] = self.diffusion(self.am[lb:ub], self.params.Da)

        r[4] = self.params.konP * self.pc
        r[5] = self.params.koffP * self.pm[lb:ub]
        r[6] = self.params.kPA * (self.am[lb:ub] ** self.params.eAneg) * self.pm[lb:ub]
        r[7] = self.diffusion(self.pm[lb:ub], self.params.Dp)

        return r

    def update_am(self, r, lb, ub):
        self.am[lb:ub] += (r[0] - r[1] - r[2] + r[3]) * self.params.deltat

    def update_pm(self, r, lb, ub):
        self.pm[lb:ub] += (r[4] - r[5] - r[6] + r[7]) * self.params.deltat

    def update_ac(self, r, lb, ub):
        x = (ub - lb) / self.params.xsteps
        self.ac += (- (x * self.params.psi) * r[0] + (x * self.params.psi) * np.mean(r[1]) + (
            x * self.params.psi) * np.mean(
            r[2])) * self.params.deltat

    def update_pc(self, r, lb, ub):
        x = (ub - lb) / self.params.xsteps
        self.pc += (- (x * self.params.psi) * r[4] + (x * self.params.psi) * np.mean(r[5]) + (
            x * self.params.psi) * np.mean(
            r[6])) * self.params.deltat

    def get_all(self):
        return [self.am, self.ac, self.pm, self.pc]

    def run(self):
        self.__init__(self.params)  # <- temporary fix
        for t in range(int(self.params.Tmax / self.params.deltat)):
            r1 = self.reactions(0, self.params.xsteps // 2)
            self.update_am(r1, 0, self.params.xsteps // 2)
            self.update_ac(r1, 0, self.params.xsteps // 2)

            r2 = self.reactions(self.params.xsteps // 2, self.params.xsteps)
            self.update_pm(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_pc(r2, self.params.xsteps // 2, self.params.xsteps)
            self.res.update(t, self.get_all())
        self.res.update(-1, self.get_all())

        for t in range(int(self.params.Tmax / self.params.deltat)):
            r = self.reactions(0, self.params.xsteps)
            self.update_am(r, 0, self.params.xsteps)
            self.update_ac(r, 0, self.params.xsteps)
            self.update_pm(r, 0, self.params.xsteps)
            self.update_pc(r, 0, self.params.xsteps)
            self.res.update(t, self.get_all())
        return self.res

    class Res:
        def __init__(self, p):
            self.params = p
            self.scores = {}

            self.am = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.ac = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pm = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pc = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

            self.atot = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.ptot = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

        def update(self, t, c):
            self.am[t + 1, :] = c[0]
            self.ac[t + 1] = c[1]
            self.pm[t + 1, :] = c[2]
            self.pc[t + 1] = c[3]

            self.aco[t + 1, :] = c[0]
            self.pco[t + 1, :] = c[2]

            self.atot[t + 1] = c[1] + self.params.psi * np.mean(c[0])
            self.ptot[t + 1] = c[3] + self.params.psi * np.mean(c[2])

        def compress(self):
            self.am = np.asarray([self.am[-1, :], ])
            self.ac = self.ac[-1]
            self.pm = np.asarray([self.pm[-1, :], ])
            self.pc = self.pc[-1]

            self.aco = np.asarray([self.aco[-1, :], ])
            self.pco = np.asarray([self.pco[-1, :], ])

            self.atot = self.atot[-1]
            self.ptot = self.ptot[-1]


###################################################################################

p0 = Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, ePneg=1, eAneg=2,
            L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.1, starts=[0, 1.56, 0, 1])

# p1 = Params(Da=1, Dp=1, konA=1, koffA=0.1, konP=1, koffP=0.1, kAP=1, kPA=1,
#             ePneg=2, eAneg=2, pA=1, pP=1, L=50, xsteps=500, psi=0.3, Tmax=1000, deltat=0.01, Aeqmin=0, Aeqmax=0.5,
#             Peqmin=0.5, Peqmax=1)
#
# p2 = Params(Da=0.1, Dp=0.1, konA=0.006, koffA=0.005, konP=0.006, koffP=0.005, kAP=1, kPA=1,
#             ePneg=2, eAneg=2, pA=0, pP=1, L=50, xsteps=500, psi=0.3, Tmax=1000, deltat=0.01, Aeqmin=0, Aeqmax=0.5,
#             Peqmin=0.5, Peqmax=1)
