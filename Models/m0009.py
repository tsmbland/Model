import numpy as np
import pickle

"""
Model with pPAR membrane binding receptor

"""


class Params:
    def __init__(self, pA, Da, konA, koffA, kAP, ePneg, pP, Dp, konP, koffP, kPA, ePA, pS, Ds, konS, koffS, kSA, eSA, L,
                 xsteps, psi, Tmax, deltat):
        ######### A ##########
        self.pA = pA  # um-3
        self.Da = Da  # um2 s-1
        self.konA = konA  # um s-1
        self.koffA = koffA  # s-1
        self.kAP = kAP
        self.ePneg = ePneg

        ######### P ##########
        self.pP = pP  # um-3
        self.Dp = Dp  # um2 s-1
        self.konP = konP  # um s-1
        self.koffP = koffP  # s-1
        self.kPA = kPA
        self.ePA = ePA

        ######### S ##########
        self.pS = pS  # um-3
        self.Ds = Ds  # um2 s-1
        self.konS = konS  # um s-1
        self.koffS = koffS  # s-1
        self.kSA = kSA
        self.eSA = eSA

        ######## Misc ########
        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s


class Model:
    def __init__(self, params):
        """
        a: Cortical aPAR
        p: Cortical pPAR
        s: Cortical scaffold
        """

        self.params = params
        self.res = self.Res(params)

        self.a = np.zeros([self.params.xsteps])
        self.p = np.zeros([self.params.xsteps])
        self.s = np.zeros([self.params.xsteps])

    def reactions(self, lb, ub):
        """
        r0: binding of a to cortex
        r1: unbinding of a from cortex
        r2: antagonism from p to a
        r3: diffusion of a
        
        r4: binding of p to scaffold
        r5: unbinding of p from scaffold
        r6: antagonism from a to p
        r7: diffusion of p
        
        r8: binding of scaffold to cortex
        r9: unbinding of scaffold from cortex
        r10: antagonism from a to s
        r11: diffusion of s
        """

        acyt = (self.params.pA - self.params.psi * np.mean(self.a))
        pcyt = (self.params.pP - self.params.psi * np.mean(self.p))
        scyt = (self.params.pS - self.params.psi * np.mean(self.s))

        r = np.zeros([12])
        r[0] = self.params.konA * acyt
        r[1] = self.params.koffA * self.a[lb:ub]
        r[2] = self.params.kAP * self.a[lb:ub] * (self.p[lb:ub] ** self.params.eAP)
        r[3] = self.params.Da * self.diffusion(self.a[lb:ub])
        r[4] = self.params.konP * pcyt
        r[5] = self.params.koffP * self.p[lb:ub]
        r[6] = self.params.kPA * self.p[lb:ub] * (self.a[lb:ub] ** self.params.ePA)
        r[7] = self.params.Dp * self.diffusion(self.p[lb:ub])
        r[8] = self.params.konS * scyt
        r[9] = self.params.koffS * self.s[lb:ub]
        r[10] = self.params.kSA * self.s[lb:ub] * (self.a[lb:ub] ** self.params.eSA)
        r[11] = self.params.Ds * self.diffusion(self.s[lb:ub])
        return r

    def diffusion(self, concs):
        diff = (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_a(self, r):
        self.a += (r[0] - r[1] - r[2] + r[3]) * self.params.deltat

    def update_p(self, r):
        self.a += (r[4] - r[5] - r[6] + r[7]) * self.params.deltat

    def update_s(self, r):
        self.a += (r[8] - r[9] - r[10] + r[11]) * self.params.deltat

    def get_all(self):
        return [self.a, self.p, self.s]

    def run(self):
        for t in range(int(self.params.Tmax / self.params.deltat)):
            r1 = self.reactions(0, self.params.xsteps // 2)
            self.update_a(r1)
            self.update_p(r1)
            self.update_s(r1)

            r2 = self.reactions(self.params.xsteps // 2, self.params.xsteps)
            self.update_a(r2)
            self.update_p(r2)
            self.update_s(r2)
        self.res.update(-1, self.get_all())

        for t in range(int(self.params.Tmax / self.params.deltat)):
            r = self.reactions(0, self.params.xsteps)
            self.update_a(r)
            self.update_p(r)
            self.update_s(r)
            self.res.update(t, self.get_all())

        return self.res

    class Res:
        def __init__(self, params):
            self.params = params
            self.a = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.p = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.s = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

        def update(self, t, c):
            self.a[t + 1, :] = c[0]
            self.p[t + 1, :] = c[1]
            self.s[t + 1, :] = c[2]
            self.aco[t + 1, :] = c[0]
            self.pco[t + 1, :] = c[1]

        def compress(self):
            self.a = np.asarray([self.a[-1, :], ])
            self.p = np.asarray([self.a[-1, :], ])
            self.s = np.asarray([self.a[-1, :], ])
            self.aco = np.asarray([self.aco[-1, :], ])
            self.pco = np.asarray([self.pco[-1, :], ])


# p0 = Params(pA=1, Da=1, konA=, koffA=, kAP=, ePneg=1, pP=1, Dp=1, konP=, koffP=, kPA=, ePA=1, pS=1, Ds=1, konS=, koffS=,
#             kSA=, eSA=2, L=67.3, xsteps=500, psi=0.174, Tmax=10, deltat=0.1)

