import numpy as np

"""
Model with pPAR membrane binding receptor
INCOMPLETE

"""


class Params:
    def __init__(self, pA, Da, konA, koffA, kAP, eAP, pP, Dp, konP, koffP, kPA, ePA, pS, Ds, konS, koffS, kSA, eSA, L,
                 xsteps, psi, Tmax, deltat):
        ######### A ##########
        self.pA = pA  # um-3
        self.Da = Da  # um2 s-1
        self.konA = konA  # um s-1
        self.koffA = koffA  # s-1
        self.kAP = kAP
        self.eAP = eAP

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
        am: Membrane aPAR
        ac: Cytoplasmic aPAR
        pc: Cytoplasmic pPAR
        sm: Membrane scaffold (unbound to p)
        sc: Cytoplasmic scaffold
        spm0: Membrane scaffold with non-membrane-bound p
        spm1: Membrane scaffold with membrane-bound p
        spc: Cytoplasmic scaffold-p

        """

        self.params = params
        self.res = self.Res(params)

        self.am = self.params.am_0 * np.zeros([self.params.xsteps])
        self.ac = self.params.ac_0
        self.am = self.params.am_0 * np.zeros([self.params.xsteps])
        self.pc = self.params.pc_0
        self.sm = self.params.sm_0 * np.zeros([self.params.xsteps])
        self.sc = self.params.sc_0
        self.spm0 = self.params.spm0_0 * np.zeros([self.params.xsteps])
        self.spm1 = self.params.spm1_0 * np.zeros([self.params.xsteps])
        self.spc = self.params.spc_0

    def reactions(self, lb, ub):
        """

        r0: binding of ac to cortex
        r1: unbinding of am from cortex
        r2: binding of sc to cortex
        r3: unbinding of sm from cortex
        r4: binding of pc to cortex
        r5: unbinding of pm from cortex

        r6: Dimerisation of sm with pc
        r7: Separation of spm0
        r8: Single binding of spc to cortex
        r9: unbinding of spm0 from cortex
        r10: binding of p in spm0 to cortex
        r11: unbinding of p in spm1
        r12: separation of spm1
        r13: dimerisation of pm and sm
        r14: dimerisation of pc and sc
        r15: separation of spc

        r16: antagonism from p to am
        r17: antagonism from am to sm
        r18: antagonism from am to pm
        r19: Antagonism from a to smp0

        r20: diffusion of a
        r21: diffusion of s
        r22: diffusion of p
        r23: diffusion of smp0


        """

        acyt = (self.params.pA - self.params.psi * np.mean(self.a))
        pcyt = (self.params.pP - self.params.psi * np.mean(self.p))
        scyt = (self.params.pS - self.params.psi * np.mean(self.s))

        r = [None] * 12
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

    def update_a(self, r, lb, ub):
        self.a[lb:ub] += (r[0] - r[1] - r[2] + r[3]) * self.params.deltat

    def update_p(self, r, lb, ub):
        self.p[lb:ub] += (r[4] - r[5] - r[6] + r[7]) * self.params.deltat

    def update_s(self, r, lb, ub):
        self.s[lb:ub] += (r[8] - r[9] - r[10] + r[11]) * self.params.deltat

    def get_all(self):
        return [self.a, self.p, self.s]

    def run(self):
        for t in range(int(self.params.Tmax / self.params.deltat)):
            r1 = self.reactions(0, self.params.xsteps // 2)
            self.update_a(r1, 0, self.params.xsteps // 2)

            r2 = self.reactions(self.params.xsteps // 2, self.params.xsteps)
            self.update_p(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_s(r2, self.params.xsteps // 2, self.params.xsteps)
        self.res.update(-1, self.get_all())

        for t in range(int(self.params.Tmax / self.params.deltat)):
            r = self.reactions(0, self.params.xsteps)
            self.update_a(r, 0, self.params.xsteps)
            self.update_p(r, 0, self.params.xsteps)
            self.update_s(r, 0, self.params.xsteps)
            self.res.update(t, self.get_all())

        return self.res

    class Res:
        def __init__(self, params):
            self.params = params
            self.scores = {}
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


p0 = Params(pA=1.56, Da=0.28, konA=1, koffA=1, kAP=1, eAP=1, pP=1, Dp=0.15, konP=1, koffP=1, kPA=1, ePA=2, pS=1,
            Ds=0, konS=1, koffS=1, kSA=1, eSA=2, L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.01)
