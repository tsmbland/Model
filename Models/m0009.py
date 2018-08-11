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
        self.pm = self.params.pm_0 * np.zeros([self.params.xsteps])
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
        r19: Antagonism from am to spm0
        r20: Antagonism from am to spm1

        r21: diffusion of am
        r22: diffusion of sm
        r23: diffusion of pm
        r24: diffusion of spm0
        r24: diffusion of spm1


        """

        r = [None] * 12
        r[0] = self.params. * self.ac
        r[1] = self.params. * self.am[lb:ub]

        r[2] = self.params. * self.sc
        r[3] = self.params. * self.sm[lb:ub]

        r[4] = self.params. * self.pc
        r[5] = self.params. * self.pm[lb:ub]

        r[6] = self.params. * self.sm * self.pc
        r[7] = self.params. * self.spm0
        r[8] = self.params. * self.spc
        r[9] = self.params. * self.spm0
        r[10] = self.params. * self.spm0
        r[11] = self.params. * self.spm1
        r[12] = self.params. * self.spm1
        r[13] = self.params. * self.pm * self.sm
        r[14] = self.params. * self.pc * self.sc
        r[15] = self.params. * self.spc
        r[16] = self.params. * (...) * self.am
        r[17] = self.params. * self.am * self.sm
        r[18] = self.params. * self.am * self.pm
        r[19] = self.params. * self.am * self.spm0
        r[20] = self.params. * self.am * self.spm1

        r[21] = self.params. * self.diffusion(self.am)
        r[22] = self.params. * self.diffusion(self.sm)
        r[23] = self.params. * self.diffusion(self.pm)
        r[24] = self.params. * self.diffusion(self.spm0)
        r[25] = self.params. * self.diffusion(self.spm1)

        return r

    def diffusion(self, concs):
        diff = (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_am(self, r, lb, ub):
        self.am[lb:ub] += (...) * self.params.deltat

    def update_ac(self, r, lb, ub):
        self.ac += (...) * self.params.deltat

    def update_pm(self, r, lb, ub):
        self.pm[lb:ub] += (...) * self.params.deltat

    def update_pc(self, r, lb, ub):
        self.pc += (...) * self.params.deltat

    def update_sm(self, r, lb, ub):
        self.sm[lb:ub] += (...) * self.params.deltat

    def update_sc(self, r, lb, ub):
        self.sc += (...) * self.params.deltat

    def update_spm0(self, r, lb, ub):
        self.spm0[lb:ub] += (...) * self.params.deltat

    def update_spm1(self, r, lb, ub):
        self.spm1[lb:ub] += (...) * self.params.deltat

    def update_spc(self, r, lb, ub):
        self.spc += (...) * self.params.deltat

    def get_all(self):
        return [self.am, self.ac, self.pm, self.pc, self.sm, self.sc, self.spm0, self.spm1, self.spc]

    def run(self):
        for t in range(int(self.params.Tmax / self.params.deltat)):
            r1 = self.reactions(0, self.params.xsteps // 2)
            self.update_am(r1, 0, self.params.xsteps // 2)
            self.update_ac(r1, 0, self.params.xsteps // 2)

            r2 = self.reactions(self.params.xsteps // 2, self.params.xsteps)
            self.update_pm(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_pc(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_sm(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_sc(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_spm0(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_spm1(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_spc(r2, self.params.xsteps // 2, self.params.xsteps)
        self.res.update(-1, self.get_all())

        for t in range(int(self.params.Tmax / self.params.deltat)):
            r = self.reactions(0, self.params.xsteps)
            self.update_am(r1, 0, self.params.xsteps)
            self.update_ac(r1, 0, self.params.xsteps)
            self.update_pm(r2, 0, self.params.xsteps)
            self.update_pc(r2, 0, self.params.xsteps)
            self.update_sm(r2, 0, self.params.xsteps)
            self.update_sc(r2, 0, self.params.xsteps)
            self.update_spm0(r2, 0, self.params.xsteps)
            self.update_spm1(r2, 0, self.params.xsteps)
            self.update_spc(r2, 0, self.params.xsteps)
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
