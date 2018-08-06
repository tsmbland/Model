import numpy as np

"""
pPAR dimerisation model

"""


class Params:
    def __init__(self, Da, kon_a, koff_a, ra, Dp, kon_p, koff_p, kon_p_2, kd_f, kd_b, rp, L, xsteps, psi, Tmax, deltat,
                 starts):
        ########### aPARS ###########

        # Diffusion
        self.Da = Da  # um2 s-1

        # Membrane exchange
        self.kon_a = kon_a
        self.koff_a = koff_a

        # Antagonism
        self.ra = ra

        ########### pPARS ###########

        # Diffusion
        self.Dp = Dp  # um2 s-1

        # Membrane exchange
        self.kon_p = kon_p
        self.koff_p = koff_p
        self.kon_p_2 = kon_p_2

        # Dimerisation
        self.kd_f = kd_f
        self.kd_b = kd_b

        # Antagonism
        self.rp = rp

        ########### Misc ###########

        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s

        ########## Starts ##########

        self.am_0 = starts[0]
        self.ac_0 = starts[1]
        self.pm1_0 = starts[2]
        self.pm2s_0 = starts[3]
        self.pm2d_0 = starts[4]
        self.pc1_0 = starts[5]
        self.pc2_0 = starts[6]


class Model:
    def __init__(self, params):
        """
        am:
        ac:
        pm1:
        pm2s:
        pm2d:
        pc1:
        pc2:
        """

        self.params = params
        self.res = self.Res(params)
        self.am = self.params.am_0 * np.ones([self.params.xsteps])
        self.ac = self.params.ac_0
        self.pm1 = self.params.pm1_0 * np.ones([self.params.xsteps])
        self.pm2s = self.params.pm2s_0 * np.ones([self.params.xsteps])
        self.pm2d = self.params.pm2d_0 * np.ones([self.params.xsteps])
        self.pc1 = self.params.pc1_0
        self.pc2 = self.params.pc2_0

    def reactions(self, lb, ub):
        """
        r1: binding of monomers to cortex
        r2: unbinding of monomers from cortex
        r3: Dimerisation of pm1
        r4: Separation of pm2d
        r5: Dimerisation of pm1 with pc1
        r6: Separation of pm2s
        r7: Antagonism from a to pm1
        r9: Single binding of pc2 to cortex
        r10: unbinding of pm2s from cortex
        r11: binding of second component of dimer to cortex
        r12: unbinding of one component of pmd2
        r13: antagonism from a to pm2s
        r15: antagonism from a to pm2d
        r21: dimerisation of pc1
        r22: separation of pc2

        r17: binding of ac to cortex
        r18: unbinding of am from cortex
        r19: antagonism from p to am

        """

        r = [None] * 28
        r[1] = self.params.kon_p * self.pc1
        r[2] = self.params.koff_p * self.pm1[lb:ub]
        r[3] = self.params.kd_f * self.pm1[lb:ub] ** 2
        r[4] = self.params.kd_b * self.pm2d[lb:ub]
        r[5] = self.params.kd_f * self.pc1 * self.pm1[lb:ub]
        r[6] = self.params.kd_b * self.pm2s[lb:ub]
        r[7] = self.params.rp * self.am[lb:ub] * self.pm1[lb:ub]
        r[9] = self.params.kon_p * self.pc2
        r[10] = self.params.koff_p * self.pm2s[lb:ub]
        r[11] = self.params.kon_p_2 * self.pm2s[lb:ub]
        r[12] = self.params.koff_p * self.pm2d[lb:ub]
        r[13] = self.params.rp * self.am[lb:ub] * self.pm2s[lb:ub]
        r[15] = self.params.rp * self.am[lb:ub] * self.pm2d[lb:ub]
        r[21] = self.params.kd_f * self.pc1 ** 2
        r[22] = self.params.kd_b * self.pc2
        r[17] = self.params.kon_a * self.ac
        r[18] = self.params.koff_a * self.am[lb:ub]
        r[19] = self.params.ra * (self.pm1 + self.pm2s + self.pm2d)[lb:ub] * self.am[lb:ub]
        r[24] = self.diffusion(self.am[lb:ub], self.params.Da)
        r[25] = self.diffusion(self.pm1[lb:ub], self.params.Dp)
        r[26] = self.diffusion(self.pm2s[lb:ub], self.params.Dp)
        r[27] = self.diffusion(self.pm2d[lb:ub], self.params.Dp)
        return r

    def diffusion(self, concs, coeff):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_am(self, r, lb, ub):
        """
        Cortical aPAR
        """
        self.am[lb:ub] += (r[24] + r[17] - r[18] - r[19]) * self.params.deltat

    def update_ac(self, r, lb, ub):
        """
        Cytoplasmic aPAR

        """
        x = (ub - lb) / self.params.xsteps

        self.ac += (- (x * self.params.psi) * r[17] + (x * self.params.psi) * np.mean(r[18]) + (
            x * self.params.psi) * np.mean(
            r[19])) * self.params.deltat

    def update_pm1(self, r, lb, ub):
        """
        Cortical monomer
        """
        self.pm1[lb:ub] += (r[25] + r[1] - r[2] - 2 * r[3] + 2 * r[4] - r[5] + r[6] - r[7]) * self.params.deltat

    def update_pm2s(self, r, lb, ub):
        """
        Cortical dimer (singly bound)
        """
        self.pm2s[lb:ub] += (r[26] + r[5] - r[6] + r[9] - r[10] - r[11] + r[12] - r[13]) * self.params.deltat

    def update_pm2d(self, r, lb, ub):
        """
        Cortical dimer (doubly bound)
        """
        self.pm2d[lb:ub] += (r[27] + r[3] - r[4] + r[11] - r[12] - r[15]) * self.params.deltat

    def update_pc1(self, r, lb, ub):
        """
        Cytoplasmic monomer
        """
        x = (ub - lb) / self.params.xsteps
        self.pc1 += ((- (x * self.params.psi) * r[1]) + ((x * self.params.psi) * np.mean(r[2])) - (
            (x * self.params.psi) * np.mean(
                r[5])) + ((x * self.params.psi) * np.mean(r[6])) + (x * self.params.psi) * np.mean(r[7]) - 2 * r[
                         21] + 2 *
                     r[
                         22]) * self.params.deltat

    def update_pc2(self, r, lb, ub):
        """
        Cytoplasmic dimer
        """
        x = (ub - lb) / self.params.xsteps
        self.pc2 += (- (x * self.params.psi) * r[9] + (x * self.params.psi) * np.mean(r[10]) + (
            x * self.params.psi) * np.mean(
            r[13]) + (x * self.params.psi) * np.mean(r[15]) + r[21] - r[22]) * self.params.deltat

    def get_all(self):
        return [self.am, self.ac, self.pm1, self.pm2s, self.pm2d, self.pc1, self.pc2]

    def run(self):
        self.__init__(self.params)  # <- temporary fix
        for t in range(int(self.params.Tmax / self.params.deltat)):
            r1 = self.reactions(0, self.params.xsteps // 2)
            self.update_am(r1, 0, self.params.xsteps // 2)
            self.update_ac(r1, 0, self.params.xsteps // 2)

            r2 = self.reactions(self.params.xsteps // 2, self.params.xsteps)
            self.update_pm1(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_pm2s(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_pm2d(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_pc1(r2, self.params.xsteps // 2, self.params.xsteps)
            self.update_pc2(r2, self.params.xsteps // 2, self.params.xsteps)
        self.res.update(-1, self.get_all())

        for t in range(int(self.params.Tmax / self.params.deltat)):
            r = self.reactions(0, self.params.xsteps)
            self.update_am(r, 0, self.params.xsteps)
            self.update_ac(r, 0, self.params.xsteps)
            self.update_pm1(r, 0, self.params.xsteps)
            self.update_pm2s(r, 0, self.params.xsteps)
            self.update_pm2d(r, 0, self.params.xsteps)
            self.update_pc1(r, 0, self.params.xsteps)
            self.update_pc2(r, 0, self.params.xsteps)
            self.res.update(t, self.get_all())
        return self.res

    class Res:
        def __init__(self, params):
            self.params = params
            self.scores = {}
            self.am = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.ac = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pm1 = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pm2s = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pm2d = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pc1 = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pc2 = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

            self.pc1 = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pc2 = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

            self.atot = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.ptot = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

        def update(self, t, c):
            self.am[t + 1, :] = c[0]
            self.ac[t + 1] = c[1]
            self.pm1[t + 1, :] = c[2]
            self.pm2s[t + 1, :] = c[3]
            self.pm2d[t + 1, :] = c[4]
            self.pc1[t + 1] = c[5]
            self.pc2[t + 1] = c[6]

            self.aco[t + 1, :] = c[0]
            self.pco[t + 1, :] = c[2] + c[3] + c[4]

            self.atot[t + 1] = c[1] + self.params.psi * np.mean(c[0])
            self.ptot[t + 1] = c[5] + 2 * c[6] + self.params.psi * np.mean(c[2] + 2 * c[3] + 2 * c[4])

        def compress(self):
            self.am = np.asarray([self.am[-1, :], ])
            self.ac = self.ac[-1]
            self.pm1 = np.asarray([self.pm1[-1, :], ])
            self.pm2s = np.asarray([self.pm2s[-1, :], ])
            self.pm2d = np.asarray([self.pm2d[-1, :], ])
            self.pc1 = self.pc1[-1]
            self.pc2 = self.pc2[-1]

            self.aco = np.asarray([self.aco[-1, :], ])
            self.pco = np.asarray([self.pco[-1, :], ])


p0 = Params(Da=1, kon_a=1, koff_a=1, ra=1, Dp=1, kon_p=1, koff_p=1, kon_p_2=5, kd_f=2, kd_b=1, rp=1, L=50, xsteps=500,
            psi=0.3, Tmax=100, deltat=0.01, starts=[0, 1, 0, 0, 0, 1, 0])
