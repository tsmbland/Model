import numpy as np


class Model3:
    """
    pPAR dimerisation model

    """

    class Params:
        def __init__(self, Da, kon_a, koff_a, ra, Dp, kon_p, koff_p, kon_p_2, kd_f, kd_b, rp, pgen):
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

            self.L = pgen.L  # um
            self.xsteps = pgen.xsteps
            self.psi = pgen.psi  # um-1
            self.Tmax = pgen.Tmax  # s
            self.deltat = pgen.deltat  # s

            ######################

    def __init__(self, params, starts):
        self.params = params
        self.am = starts[0] * np.ones([self.params.xsteps])
        self.ac = starts[1]
        self.pm1 = starts[2] * np.ones([self.params.xsteps])
        self.pm2s = starts[3] * np.ones([self.params.xsteps])
        self.pm2d = starts[4] * np.ones([self.params.xsteps])
        self.pc1 = starts[5]
        self.pc2 = starts[6]

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

    def update_am(self, lb, ub):
        """
        Cortical aPAR
        """
        diff = diffusion(self.am[lb:ub], self.params.Da, self.params)
        r17 = self.params.kon_a * self.ac
        r18 = - self.params.koff_a * self.am[lb:ub]
        r19 = - self.params.ra * (self.pm1 + self.pm2s + self.pm2d) * self.am
        self.am[lb:ub] += (diff + r17 + r18 + r19) * self.params.deltat

    def update_ac(self, lb, ub):
        """
        Cytoplasmic aPAR
        """
        r17 = - (1 / self.params.psi) * self.params.kon_a * self.ac
        r18 = (1 / self.params.psi) * np.mean(self.params.koff_a * self.am[lb:ub])
        r19 = (1 / self.params.psi) * np.mean(self.params.ra * (self.pm1 + self.pm2s + self.pm2d) * self.am)
        self.ac += (r17 + r18 + r19) * self.params.deltat

    def update_pm1(self, lb, ub):
        """
        Cortical monomer
        """
        diff = diffusion(self.pm1[lb:ub], self.params.Dp, self.params)
        r1 = self.params.kon_p * self.pc1
        r2 = - self.params.koff_p * self.pm1[lb:ub]
        r3 = - 2 * self.params.kd_f * self.pm1[lb:ub] ** 2
        r4 = 2 * self.params.kd_b * self.pm2d[lb:ub]
        r5 = - self.params.kd_f * self.pc1 * self.pm1[lb:ub]
        r6 = self.params.kd_b * self.pm2s[lb:ub]
        r7 = - self.params.rp * self.am[lb:ub] * self.pm1[lb:ub]
        self.pm1[lb:ub] += (diff + r1 + r2 + r3 + r4 + r5 + r6 + r7) * self.params.deltat

    def update_pm2s(self, lb, ub):
        """
        Cortical dimer (singly bound)
        """
        diff = diffusion(self.pm2s[lb:ub], self.params.Dp, self.params)
        r5 = self.params.kd_f * self.pc1 * self.pm1[lb:ub]
        r6 = - self.params.kd_b * self.pm2s[lb:ub]
        r9 = self.params.kon_p * self.pc2
        r10 = - self.params.koff_p * self.pm2s[lb:ub]
        r11 = - self.params.kon_p_2 * self.pm2s[lb:ub]
        r12 = self.params.koff_p * self.pm2d[lb:ub]
        r13 = - self.params.rp * self.am[lb:ub] * self.pm2s[lb:ub]
        self.pm2s[lb:ub] += (diff + r5 + r6 + r9 + r10 + r11 + r12 + r13) * self.params.deltat

    def update_pm2d(self, lb, ub):
        """
        Cortical dimer (doubly bound)
        """
        diff = diffusion(self.pm2d[lb:ub], self.params.Dp, self.params)
        r3 = self.params.kd_f * self.pm1[lb:ub] ** 2
        r4 = - self.params.kd_b * self.pm2d[lb:ub]
        r11 = self.params.kon_p_2 * self.pm2s[lb:ub]
        r12 = - self.params.koff_p * self.pm2d[lb:ub]
        r15 = - self.params.rp * self.am[lb:ub] * self.pm2d[lb:ub]
        self.pm2d[lb:ub] += (diff + r3 + r4 + r11 + r12 + r15) * self.params.deltat

    def update_pc1(self, lb, ub):
        """
        Cytoplasmic monomer
        """
        r1 = - (1 / self.params.psi) * self.params.kon_p * self.pc1
        r2 = (1 / self.params.psi) * np.mean(self.params.koff_p * self.pm1[lb:ub])
        r5 = - (1 / self.params.psi) * np.mean(self.params.kd_f * self.pc1 * self.pm1[lb:ub])
        r6 = (1 / self.params.psi) * np.mean(self.params.kd_b * self.pm2s[lb:ub])
        r7 = (1 / self.params.psi) * np.mean(self.params.rp * self.am[lb:ub] * self.pm1[lb:ub])
        r21 = - 2 * self.params.kd_f * self.pc1 ** 2
        r22 = 2 * self.params.kd_b * self.pc2
        self.pc1 += (r1 + r2 + r5 + r6 + r7 + r21 + r22) * self.params.deltat

    def update_pc2(self, lb, ub):
        """
        Cytoplasmic dimer
        """
        r9 = - (1 / self.params.psi) * (self.params.kon_p * self.pc2)
        r10 = (1 / self.params.psi) * np.mean(self.params.koff_p * self.pm2s[lb:ub])
        r13 = (1 / self.params.psi) * np.mean(self.params.rp * self.am[lb:ub] * self.pm2s[lb:ub])
        r15 = (1 / self.params.psi) * np.mean(self.params.rp * self.am[lb:ub] * self.pm2d[lb:ub])
        r21 = self.params.kd_f * self.pc1 ** 2
        r22 = - self.params.kd_b * self.pc2
        self.pc2 += (r9 + r10 + r13 + r15 + r21 + r22) * self.params.deltat

    def get_all(self):
        return [self.am, self.ac, self.pm1, self.pm2s, self.pm2d, self.pc1, self.pc2]

    def update_all(self):
        self.update_am(0, self.params.xsteps)
        self.update_ac(0, self.params.xsteps)
        self.update_pm1(0, self.params.xsteps)
        self.update_pm2s(0, self.params.xsteps)
        self.update_pm2d(0, self.params.xsteps)
        self.update_pc1(0, self.params.xsteps)
        self.update_pc2(0, self.params.xsteps)

    def update_all_e(self):
        self.update_am(0, self.params.xsteps // 2)
        self.update_ac(0, self.params.xsteps // 2)
        self.update_pm1(self.params.xsteps // 2, self.params.xsteps)
        self.update_pm2s(self.params.xsteps // 2, self.params.xsteps)
        self.update_pm2d(self.params.xsteps // 2, self.params.xsteps)
        self.update_pc1(self.params.xsteps // 2, self.params.xsteps)
        self.update_pc2(self.params.xsteps // 2, self.params.xsteps)

    class Res:
        def __init__(self, params):
            self.params = params
            self.am = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.ac = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pm1 = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pm2s = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pm2d = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pc1 = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pc2 = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

        def update(self, t, res):
            self.am[t + 1, :] = res[0]
            self.ac[t + 1] = res[1]
            self.pm1[t + 1, :] = res[2]
            self.pm2s[t + 1, :] = res[3]
            self.pm2d[t + 1, :] = res[4]
            self.pc1[t + 1] = res[5]
            self.pc2[t + 1] = res[6]


def diffusion(concs, coeff, p):
    diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
        np.append([1], np.array(range(len(concs) - 1)))]) / (p.L / p.xsteps)

    return diff