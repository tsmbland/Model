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

        r = np.zeros([28])
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
        r[19] = self.params.ra * (self.pm1 + self.pm2s + self.pm2d) * self.am
        r[24] = diffusion(self.am[lb:ub], self.params.Da, self.params)
        r[25] = diffusion(self.pm1[lb:ub], self.params.Dp, self.params)
        r[26] = diffusion(self.pm2s[lb:ub], self.params.Dp, self.params)
        r[27] = diffusion(self.pm2d[lb:ub], self.params.Dp, self.params)

        return r

    def update_am(self, r):
        """
        Cortical aPAR
        """
        self.am += (r[24] + r[17] - r[18] - r[19]) * self.params.deltat

    def update_ac(self, r):
        """
        Cytoplasmic aPAR
        """
        self.ac += (- (1 / self.params.psi) * r[17] + (1 / self.params.psi) * np.mean(r[18]) + (
            1 / self.params.psi) * np.mean(r[19])) * self.params.deltat

    def update_pm1(self, r):
        """
        Cortical monomer
        """
        self.pm1 += (r[25] + r[1] - r[2] - 2 * r[3] + 2 * r[4] - r[5] + r[6] - r[7]) * self.params.deltat

    def update_pm2s(self, r):
        """
        Cortical dimer (singly bound)
        """
        self.pm2s += (r[26] + r[5] - r[6] + r[9] - r[10] - r[11] + r[12] - r[13]) * self.params.deltat

    def update_pm2d(self, r):
        """
        Cortical dimer (doubly bound)
        """
        self.pm2d += (r[27] + r[3] - r[4] + r[11] - r[12] - r[15]) * self.params.deltat

    def update_pc1(self, r):
        """
        Cytoplasmic monomer
        """
        self.pc1 += (- (1 / self.params.psi) * r[1] + (1 / self.params.psi) * np.mean(r[2]) - (
            1 / self.params.psi) * np.mean(r[5]) + (1 / self.params.psi) * np.mean(r[6]) + (
                         1 / self.params.psi) * np.mean(
            r[7]) - 2 * r[21] + 2 * r[22]) * self.params.deltat

    def update_pc2(self, r):
        """
        Cytoplasmic dimer
        """
        self.pc2 += (- (1 / self.params.psi) * r[9] + (1 / self.params.psi) * np.mean(r[10]) + (
            1 / self.params.psi) * np.mean(r[13]) + (1 / self.params.psi) * np.mean(r[15]) + r[21] - r[
                         22]) * self.params.deltat

    def get_all(self):
        return [self.am, self.ac, self.pm1, self.pm2s, self.pm2d, self.pc1, self.pc2]

    def update_all(self, r):
        self.update_am(r)
        self.update_ac(r)
        self.update_pm1(r)
        self.update_pm2s(r)
        self.update_pm2d(r)
        self.update_pc1(r)
        self.update_pc2(r)

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

            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

        def update(self, t, res):
            self.am[t + 1, :] = res[0]
            self.ac[t + 1] = res[1]
            self.pm1[t + 1, :] = res[2]
            self.pm2s[t + 1, :] = res[3]
            self.pm2d[t + 1, :] = res[4]
            self.pc1[t + 1] = res[5]
            self.pc2[t + 1] = res[6]

            self.aco[t + 1, :] = res[0]
            self.pco[t + 1, :] = res[2] + res[3] + res[4]


def diffusion(concs, coeff, p):
    diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
        np.append([1], np.array(range(len(concs) - 1)))]) / (p.L / p.xsteps)

    return diff
