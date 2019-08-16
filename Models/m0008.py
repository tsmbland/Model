import numpy as np
import pandas as pd

"""
pPAR dimerisation model

"""


class Model:
    def __init__(self, Da, kon_a, koff_a, ra, Dp, kon_p, koff_p, kon_p_2, kd_f, kd_b, rp, xsteps, psi, Tmax,
                 deltat, deltal, radii, am_0, ac_0, pm1_0, pm2s_0, pm2d_0, pc1_0, pc2_0, spatial=True):
        ######### Species ###########
        self.am = am_0
        self.ac = ac_0
        self.pm1 = pm1_0
        self.pm2s = pm2s_0
        self.pm2d = pm2d_0
        self.pc1 = pc1_0
        self.pc2 = pc2_0
        self.time = 0

        ########### aPARS ###########

        # Diffusion
        self.Da = Da / deltal

        # Membrane exchange
        self.kon_a = kon_a
        self.koff_a = koff_a

        # Antagonism
        self.ra = ra

        ########### pPARS ###########

        # Diffusion
        self.Dp = Dp / deltal

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

        self.spatial = spatial
        self.xsteps = int(xsteps)
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s
        self.deltal = deltal  # um
        self.radii = radii  # um
        self.L = 500

    def diffusion(self, concs):
        return concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]

    def reactions(self, ):
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
        r[1] = self.kon_p * self.pc1
        r[2] = self.koff_p * self.pm1
        r[3] = self.kd_f * self.pm1 ** 2
        r[4] = self.kd_b * self.pm2d
        r[5] = self.kd_f * self.pc1 * self.pm1
        r[6] = self.kd_b * self.pm2s
        r[7] = self.rp * self.am * self.pm1
        r[9] = self.kon_p * self.pc2
        r[10] = self.koff_p * self.pm2s
        r[11] = self.kon_p_2 * self.pm2s
        r[12] = self.koff_p * self.pm2d
        r[13] = self.rp * self.am * self.pm2s
        r[15] = self.rp * self.am * self.pm2d
        r[21] = self.kd_f * self.pc1 ** 2
        r[22] = self.kd_b * self.pc2
        r[17] = self.kon_a * self.ac
        r[18] = self.koff_a * self.am
        r[19] = self.ra * (self.pm1 + self.pm2s + self.pm2d) * self.am

        if self.spatial:
            r[24] = self.Da * self.diffusion(self.am)
            r[25] = self.Dp * self.diffusion(self.pm1)
            r[26] = self.Dp * self.diffusion(self.pm2s)
            r[27] = self.Dp * self.diffusion(self.pm2d)
        else:
            r[24] = 0
            r[25] = 0
            r[26] = 0
            r[27] = 0
        return r

    def update_am(self, r):
        self.am += (r[24] + r[17] - r[18] - r[19]) * self.deltat

    def update_ac(self, r):
        self.ac += (- (self.psi * r[17]) + (self.psi * np.mean(r[18])) + (self.psi * np.mean(r[19]))) * self.deltat

    def update_pm1(self, r):
        self.pm1 += (r[25] + r[1] - r[2] - 2 * r[3] + 2 * r[4] - r[5] + r[6] - r[7]) * self.deltat

    def update_pm2s(self, r):
        self.pm2s += (r[26] + r[5] - r[6] + r[9] - r[10] - r[11] + r[12] - r[13]) * self.deltat

    def update_pm2d(self, r):
        self.pm2d += (r[27] + r[3] - r[4] + r[11] - r[12] - r[15]) * self.deltat

    def update_pc1(self, r):
        self.pc1 += ((- self.psi * r[1]) + (self.psi * np.mean(r[2])) - (self.psi * np.mean(r[5])) + (
            self.psi * np.mean(r[6])) + (self.psi * np.mean(r[7])) - (2 * r[21]) + (2 * r[22])) * self.deltat

    def update_pc2(self, r):
        self.pc2 += (- (self.psi * r[9]) + (self.psi * np.mean(r[10])) + (self.psi * np.mean(r[13])) + (
            self.psi * np.mean(r[15])) + r[21] - r[22]) * self.deltat

    def react(self):
        r = self.reactions()
        self.update_am(r)
        self.update_ac(r)
        self.update_pm1(r)
        self.update_pm2s(r)
        self.update_pm2d(r)
        self.update_pc1(r)
        self.update_pc2(r)

    def run(self):
        for t in range(int(self.Tmax / self.deltat)):
            self.react()
            self.time = (t + 1) * self.deltat

    def save(self, direc):
        np.savetxt(direc + 'am.txt', self.am)
        np.savetxt(direc + 'ac.txt', [self.ac])
        np.savetxt(direc + 'pm1.txt', [self.pm1])
        np.savetxt(direc + 'pm2s.txt', self.pm2s)
        np.savetxt(direc + 'pm2d.txt', self.pm2d)
        np.savetxt(direc + 'pc1.txt', [self.pc1])
        np.savetxt(direc + 'pc2.txt', [self.pc2])
        np.savetxt(direc + 'time.txt', [self.time])


# def func()
#     r = [None] * 28
#     r[1] = self.kon_p * self.pc1
#     r[2] = self.koff_p * self.pm1
#     r[3] = self.kd_f * self.pm1 ** 2
#     r[4] = self.kd_b * self.pm2d
#     r[5] = self.kd_f * self.pc1 * self.pm1
#     r[6] = self.kd_b * self.pm2s
#     r[7] = self.rp * self.am * self.pm1
#     r[9] = self.kon_p * self.pc2
#     r[10] = self.koff_p * self.pm2s
#     r[11] = self.kon_p_2 * self.pm2s
#     r[12] = self.koff_p * self.pm2d
#     r[13] = self.rp * self.am * self.pm2s
#     r[15] = self.rp * self.am * self.pm2d
#     r[21] = self.kd_f * self.pc1 ** 2
#     r[22] = self.kd_b * self.pc2
#     r[17] = self.kon_a * self.ac
#     r[18] = self.koff_a * self.am
#     r[19] = self.ra * (self.pm1 + self.pm2s + self.pm2d) * self.am
#
#     self.pm1 += (r[1] - r[2] - 2 * r[3] + 2 * r[4] - r[5] + r[6] - r[7])
#     self.pm2s += (r[5] - r[6] + r[9] - r[10] - r[11] + r[12] - r[13])
#     self.pm2d += (r[3] - r[4] + r[11] - r[12] - r[15])
#     self.pc1 += ((- self.psi * r[1]) + (self.psi * np.mean(r[2])) - (self.psi * np.mean(r[5])) + (
#         self.psi * np.mean(r[6])) + (self.psi * np.mean(r[7])) - (2 * r[21]) + (2 * r[22]))
#     self.pc2 += (- (self.psi * r[9]) + (self.psi * np.mean(r[10])) + (self.psi * np.mean(r[13])) + (
#         self.psi * np.mean(r[15])) + r[21] - r[22])
