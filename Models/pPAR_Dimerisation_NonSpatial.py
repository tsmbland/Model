import numpy as np

"""
pPAR dimerisation model

No aPAR

"""


class Model:
    def __init__(self, kon_p, kon_p2, koff_p, kd_f, kd_b, psi, Tmax, deltat, pm1_0, pm2s_0, pm2d_0, pc1_0, pc2_0):
        # Start concs
        self.pm1_0 = pm1_0
        self.pm2s_0 = pm2s_0
        self.pm2d_0 = pm2d_0
        self.pc1_0 = pc1_0
        self.pc2_0 = pc2_0

        # Species
        self.pm1 = pm1_0
        self.pm2s = pm2s_0
        self.pm2d = pm2d_0
        self.pc1 = pc1_0
        self.pc2 = pc2_0

        self.pco = 0
        self.pcy = 0
        self.ptot = 0
        self.time = 0
        self.time = 0

        # Membrane exchange
        self.kon_p = kon_p
        self.kon_p2 = kon_p2
        self.koff_p = koff_p

        # Dimerisation
        self.kd_f = kd_f
        self.kd_b = kd_b

        # Misc
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s

    def reactions(self):
        """
        r1: binding of monomers to cortex
        r2: unbinding of monomers from cortex
        r3: Dimerisation of pm1
        r4: Separation of pm2d
        r5: Dimerisation of pm1 with pc1
        r6: Separation of pm2s
        r7: Single binding of pc2 to cortex
        r8: unbinding of pm2s from cortex
        r9: binding of second component of dimer to cortex
        r10: unbinding of one component of pmd2
        r11: dimerisation of pc1
        r12: separation of pc2

        """

        r = [None] * 13
        r[1] = self.kon_p * self.pc1
        r[2] = self.koff_p * self.pm1
        r[3] = self.kd_f * self.pm1 ** 2
        r[4] = self.kd_b * self.pm2d
        r[5] = self.kd_f * self.pc1 * self.pm1
        r[6] = self.kd_b * self.pm2s
        r[7] = self.kon_p * self.pc2
        r[8] = self.koff_p * self.pm2s
        r[9] = self.kon_p2 * self.pm2s
        r[10] = self.koff_p * self.pm2d
        r[11] = self.kd_f * self.pc1 ** 2
        r[12] = self.kd_b * self.pc2
        return r

    def update_pm1(self, r):
        self.pm1 += (r[1] - r[2] - 2 * r[3] + 2 * r[4] - r[5] + r[6]) * self.deltat

    def update_pm2s(self, r):
        self.pm2s += (r[5] - r[6] + r[7] - r[8] - r[9] + r[10]) * self.deltat

    def update_pm2d(self, r):
        self.pm2d += (r[3] - r[4] + r[9] - r[10]) * self.deltat

    def update_pc1(self, r):
        self.pc1 += (- self.psi * r[1] + self.psi * r[2] - self.psi * r[5] + self.psi * r[6] - 2 * r[11] + 2 * r[
            12]) * self.deltat

    def update_pc2(self, r):
        self.pc2 += (- self.psi * r[7] + self.psi * r[8] + r[11] - r[12]) * self.deltat

    def pool(self):
        self.pco = self.pm1 + 2 * self.pm2s + 2 * self.pm2d
        self.pcy = self.pc1 + 2 * self.pc2
        self.ptot = self.pcy + self.psi * self.pco

    def react(self):
        r = self.reactions()
        self.update_pm1(r)
        self.update_pm2s(r)
        self.update_pm2d(r)
        self.update_pc1(r)
        self.update_pc2(r)

    def save(self, direc):
        np.savetxt(direc + '/pm1.txt', self.pm1)
        np.savetxt(direc + '/pm2s.txt', self.pm2s)
        np.savetxt(direc + '/pm2d.txt', self.pm2d)
        np.savetxt(direc + '/pc1.txt', self.pc1)
        np.savetxt(direc + '/pc2.txt', self.pc2)

        np.savetxt(direc + '/pco.txt', self.pco)
        np.savetxt(direc + '/pcy.txt', self.pcy)
        np.savetxt(direc + '/ptot.txt', self.ptot)
