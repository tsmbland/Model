import numpy as np
from scipy.integrate import odeint

###############################################################################


"""
Single species dimerisation model

"""


class Model:
    def __init__(self, kon, kon2, koff, kd_f, kd_b, psi, dosage):
        # Dosage
        self.dosage = dosage

        # Membrane exchange
        self.kon = kon
        self.kon2 = kon2
        self.koff = koff

        # Dimerisation
        self.kd_f = kd_f
        self.kd_b = kd_b

        # Misc
        self.psi = psi  # um-1

    def dxdt(self, X, t):
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

        pm1 = X[0]
        pm2s = X[1]
        pm2d = X[2]
        pc1 = X[3]
        pc2 = X[4]

        r = [None] * 13
        r[1] = self.kon * pc1
        r[2] = self.koff * pm1
        r[3] = self.kd_f * pm1 ** 2
        r[4] = self.kd_b * pm2d
        r[5] = self.kd_f * pc1 * pm1
        r[6] = self.kd_b * pm2s
        r[7] = self.kon * pc2
        r[8] = self.koff * pm2s
        r[9] = self.kon2 * pm2s
        r[10] = self.koff * pm2d
        r[11] = self.kd_f * pc1 ** 2
        r[12] = self.kd_b * pc2

        dpm1 = (r[1] - r[2] - 2 * r[3] + 2 * r[4] - r[5] + r[6])
        dpm2s = (r[5] - r[6] + r[7] - r[8] - r[9] + r[10])
        dpm2d = (r[3] - r[4] + r[9] - r[10])
        dpc1 = (- self.psi * r[1] + self.psi * r[2] - self.psi * r[5] + self.psi * r[6] - 2 * r[11] + 2 * r[12])
        dpc2 = (- self.psi * r[7] + self.psi * r[8] + r[11] - r[12])

        return [dpm1, dpm2s, dpm2d, dpc1, dpc2]

    def solve(self):
        sol = odeint(self.dxdt, (0, 0, 0, self.dosage, 0), t=np.linspace(0, 1000, 10000))
        return sol[-1]

