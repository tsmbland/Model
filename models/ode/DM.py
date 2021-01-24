import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc + '/../..')

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
        self.psi = psi  # depletion factor from the paper

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

        # Species
        pm1 = X[0]
        pm2s = X[1]
        pm2d = X[2]
        pc1 = X[3]
        pc2 = X[4]

        # Reactions
        r = np.zeros(12)
        r[0] = self.kon * pc1
        r[1] = self.koff * pm1
        r[2] = self.kd_f * pm1 ** 2
        r[3] = self.kd_b * pm2d
        r[4] = self.kd_f * pc1 * pm1
        r[5] = self.kd_b * pm2s
        r[6] = self.kon * pc2
        r[7] = self.koff * pm2s
        r[8] = self.kon2 * pm2s
        r[9] = self.koff * pm2d
        r[10] = self.kd_f * pc1 ** 2
        r[11] = self.kd_b * pc2

        # Equations
        dpm1 = r[0] - r[1] - 2 * r[2] + 2 * r[3] - r[4] + r[5]
        dpm2s = r[4] - r[5] + r[6] - r[7] - r[8] + r[9]
        dpm2d = r[2] - r[3] + r[8] - r[9]
        dpc1 = self.psi * (- r[0] + r[1] - r[4] + r[5]) - 2 * r[10] + 2 * r[11]
        dpc2 = self.psi * (- r[6] + r[7]) + r[10] - r[11]
        return [dpm1, dpm2s, dpm2d, dpc1, dpc2]

    def solve(self):
        sol = odeint(self.dxdt, (0, 0, 0, self.dosage, 0), t=np.linspace(0, 1000, 10000))
        return sol[-1]

# import matplotlib.pyplot as plt
#
# for d in np.linspace(0, 1, 100):
#     print(d)
#     m = Model(kon=1, kon2=10, koff=0.1, kd_f=1, kd_b=1, psi=0.1, dosage=d)
#     a = m.solve()
#     m = (a[0] + 2 * a[1] + 2 * a[2])
#     c = (a[3] + 2 * a[4])
#     plt.scatter(c, m)
# plt.show()
