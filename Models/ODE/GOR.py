import numpy as np
from scipy.integrate import odeint

"""
PARAMETER SETS

JILKINE 2011

a1   = ?
a2   = ?
a3   = 0.01733
p0   = ?


HUBATSCH 2019

a1   = 0.0067
a2   = 0.0033
a3   = 0.01
p0   = ?

"""


###############################################################################

class GOR:
    def __init__(self, a1, a2, a3, p0):
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.p0 = p0

    def dUdt(self, X, t):
        U = X[0]
        V = self.p0 - X[0]
        dU = (self.a1 * (U ** 2) * V) + (self.a2 * U * V) - (self.a3 * U)
        return [dU]

    def dUdU(self, X, step=0.001):
        U = X[0]
        V = self.p0 - X[0]
        dUdU = (self.a1 * ((U + step) ** 2) * V) + (self.a2 * (U + step) * V) - (self.a3 * (U + step))
        return dUdU / step

    def bistability(self):
        """
        Bistability - solved using odeint

        0: monostable
        1: bistable

        """
        sol1 = odeint(self.dUdt, [self.p0], t=np.linspace(0, 10000, 100000))[-1]
        sol2 = odeint(self.dUdt, [0], t=np.linspace(0, 10000, 100000))[-1]

        # Distance measure
        dist = abs(sol1 - sol2)

        if dist > 0.0001:
            return 2
        else:
            return 1

    def bistability_instability(self):
        """
        Bistability + instability regions - solved using odeint

        0: monostable, stable
        1: monostable, unstable
        2: bistable, both stable
        3: bistable, sol1 unstable
        4: bistable, sol2 unstable
        5: bistable, both unstable

        """

        sol1 = odeint(self.dUdt, [self.p0], t=np.linspace(0, 10000, 100000))[-1]
        sol2 = odeint(self.dUdt, [0], t=np.linspace(0, 10000, 100000))[-1]

        # Distance measure
        dist = abs(sol1 - sol2)

        if dist < 0.0001:
            # Monostable
            point = (sol1 + sol2) / 2
            dUdU = self.dUdU(point, step=0.001)
            if dUdU < 0:
                # Stable
                return 1
            else:
                # Unstable
                return 2

        else:
            # Bistable
            dUdU1 = self.dUdU(sol1, step=0.001)
            dUdU2 = self.dUdU(sol2, step=0.001)

            if dUdU1 < 0 and dUdU2 < 0:
                # Both stable
                return 3
            elif not dUdU1 < 0 and dUdU2 < 0:
                # sol1 unstable
                return 4
            elif dUdU1 < 0 and not dUdU2 < 0:
                # sol2 unstable
                return 5
            else:
                # Both unstable
                return 6
