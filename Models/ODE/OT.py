import numpy as np
from scipy.integrate import odeint
from Funcs import ParamSpaceQual2D
import matplotlib.pyplot as plt

"""
PARAMETER SETS

JILKINE 2011

a1   = 25
a2   = 0.7
s    = 1
p0   = 2


HUBATSCH 2019

a1   = 1
a2   = 0.7
s    = 1
p0   = ?


"""


###############################################################################

class OT:
    def __init__(self, a1, a2, s, p0):
        self.a1 = a1
        self.a2 = a2
        self.s = s
        self.p0 = p0

    def dUdt(self, X, t):
        U = X[0]
        V = self.p0 - U
        dUdt = self.a1 * (V - (U + V) / ((self.a2 * self.s * (U + V) + 1) ** 2))
        return [dUdt]

    def dUdU(self, X, step=0.001):
        U = X[0]
        V = self.p0 - U
        dUdU = self.a1 * (V - ((U + step) + V) / ((self.a2 * self.s * ((U + step) + V) + 1) ** 2))
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


###############################################################################


# import time
#
# t = time.time()
#
#
# def evaluate1(p0, a2):
#     return OT(a1=1, a2=a2, s=1, p0=p0).bistability_instability()
#
#
# a = Bifurcation2D(evaluate1, p1_range=(0, 5), p2_range=(0, 5), log=False, cores=4, resolution0=50,
#                   resolution_step=2, n_iterations=5, direc='_test2', parallel=False, crange=[1, 6])
# a.run()
#
# print(time.time() - t)
