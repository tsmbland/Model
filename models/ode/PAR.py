import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc + '/../..')

import numpy as np
from scipy.integrate import odeint


###############################################################################


class PAR:
    def __init__(self, konA, koffA, kposA, konP, koffP, kposP, ePneg, eAneg, psi, pA, pP, kAP, kPA, eApos=1, ePpos=1):

        self.konA = konA
        self.koffA = koffA
        self.kposA = kposA
        self.konP = konP
        self.koffP = koffP
        self.kposP = kposP
        self.ePneg = ePneg
        self.eAneg = eAneg
        self.psi = psi
        self.pA = pA
        self.pP = pP
        self.kAP = kAP
        self.kPA = kPA
        self.eApos = eApos
        self.ePpos = ePpos

    def dxdt(self, X, t):
        A = X[0]
        P = X[1]
        Acyt = self.pA - self.psi * A
        Pcyt = self.pP - self.psi * P
        dA = (self.konA * Acyt) - (self.koffA * A) - (self.kAP * (P ** self.ePneg) * A) + (
                self.kposA * (A ** self.eApos) * Acyt)
        dP = (self.konP * Pcyt) - (self.koffP * P) - (self.kPA * (A ** self.eAneg) * P) + (
                self.kposP * (P ** self.ePpos) * Pcyt)
        return [dA, dP]

    def jacobian(self, X, step=0.0001):
        A = X[0]
        P = X[1]
        Acyt = self.pA - self.psi * A
        Pcyt = self.pP - self.psi * P
        dPdA = (self.konP * Pcyt) - (self.koffP * P) - (self.kPA * ((A + step) ** self.eAneg) * P) + (
                self.kposP * (P ** self.ePpos) * Pcyt)
        dAdP = (self.konA * Acyt) - (self.koffA * A) - (self.kAP * ((P + step) ** self.ePneg) * A) + (
                self.kposA * (A ** self.eApos) * Acyt)
        dAdA = (self.konA * Acyt) - (self.koffA * (A + step)) - (self.kAP * (P ** self.ePneg) * (A + step)) + (
                self.kposA * ((A + step) ** self.eApos) * Acyt)
        dPdP = (self.konP * Pcyt) - (self.koffP * (P + step)) - (self.kPA * (A ** self.eAneg) * (P + step)) + (
                self.kposP * ((P + step) ** self.ePpos) * Pcyt)
        return np.r_[np.c_[dAdA, dAdP], np.c_[dPdA, dPdP]] / step

    def bistability(self):
        """
        Bistability - solved using odeint

        0: monostable
        1: bistable

        """

        sol1 = odeint(self.dxdt, (self.pA / self.psi, 0), t=np.linspace(0, 10000, 100000))[-1]
        sol2 = odeint(self.dxdt, (0, self.pP / self.psi), t=np.linspace(0, 10000, 100000))[-1]

        # Distance measure
        dist = (((sol1[0] - sol2[0]) ** 2) + ((sol1[1] - sol2[1]) ** 2)) ** 0.5

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

        sol1 = odeint(self.dxdt, (self.pA / self.psi, 0), t=np.linspace(0, 10000, 100000))[-1]
        sol2 = odeint(self.dxdt, (0, self.pP / self.psi), t=np.linspace(0, 10000, 100000))[-1]

        # Distance measure
        dist = (((sol1[0] - sol2[0]) ** 2) + ((sol1[1] - sol2[1]) ** 2)) ** 0.5

        if dist < 0.0001:
            # Monostable
            point = (sol1 + sol2) / 2
            w, v = np.linalg.eig(self.jacobian(point, step=0.0001))
            if np.all(w < 0):
                # Stable
                return 1
            else:
                # Unstable
                return 2

        else:
            # Bistable
            w1, v1 = np.linalg.eig(self.jacobian(sol1, step=0.0001))
            w2, v2 = np.linalg.eig(self.jacobian(sol2, step=0.0001))
            if np.all(w1 < 0) and np.all(w2 < 0):
                # Both stable
                return 3
            elif not np.all(w1 < 0) and np.all(w2 < 0):
                # sol1 unstable
                return 4
            elif np.all(w1 < 0) and not np.all(w2 < 0):
                # sol2 unstable
                return 5
            else:
                # Both unstable
                return 6


###############################################################################


"""
PARAMETER SETS

GOEHRING 2011 symmetric (misreported) / TRONG 2014

konA    = 1
koffA   = 0.3
konP    = 1
koffP   = 0.3
ePneg   = 2
eAneg   = 2
kAP     = 1
kPA     = 1
psi     = 0.3
pA      = 1
pP      = 1


GOEHRING 2011 measured

konA    = 0.00858
koffA   = 0.0054
konP    = 0.0474
koffP   = 0.0073
ePneg   = 1
eAneg   = 2
kAP     = 0.19
kPA     = 2
psi     = 0.174
pA      = 1.56
pP      = 1


GROSS 2018

konA    = 0.02115
koffA   = 0.0092
konP    = 0.13012
koffP   = 0.0073
ePneg   = -
eAneg   = -
kAP     = -
kPA     = -
psi     = -
pA      = -
pP      = -


HUBATSCH 2019 

konA    = 0.006
koffA   = 0.005
konP    = 0.006
koffP   = 0.005
ePneg   = 2
eAneg   = 2
kAP     = 1
kPA     = 1
psi     = variable
pA      = ?
pP      = ?



"""

# import time
#
# t = time.time()
#
#
# def evaluate1(kAP, kPA):
#     return PAR(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, ePneg=2, eAneg=2, psi=0.1,
#                pA=1, pP=1, kAP=kAP, kPA=kPA).bistability_instability()
#

#
# def evaluate2(pA, pP):
#     return PAR(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, ePneg=2, eAneg=2, psi=0.1,
#                pA=pA, pP=pP, kAP=0.01, kPA=0.01).bistability_instability()
#
#
# def evaluate3(pP, kAP):
#     return PAR(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, ePneg=2, eAneg=2, psi=0.1,
#                pA=1, pP=pP, kAP=kAP, kPA=0.01).bistability_instability()
#
#
# def evaluate4(pP, kPA):
#     return PAR(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, ePneg=2, eAneg=2, psi=0.1,
#                pA=1, pP=pP, kAP=0.01, kPA=kPA).bistability_instability()
#
#
# def evaluate5(pP, k):
#     return PAR(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, ePneg=2, eAneg=2, psi=0.1,
#                pA=1, pP=pP, kAP=k, kPA=k).bistability_instability()
#
#
# a = ParamSweep2D(evaluate1, p1_range=(0, 0.01), p2_range=(0, 0.01), log=False, cores=4, resolution0=50,
#                  resolution_step=2, n_iterations=4, direc='_test2', parallel=True, crange=[1, 6])
# a.run()
#
# print(time.time() - t)

# import time
#
#
# def evaluate1(pP, kposP):
#     konP0 = 0.1
#     y0 = (0.1 * 1) / (0.1 * 0.1 + 0.01)
#     konP = konP0 - (kposP * y0)
#
#     return PAR(konA=0.1, koffA=0.01, kposA=0, konP=konP, koffP=0.01, kposP=kposP, ePneg=1, eAneg=1, psi=0.1,
#                pA=1, pP=pP, kAP=0.1, kPA=0.1).bistability_instability()
#
#

# from ModelFuncs import ParamSpace2D
#
# t = time.time()
# a = ParamSpace2D(evaluate1, p1_range=(0, 5), p2_range=(0, 0.02), cores=4, resolution0=50,
#                  resolution_step=2, n_iterations=2, direc='_test', parallel=True, crange=[1, 6])
# a.run()
# print(time.time() - t)

#
# t = time.time()
# a = Bifurcation2D(evaluate1, p1_range=(0, 5), p2_range=(0, 0.02), cores=4, resolution0=50,
#                   resolution_step=2, n_iterations=5, direc='_test2', parallel=True, crange=[1, 6])
# a.run()
# print(time.time() - t)
