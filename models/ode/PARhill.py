import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc + '/../..')

import numpy as np
from scipy.integrate import odeint


###############################################################################


class PARhill:
    def __init__(self, konA, koffA, kposA, konP, koffP, kposP, ePneg, eAneg, psi, pA, pP, kAP, kPA, khillA, khillP,
                 ehillA=1, ehillP=1):

        self.konA = konA
        self.koffA = koffA
        self.kposA = kposA
        self.konP = konP
        self.koffP = koffP
        self.kposP = kposP
        self.khillA = khillA
        self.khillP = khillP
        self.ehillA = ehillA
        self.ehillP = ehillP
        self.ePneg = ePneg
        self.eAneg = eAneg
        self.psi = psi
        self.pA = pA
        self.pP = pP
        self.kAP = kAP
        self.kPA = kPA

    def dxdt(self, X, t):
        A = X[0]
        P = X[1]
        Acyt = self.pA - self.psi * A
        Pcyt = self.pP - self.psi * P
        dA = (self.konA * Acyt) - (self.koffA * A) - (self.kAP * (P ** self.ePneg) * A) + (
                self.kposA * Acyt * (A ** self.ehillA) / ((self.khillA ** self.ehillA) + (A ** self.ehillA)))
        dP = (self.konP * Pcyt) - (self.koffP * P) - (self.kPA * (A ** self.eAneg) * P) + (
                self.kposP * Pcyt * (P ** self.ehillP) / ((self.khillP ** self.ehillP) + (P ** self.ehillP)))
        return [dA, dP]

    def jacobian(self, X, step=0.0001):
        A = X[0]
        P = X[1]
        Acyt = self.pA - self.psi * A
        Pcyt = self.pP - self.psi * P
        dPdA = (self.konP * Pcyt) - (self.koffP * P) - (self.kPA * ((A + step) ** self.eAneg) * P) + (
                self.kposP * Pcyt * (P ** self.ehillP) / ((self.khillP ** self.ehillP) + (P ** self.ehillP)))
        dAdP = (self.konA * Acyt) - (self.koffA * A) - (self.kAP * ((P + step) ** self.ePneg) * A) + (
                self.kposA * Acyt * (A ** self.ehillA) / ((self.khillA ** self.ehillA) + (A ** self.ehillA)))
        dAdA = (self.konA * Acyt) - (self.koffA * (A + step)) - (self.kAP * (P ** self.ePneg) * (A + step)) + (
                self.kposA * Acyt * ((A + step) ** self.ehillA) / (
                (self.khillA ** self.ehillA) + ((A + step) ** self.ehillA)))
        dPdP = (self.konP * Pcyt) - (self.koffP * (P + step)) - (self.kPA * (A ** self.eAneg) * (P + step)) + (
                self.kposP * Pcyt * ((P + step) ** self.ehillP) / (
                (self.khillP ** self.ehillP) + ((P + step) ** self.ehillP)))
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
                # sol1 (aPAR dominant) unstable
                return 4
            elif np.all(w1 < 0) and not np.all(w2 < 0):
                # sol2 (pPAR dominant) unstable
                return 5
            else:
                # Both unstable
                return 6

###############################################################################
