import numpy as np
from scipy.integrate import odeint


class Bistability:
    def __init__(self, konA, koffA, kposA, konP, koffP, kposP, ePneg, eAneg, psi, pA, pP):

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

    def dxdt(self, X, t, kAP, kPA):
        A = X[0]
        P = X[1]
        dA = self.konA * (self.pA - self.psi * A) - (self.koffA * A) - (kAP * (P ** self.ePneg) * A)
        dP = self.konP * (self.pP - self.psi * P) - (self.koffP * P) - (kPA * (A ** self.eAneg) * P)
        return np.array([dA, dP])

    def jacobian(self, X, kAP, kPA, step):
        A = X[0]
        P = X[1]

        dPdA = self.konP * (self.pP - self.psi * P) - (self.koffP * P) - (kPA * ((A + step) ** self.eAneg) * P)
        dAdP = self.konA * (self.pA - self.psi * A) - (self.koffA * A) - (kAP * ((P + step) ** self.ePneg) * A)
        dAdA = self.konA * (self.pA - self.psi * A) - (self.koffA * (A + step)) - (kAP * (P ** self.ePneg) * (A + step))
        dPdP = self.konP * (self.pP - self.psi * P) - (self.koffP * (P + step)) - (kPA * (A ** self.eAneg) * (P + step))

        return np.r_[np.c_[dAdA, dAdP], np.c_[dPdA, dPdP]] / step

    # def func(self, kAP, kPA):
    #     """
    #     Bistability - solved using odeint
    #
    #     0: monostable
    #     1: bistable
    #
    #     """
    #
    #     sol1 = odeint(self.dxdt, (self.pA / self.psi, 0), t=np.linspace(0, 10000, 100000), args=(kAP, kPA))[-1]
    #     sol2 = odeint(self.dxdt, (0, self.pP / self.psi), t=np.linspace(0, 10000, 100000), args=(kAP, kPA))[-1]
    #
    #     # Distance measure
    #     dist = (((sol1[0] - sol2[0]) ** 2) + ((sol1[1] - sol2[1]) ** 2)) ** 0.5
    #
    #     if dist > 0.0001:
    #         return 1
    #     else:
    #         return 0

    def func(self, kAP, kPA):
        """
        Bistability + instability regions - solved using odeint

        0: monostable, stable
        1: monostable, unstable
        2: bistable, both stable
        3: bistable, sol1 unstable
        4: bistable, sol2 unstable
        5: bistable, both unstable

        """

        sol1 = odeint(self.dxdt, (self.pA / self.psi, 0), t=np.linspace(0, 10000, 100000), args=(kAP, kPA))[-1]
        sol2 = odeint(self.dxdt, (0, self.pP / self.psi), t=np.linspace(0, 10000, 100000), args=(kAP, kPA))[-1]

        # Distance measure
        dist = (((sol1[0] - sol2[0]) ** 2) + ((sol1[1] - sol2[1]) ** 2)) ** 0.5

        if dist < 0.0001:
            # Monostable
            point = (sol1 + sol2) / 2
            w, v = np.linalg.eig(self.jacobian(point, kAP, kPA, step=0.001))
            if np.all(w < 0):
                # Stable
                return 0
            else:
                # Unstable
                return 1

        else:
            # Bistable
            w1, v1 = np.linalg.eig(self.jacobian(sol1, kAP, kPA, step=0.001))
            w2, v2 = np.linalg.eig(self.jacobian(sol2, kAP, kPA, step=0.001))
            if np.all(w1 < 0) and np.all(w2 < 0):
                # Both stable
                return 2
            elif not np.all(w1 < 0) and np.all(w2 < 0):
                # sol1 unstable
                return 3
            elif np.all(w1 < 0) and not np.all(w2 < 0):
                # sol2 unstable
                return 4
            else:
                # Both unstable
                return 5
