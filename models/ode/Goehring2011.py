import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc + '/../..')

import numpy as np
from scipy.integrate import odeint


###############################################################################


class Goehring2011:
    def __init__(self, konA=0.00858, koffA=0.0054, konP=0.0474, koffP=0.0073, kPA=2, kAP=0.19,
                 alpha=1, beta=2, psi=0.174, pA=1.56, pP=1):
        self.konA = konA
        self.koffA = koffA
        self.konP = konP
        self.koffP = koffP
        self.alpha = alpha
        self.beta = beta
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
        dA = (self.konA * Acyt) - (self.koffA * A) - (self.kAP * (P ** self.alpha) * A)
        dP = (self.konP * Pcyt) - (self.koffP * P) - (self.kPA * (A ** self.beta) * P)
        return [dA, dP]

    def jacobian(self, X, step=0.0001):
        A = X[0]
        P = X[1]
        Acyt = self.pA - self.psi * A
        Pcyt = self.pP - self.psi * P
        dPdA = (self.konP * Pcyt) - (self.koffP * P) - (self.kPA * ((A + step) ** self.beta) * P)
        dAdP = (self.konA * Acyt) - (self.koffA * A) - (self.kAP * ((P + step) ** self.alpha) * A)
        dAdA = (self.konA * Acyt) - (self.koffA * (A + step)) - (self.kAP * (P ** self.alpha) * (A + step))
        dPdP = (self.konP * Pcyt) - (self.koffP * (P + step)) - (self.kPA * (A ** self.beta) * (P + step))
        return np.r_[np.c_[dAdA, dAdP], np.c_[dPdA, dPdP]] / step
