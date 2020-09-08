import numpy as np
from scipy.integrate import odeint


###############################################################################


class Goehring2011:
    def __init__(self, konA, koffA, konP, koffP, alpha, beta, psi, pA, pP, kAP, kPA):
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
        dA = (self.konA * Acyt) - (self.koffA * A) - (self.kAP * (P ** self.beta) * A)
        dP = (self.konP * Pcyt) - (self.koffP * P) - (self.kPA * (A ** self.alpha) * P)
        return [dA, dP]
