import numpy as np
import matplotlib.pyplot as plt

"""

"""


class Model:
    def __init__(self, Dm=0.1, Dc=10, k0=0.067, lamda=1, K=1, delta=1, xsteps=100, Tmax=1000, deltat=0.01, deltax=0.1,
                 m_0=..., c_0=...):
        # Species
        self.m = m_0
        self.c = c_0
        self.time = 0

        # Diffusion
        self.Dm = Dm
        self.Dc = Dc

        # Membrane exchange
        self.k0 = k0
        self.lamda = lamda
        self.K = K
        self.delta = delta

        # Misc
        self.xsteps = int(xsteps)
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s
        self.deltax = deltax  # um

    def diffusion(self, concs):
        d = concs[np.append(np.array(range(1, len(concs))), [len(concs) - 1])] - 2 * concs + concs[
            np.append([0], np.array(range(len(concs) - 1)))]
        return d / (self.deltax ** 2)

    def update(self):
        f = self.c * (self.k0 + ((self.lamda * (self.m ** 2)) / ((self.K ** 2) + (self.m ** 2)))) - (
            self.delta * self.m)
        self.m += (self.Dm * self.diffusion(self.m) + f) * self.deltat
        self.c += (self.Dc * self.diffusion(self.c) - f) * self.deltat

    def run(self):
        for t in range(int(self.Tmax / self.deltat)):
            self.update()
            self.time = (t + 1) * self.deltat

    def save(self, direc):
        np.savetxt(direc + 'm.txt', self.m)
        np.savetxt(direc + 'c.txt', self.c)
