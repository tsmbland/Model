import numpy as np
import matplotlib.pyplot as plt

"""

"""


class Model:
    def __init__(self, Dm, kon, koff, kpos, e, xsteps, Tmax, deltat, deltax, m_0, C):
        # Species
        self.C = C
        self.c = C
        self.m = m_0
        self.time = 0

        # Diffusion
        self.Dm = Dm

        # Membrane exchange
        self.kon = kon
        self.koff = koff
        self.kpos = kpos
        self.e = e

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
        on = self.kon * self.c
        off = self.koff * self.m
        pf = self.kpos * (self.m ** self.e) * self.c
        f = on - off + pf
        m_old = self.m
        self.m += (self.Dm * self.diffusion(self.m) + f) * self.deltat
        self.c = self.C - np.mean(m_old)

    def run(self):
        for t in range(int(self.Tmax / self.deltat)):
            self.update()
            self.time = (t + 1) * self.deltat

    def save(self, direc):
        np.savetxt(direc + 'm.txt', self.m)
        np.savetxt(direc + 'c.txt', [self.c])
