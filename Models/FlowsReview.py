import numpy as np
import matplotlib.pyplot as plt

"""
Looking at flows, diffusion rates, off rates

Changes:
- linear velocity gradient
- explicitly model cytoplasmic species (same a membrane, just faster diffusion)
- change koff but keep kon/koff constant

"""


class Model:
    def __init__(self, Dm, Dc, Vm, Vc, kon, koff, xsteps, Tmax, deltat, deltax, c_0, m_0, flowtype):
        # Species
        self.m = np.ones([xsteps]) * m_0
        self.c = np.ones([xsteps]) * c_0
        self.time = 0

        # Diffusion
        self.Dm = Dm
        self.Dc = Dc

        # Flow
        self.Vm = Vm
        self.Vc = Vc
        self.flowtype = flowtype

        # Membrane exchange
        self.kon = kon
        self.koff = koff

        # Misc
        self.xsteps = int(xsteps)
        self.Tmax = Tmax
        self.deltat = deltat
        self.deltal = deltax

    def diffusion(self, concs):
        return concs[np.append(np.array(range(1, len(concs))), [len(concs) - 1])] - 2 * concs + concs[
            np.append([0], np.array(range(len(concs) - 1)))]

    def flow(self, concs):

        if self.flowtype == 1:
            # Realistic
            x = np.array(range(self.xsteps)) * (100 / self.xsteps)
            v = (x / np.exp(0.00075 * (x ** 2)))[::-1]
            v[0] = 0
            return - np.diff(np.r_[concs, concs[-1]] * np.r_[v, 0])

        if self.flowtype == 2:
            # Realistic 2
            x = np.array(range(self.xsteps)) * (100 / self.xsteps)
            v = (x / np.exp(0.08 * (x ** 1)))[::-1]
            v[0] = 0
            return - np.diff(np.r_[concs, concs[-1]] * np.r_[v, 0])

        if self.flowtype == 3:
            # Linear
            v = np.arange(self.xsteps) / self.xsteps
            return - np.diff(np.r_[concs, concs[-1]] * np.r_[v, 0])

    def reactions(self):
        """
        r0: on
        r1: off
        r2: mem diffusion
        r3: cyt diffusion
        r4: mem flow
        r5: cyt flow
        """
        r = [None] * 6
        r[0] = self.kon * self.c
        r[1] = self.koff * self.m
        r[2] = (self.Dm / (self.deltal ** 2)) * self.diffusion(self.m)
        r[3] = (self.Dc / (self.deltal ** 2)) * self.diffusion(self.c)
        r[4] = (self.Vm / self.deltal) * self.flow(self.m)
        r[5] = (self.Vc / self.deltal) * self.flow(self.c)
        return r

    def update_m(self, r):
        self.m += (r[0] - r[1] + r[2] - r[4]) * self.deltat

    def update_c(self, r):
        self.c += (-r[0] + r[1] + r[3] - r[5]) * self.deltat

    def react(self):
        r = self.reactions()
        self.update_m(r)
        self.update_c(r)

    def save(self, direc):
        np.savetxt(direc + 'c.txt', self.c)
        np.savetxt(direc + 'm.txt', self.m)
        np.savetxt(direc + 'time.txt', [self.time])
