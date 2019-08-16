import numpy as np
import matplotlib.pyplot as plt

"""
Positive feedback model, with saturation


"""


class Model:
    def __init__(self, Da, Dp, konA, koffA, konP, koffP, kposP, kposP2, kAP, kPA, ePneg, eAneg, xsteps, psi, Tmax,
                 deltat, deltal, radii, am_0, ac_0, pm_0, pc_0):
        # Species
        self.am = am_0
        self.ac = ac_0
        self.pm = pm_0
        self.pc = pc_0
        self.time = 0

        # Diffusion
        self.Da = Da  # um2 s-1
        self.Dp = Dp  # um2 s-1

        # Membrane exchange
        self.konA = konA  # um s-1
        self.koffA = koffA  # s-1

        self.konP = konP  # um s-1
        self.koffP = koffP  # s-1
        self.kposP = kposP  # um3 s-1
        self.kposP2 = kposP2

        # Antagonism
        self.kAP = kAP  # um2 s-1
        self.kPA = kPA  # um4 s-1
        self.ePneg = ePneg
        self.eAneg = eAneg

        # Misc
        self.xsteps = int(xsteps)
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s
        self.deltal = deltal  # um
        self.radii = radii  # um

    def diffusion(self, concs):
        return (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.deltal ** 2)

    def reactions(self):
        """
        r0: a on
        r1: a off
        r2: p to a antagnism
        r3: a diffusion

        r4: p on
        r5: p off
        r6: a to p antagonism
        r7: p diffusion

        r8: pos feedback

        """

        r = [None] * 9

        r[0] = self.konA * self.ac
        r[1] = self.koffA * self.am
        r[2] = self.kAP * (self.pm ** self.ePneg) * self.am
        r[3] = self.Da * self.diffusion(self.am)

        r[4] = (self.konP * self.pc)
        r[5] = self.koffP * self.pm
        r[6] = self.kPA * (self.am ** self.eAneg) * self.pm
        r[7] = self.Dp * self.diffusion(self.pm)

        r[8] = self.kposP * self.pc * self.pm * np.exp(-self.kposP2 * (self.pm ** 2))

        return r

    def update_am(self, r):
        self.am += (r[0] - r[1] - r[2] + r[3]) * self.deltat

    def update_pm(self, r):
        self.pm += (r[4] - r[5] - r[6] + r[7] + r[8]) * self.deltat

    def update_ac(self, r):
        self.ac += (- self.psi * r[0] + self.psi * np.average(r[1], weights=self.radii) + self.psi * np.average(r[2],
                                                                                                                weights=self.radii)) * self.deltat

    def update_pc(self, r):
        self.pc += (- self.psi * r[4] + self.psi * np.average(r[5], weights=self.radii) + self.psi * np.average(r[6],
                                                                                                                weights=self.radii) - self.psi * np.average(
            r[8], weights=self.radii)) * self.deltat

    def react(self):
        r = self.reactions()
        self.update_am(r)
        self.update_ac(r)
        self.update_pm(r)
        self.update_pc(r)

    def run(self):
        for t in range(int(self.Tmax / self.deltat)):
            self.react()
            self.time = (t + 1) * self.deltat

    def save(self, direc):
        np.savetxt(direc + 'ac.txt', [self.ac])
        np.savetxt(direc + 'am.txt', self.am)
        np.savetxt(direc + 'pc.txt', [self.pc])
        np.savetxt(direc + 'pm.txt', self.pm)
        np.savetxt(direc + 'time.txt', [self.time])
