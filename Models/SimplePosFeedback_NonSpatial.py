import numpy as np

"""

"""


class Model:
    def __init__(self, konA, koffA, konP, koffP, kposA, kposP, kAP, kPA, ePneg, eAneg, Tmax, deltat,
                 psi, am_0, ac_0, pm_0, pc_0):
        # Start concs
        self.am_0 = am_0
        self.ac_0 = ac_0
        self.pm_0 = pm_0
        self.pc_0 = pc_0

        # Species
        self.am = am_0
        self.ac = ac_0
        self.pm = pm_0
        self.pc = pc_0
        self.time = 0

        # Membrane exchange
        self.konA = konA  # um s-1
        self.koffA = koffA  # s-1

        self.konP = konP  # um s-1
        self.koffP = koffP  # s-1

        # Positive feedback
        self.kposA = kposA
        self.kposP = kposP

        # Antagonism
        self.kAP = kAP
        self.kPA = kPA
        self.ePneg = ePneg
        self.eAneg = eAneg

        # Misc
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s
        self.psi = psi  # um-1

    def reactions(self):
        """
        r0: a on
        r1: a off
        r2: a positive feedback
        r3: p to a antagonism

        r5: p on
        r6: p off
        r7: p positive feedback
        r8: a to p antagonism

        """

        r = [None] * 10

        r[0] = self.konA * self.ac
        r[1] = self.koffA * self.am
        r[2] = self.kposA * self.am * self.ac
        r[3] = self.kAP * (self.pm ** self.ePneg) * self.am

        r[5] = self.konP * self.pc
        r[6] = self.koffP * self.pm
        r[7] = self.kposP * self.pm * self.pc
        r[8] = self.kPA * (self.am ** self.eAneg) * self.pm

        return r

    def update_am(self, r):
        self.am += (r[0] - r[1] + r[2] - r[3]) * self.deltat

    def update_pm(self, r):
        self.pm += (r[5] - r[6] + r[7] - r[8]) * self.deltat

    def update_ac(self, r):
        self.ac += (- r[0] + r[1] - r[2] + r[3]) * self.psi * self.deltat

    def update_pc(self, r):
        self.pc += (- r[5] + r[6] - r[7] + r[8]) * self.psi * self.deltat

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
        np.savetxt(direc + 'ac.txt', self.ac)
        np.savetxt(direc + 'am.txt', self.am)
        np.savetxt(direc + 'pc.txt', self.pc)
        np.savetxt(direc + 'pm.txt', self.pm)
        np.savetxt(direc + 'time.txt', [self.time])
