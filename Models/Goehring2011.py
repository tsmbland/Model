import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.widgets import Slider
import glob
import shutil

"""


"""


class Model:
    def __init__(self, Da, Dp, konA, koffA, konP, koffP, kPA, kAP, alpha, beta, xsteps, psi, Tmax, deltat, L,
                 pA, pP):

        # Species
        self.am = np.zeros([int(xsteps)])
        self.ac = pA
        self.pm = np.zeros([int(xsteps)])
        self.pc = pP
        self.time = 0

        # Diffusion
        self.Da = Da  # input is um2 s-1
        self.Dp = Dp  # um2 s-1

        # Membrane exchange
        self.konA = konA  # um s-1
        self.koffA = koffA  # s-1
        self.konP = konP  # um s-1
        self.koffP = koffP  # s-1

        # Antagonism
        self.kPA = kPA  # um4 s-1
        self.kAP = kAP  # um2 s-1
        self.alpha = alpha
        self.beta = beta

        # Misc
        self.L = L
        self.xsteps = int(xsteps)
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s
        self.deltax = self.L / xsteps  # um
        self.psi = psi  # um-1

    def diffusion(self, concs):
        d = concs[np.append(np.array(range(1, len(concs))), [len(concs) - 1])] - 2 * concs + concs[
            np.append([0], np.array(range(len(concs) - 1)))]
        return d / (self.deltax ** 2)

    def reactions(self):
        """
        r0: a on
        r1: a off
        r2: p to a antagonism
        r3: a diffusion

        r4: p on
        r5: p off
        r6: a to p antagonism
        r7: p diffusion

        """

        r = [None] * 8

        r[0] = self.konA * self.ac
        r[1] = self.koffA * self.am
        r[2] = self.kAP * (self.pm ** self.beta) * self.am
        r[3] = self.Da * self.diffusion(self.am)

        r[4] = self.konP * self.pc
        r[5] = self.koffP * self.pm
        r[6] = self.kPA * (self.am ** self.alpha) * self.pm
        r[7] = self.Dp * self.diffusion(self.pm)

        return r

    def update_am(self, r):
        self.am += (r[0] - r[1] - r[2] + r[3]) * self.deltat

    def update_pm(self, r):
        self.pm += (r[4] - r[5] - r[6] + r[7]) * self.deltat

    def update_ac(self, r):
        self.ac += (- r[0] + np.mean(r[1]) + np.mean(r[2])) * self.psi * self.deltat

    def update_pc(self, r):
        self.pc += (- r[4] + np.mean(r[5]) + np.mean(r[6])) * self.psi * self.deltat

    def react(self):
        r = self.reactions()
        self.update_am(r)
        self.update_ac(r)
        self.update_pm(r)
        self.update_pc(r)

    def initiate(self):
        # Remove antagonism
        kAP, kPA = self.kAP, self.kPA
        self.kAP = 0
        self.kPA = 0

        # Initial equilibration (no antagonism)
        for t in range(10000):
            self.react()

        # Polarise
        self.am *= 2 * np.r_[np.ones([self.xsteps // 2]), np.zeros([self.xsteps // 2])]
        self.pm *= 2 * np.r_[np.zeros([self.xsteps // 2]), np.ones([self.xsteps // 2])]

        # Add back antagonism
        self.kAP = kAP
        self.kPA = kPA

    def run(self, save_direc=None, save_gap=10):
        """
        save_direc: directory to save model data (must already exist)
        save_freq: frequency to save (units of model time)

        """

        # Set up save directory
        if save_direc is not None:
            [shutil.rmtree(f) for f in glob.glob('%s/*/' % save_direc)]

        # Run
        for t in range(int(self.Tmax / self.deltat)):
            if save_direc is not None:
                if self.time % save_gap == 0 or self.time == 0:
                    if not os.path.exists('%s/%s' % (save_direc, '{:05d}'.format(int(self.time)))):
                        os.mkdir('%s/%s' % (save_direc, '{:05d}'.format(int(self.time))))
                    self.save('%s/%s/' % (save_direc, '{:05d}'.format(int(self.time))))

            self.react()
            self.time = (t + 1) * self.deltat

    def save(self, direc):
        np.savetxt(direc + '/ac.txt', [self.ac])
        np.savetxt(direc + '/am.txt', self.am)
        np.savetxt(direc + '/pc.txt', [self.pc])
        np.savetxt(direc + '/pm.txt', self.pm)
        np.savetxt(direc + '/time.txt', [self.time])

    def plot_res(self):
        plt.plot(np.linspace(0, self.L, self.xsteps), self.am, c='tab:red')
        plt.plot(np.linspace(0, self.L, self.xsteps), self.pm, c='tab:blue')
        plt.ylim(bottom=0)
        plt.ylabel('Cortical concentration (a.u.)')
        plt.xlabel('Position (μm)')
        plt.show()

    def animate(self, direc):
        """
        direc: directory to results (as specified in run function)

        """

        direcs = sorted(glob.glob('%s/*/' % direc))
        times = np.array([int(os.path.basename(os.path.normpath(direcs[i]))) for i in range(len(direcs))])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.25, wspace=0.5)
        axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
        sframe = Slider(axframe, 'Time (s)', 0, times[-1], valinit=0, valfmt='%d')

        def update(i):
            ax.clear()
            direc = direcs[np.argmin(abs(times - int(i)))]
            am = np.loadtxt(direc + '/am.txt')
            pm = np.loadtxt(direc + '/pm.txt')
            ax.plot(np.linspace(0, self.L, self.xsteps), am, c='tab:red')
            ax.plot(np.linspace(0, self.L, self.xsteps), pm, c='tab:blue')
            ax.set_ylim(bottom=0)
            ax.set_ylabel('Cortical concentration (a.u.)')
            ax.set_xlabel('Position (μm)')

        sframe.on_changed(update)
        plt.show()
