from Models.SimplePosFeedback import Model
import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
import M as x
import copy

base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
                   koffP=0.01, kposP=0.02, kAP=10 ** -1.5, kPA=10 ** -0.75, ePneg=1, eAneg=2, xsteps=100, Tmax=1000,
                   deltat=0.01,
                   deltax=0.5, psi=0.1, am_0=0 * np.ones([100]), ac_0=1, pm_0=10 * np.ones([100]), pc_0=0)

"""
Run

"""

deltas = 1
sdirec = '_temp'


def run():
    if os.path.isdir(sdirec):
        shutil.rmtree(sdirec)
    os.mkdir(sdirec)

    # Set up model
    model = copy.deepcopy(base_model)

    # Remove antagonism
    # model.kAP = 0
    # model.kPA = 0

    # Initial equilibration (no antagonism)
    # for t in range(10000):
    #     model.react()

    # Polarise
    # model.am *= 2 * np.r_[np.ones([50]), np.zeros([50])]
    # model.pm *= 2 * np.r_[np.zeros([50]), np.ones([50])]

    # Add antagonism
    # model.kAP = 10 ** -1.5
    # model.kPA = 10 ** -0.75

    # Simulation + saving
    for t in range(int(model.Tmax / model.deltat)):
        if model.time % deltas == 0 or model.time == 0:
            if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
                os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
            model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))

        model.react()
        model.time = (t + 1) * model.deltat


run()

"""
Slider plot

"""

from matplotlib.widgets import Slider

direcs = sorted(x.direcslist(sdirec))
fig = plt.figure()
ax = fig.add_subplot(111)
plt.subplots_adjust(bottom=0.25, wspace=0.5)
axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
sframe = Slider(axframe, 'Iteration', 0, len(direcs), valinit=0, valfmt='%d')


def update(i):
    ax.clear()
    direc = direcs[int(i)]
    print(direc)

    am = np.loadtxt(direc + '/am.txt')
    pm = np.loadtxt(direc + '/pm.txt')
    ax.plot(am, c='k')
    ax.plot(pm, c='r')
    ax.set_ylim(bottom=0)
    ax.set_ylabel('Cortical concentration')


sframe.on_changed(update)
plt.show()
