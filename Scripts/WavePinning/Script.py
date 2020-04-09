from Models.WavePinning import Model
import numpy as np
import os
import shutil
import copy
import matplotlib.pyplot as plt
from M import direcslist

base_model = Model(Dm=0.1, Dc=2, k0=0.067, lamda=1, K=1, delta=1, xsteps=100, Tmax=500, deltat=0.001, deltax=0.1,
                   m_0=np.zeros([100]), c_0=2.3 * np.ones([100]))

deltas = 1
sdirec = '_temp'

"""
Run

"""

if os.path.isdir(sdirec):
    shutil.rmtree(sdirec)
os.mkdir(sdirec)

# Set up model
model = copy.deepcopy(base_model)

# Initial equilibration (no antagonism)
for t in range(10000):
    model.update()

# Polarise
model.m *= np.r_[np.ones([95]), 4 * np.ones([5])]

# Simulation + saving
for t in range(int(model.Tmax / model.deltat)):
    if model.time % deltas == 0 or model.time == 0:
        if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
            os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
        model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))

    model.update()
    model.time = (t + 1) * model.deltat

"""
View

"""

from matplotlib.widgets import Slider

direcs = sorted(direcslist(sdirec))
fig = plt.figure()
ax = fig.add_subplot(111)
plt.subplots_adjust(bottom=0.25, wspace=0.5)
axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
sframe = Slider(axframe, 'Iteration', 0, len(direcs), valinit=0, valfmt='%d')


def update(i):
    ax.clear()
    direc = direcs[int(i)]
    print(direc)

    m = np.loadtxt(direc + '/m.txt')
    c = np.loadtxt(direc + '/c.txt')
    ax.plot(m, c='k')
    ax.plot(c, c='r')
    ax.set_ylim(bottom=0)
    ax.set_ylabel('Concentration')


sframe.on_changed(update)
plt.show()
