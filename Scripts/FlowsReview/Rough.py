import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
import M as x
import Models.FlowsReview as m
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sdirec = '_temp'
deltas = 1

"""
Run

"""

model = m.Model(Dm=1, Dc=10, Vm=0.01, Vc=0, kon=0.002, koff=0.001, xsteps=100, Tmax=100, deltat=0.001, deltax=0.1,
                c_0=1, m_0=1, flowtype=3)
for t in range(int(model.Tmax / model.deltat)):
    model.react()
    model.time = (t + 1) * model.deltat

    if model.time % deltas == 0 or model.time == model.deltat:
        if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
            os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
        model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))

"""
Slider plot

"""

from matplotlib.widgets import Slider

direcs = x.direcslist(sdirec)

ax = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))
plt.subplots_adjust(bottom=0.25, wspace=0.5)
axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
sframe = Slider(axframe, 'Iteration', 0, len(direcs), valinit=0, valfmt='%d')


def update(i):
    ax.clear()
    ax2.clear()
    direc = direcs[int(i)]
    print(direc)

    # Cortical
    m = np.loadtxt(direc + '/m.txt')
    ax.plot(m, c='k')
    ax.set_ylim(bottom=0)
    ax.set_ylabel('Cortical concentration')

    # Cytoplasmic
    c = np.loadtxt(direc + '/c.txt')
    ax2.plot(c, c='k')
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel('Cytoplasmic concentration')

    sns.despine()


sframe.on_changed(update)
plt.show()
