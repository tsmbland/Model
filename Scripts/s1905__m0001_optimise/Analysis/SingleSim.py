import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '../')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '../../')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '../../..')
sys.path.append('../..')
import InVivo.s1905__par2_ng as a
import M as x
import Models.m0001 as m
import numpy as np
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import shutil
import copy

"""


"""

params = {'Da': 0, 'Dp': 0.15, 'konA': 0.0, 'koffA': 0, 'konP': 0.09489489, 'koffP': 0.0073,
          'kposP': 0.0968969,
          'kAP': 0.0, 'kPA': 0.01, 'ePneg': 1, 'eAneg': 1, 'xsteps': 500, 'psi': a.svr, 'Tmax': 10000,
          'deltat': 0.01, 'deltal': a.spres, 'radii': a.r_local, 'am_0': a.p6_mem * 0 + 0, 'ac_0': a.p6_cyt * 0,
          'pm_0': a.p2_mem, 'pc_0': a.p2_cyt}

sdirec = '_temp'
deltas = 10

"""
Run

"""

if os.path.exists(sdirec):
    shutil.rmtree(sdirec)
os.mkdir(sdirec)

model = m.Model(**dict(copy.deepcopy(params)))
model.Tmax = 100

if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
    os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))

for t in range(int(model.Tmax / model.deltat)):
    model.react()
    model.time = (t + 1) * model.deltat

    if model.time % deltas == 0:
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
    am = np.loadtxt(direc + '/am.txt')
    pm = np.loadtxt(direc + '/pm.txt')
    ax.plot(am, c='k')
    ax.plot(pm, c='r')
    ax.plot(params['am_0'], c='k', linestyle='--')
    ax.plot(params['pm_0'], c='r', linestyle='--')
    ax.set_ylim(bottom=-0.5, top=3)
    ax.set_ylabel('Cortical concentration')

    # Cytoplasmic
    ac = np.loadtxt(direc + '/ac.txt')
    pc = np.loadtxt(direc + '/pc.txt')
    ax2.bar(1, ac, color='k', alpha=0.2)
    ax2.bar(2, pc, color='r', alpha=0.2)
    # ax2.bar(1, params['ac_0'], color='k', alpha=0.2, linewidth=5, linestyle='--')
    # ax2.bar(2, params['pc_0'], color='r', alpha=0.2, linewidth=5, linestyle='--')
    ax2.set_xticks([])
    ax2.set_ylabel('Cytoplasmic concentration')
    ax2.set_ylim(bottom=0, top=0.2)

    sns.despine()


sframe.on_changed(update)
plt.show()
