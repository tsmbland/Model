import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '../')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '../../..')
import M as x
import Models.m0008 as m
import numpy as np
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import shutil
import copy

pdirec = x.ddirec + 's1904__m0008_optimise/g007/233'
params = pickle.load(open('%s/_params.pkl' % pdirec, 'rb'))
sdirec = '_temp'
deltas = 10

if os.path.exists(sdirec):
    shutil.rmtree(sdirec)
os.mkdir(sdirec)

"""
Run

"""

model = m.Model(**dict(copy.deepcopy(params)))
model.Tmax = 1000

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
    pm = np.loadtxt(direc + '/pm1.txt') + np.loadtxt(direc + '/pm2d.txt') + np.loadtxt(direc + '/pm2s.txt')
    pm_0 = params['pm1_0'] + params['pm2d_0'] + params['pm2s_0']
    ax.plot(am, c='k')
    ax.plot(pm, c='r')
    ax.plot(params['am_0'], c='k', linestyle='--')
    ax.plot(pm_0, c='r', linestyle='--')
    ax.set_ylim(bottom=-0.5, top=5)
    ax.set_ylabel('Cortical concentration')

    # Cytoplasmic
    ac = np.loadtxt(direc + '/ac.txt')
    pc = np.loadtxt(direc + '/pc1.txt') + np.loadtxt(direc + '/pc2.txt')
    pc_0 = params['pc1_0'] + params['pc2_0']
    ax2.bar(1, ac, color='k', alpha=0.2)
    ax2.bar(2, pc, color='r', alpha=0.2)
    ax2.bar(1, params['ac_0'], color='k', alpha=0.2, linewidth=5, linestyle='--')
    ax2.bar(2, pc_0, color='r', alpha=0.2, linewidth=5, linestyle='--')
    ax2.set_xticks([])
    ax2.set_ylabel('Cytoplasmic concentration')
    ax2.set_ylim(bottom=0, top=0.2)

    sns.despine()


sframe.on_changed(update)
plt.show()
