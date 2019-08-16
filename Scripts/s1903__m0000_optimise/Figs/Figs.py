import matplotlib.pyplot as plt
import numpy as np
import M as x
import pickle

direc = '/Users/blandt/Documents/PhDWork/Code/ModelData/22222222'

"""
Plot distribution

"""

# for d in x.direcslist(direc):
#     p = pickle.load(open('%s/_params.pkl' % d, 'rb'))
#
#     am = np.loadtxt(d + '/am.txt')
#     pm = np.loadtxt(d + '/pm.txt')
#
#     plt.plot(p['am_0'], c='k', alpha=0.2)
#     plt.plot(p['pm_0'], c='r', alpha=0.2)
#
#     plt.plot(am, c='k')
#     plt.plot(pm, c='r')
#
#     plt.title(np.mean([np.mean((p['am_0'] - am) ** 2), np.mean((p['pm_0'] - pm) ** 2)]))
#     plt.show()

"""
Ordered slider

"""

from matplotlib.widgets import Slider

direcs = x.direcslist(direc)
mses = np.zeros([len(direcs)])

for i, d in enumerate(direcs):
    p = pickle.load(open('%s/_params.pkl' % d, 'rb'))
    am = np.loadtxt(d + '/am.txt')
    pm = np.loadtxt(d + '/pm.txt')
    mses[i] = np.mean([np.mean((p['am_0'] - am) ** 2), np.mean((p['pm_0'] - pm) ** 2)])

direcs = [direcs[i] for i in np.argsort(mses)]

fig = plt.figure()
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.25, bottom=0.25)
axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
sframe = Slider(axframe, 'Iteration', 0, len(direcs), valinit=0, valfmt='%d')


def update(i):
    ax.clear()
    direc = direcs[int(i)]
    print(direc)

    p = pickle.load(open('%s/_params.pkl' % direc, 'rb'))
    am = np.loadtxt(direc + '/am.txt')
    pm = np.loadtxt(direc + '/pm.txt')
    ax.plot(p['am_0'], c='k', alpha=0.2)
    ax.plot(p['pm_0'], c='r', alpha=0.2)
    ax.plot(am, c='k')
    ax.plot(pm, c='r')
    ax.set_title(np.mean([np.mean((p['am_0'] - am) ** 2), np.mean((p['pm_0'] - pm) ** 2)]))


sframe.on_changed(update)
plt.show()
