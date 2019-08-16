import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '../')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '../../..')
import M as x
import Models.m0000 as m
import numpy as np
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import shutil
import copy
from matplotlib import animation

pdirec = x.ddirec + 's1904__m0000_optimise_B/g081/339'
# pdirec = x.ddirec + 's1904__m0000_optimise_B/g053/902'
# pdirec = x.ddirec + 's1904__m0000_optimise_B/g034/634'
params = pickle.load(open('%s/_params.pkl' % pdirec, 'rb'))

print('konA = ' + str(params['konA']))
print('konP = ' + str(params['konP']))
print('kAP = ' + str(params['kAP']))
print('kPA = ' + str(params['kPA']))

sdirec = '_temp'
deltas = 10

# """
# Run
#
# """
#
# if os.path.exists(sdirec):
#     shutil.rmtree(sdirec)
# os.mkdir(sdirec)
#
# model = m.Model(**dict(copy.deepcopy(params)))
# model.Tmax = 1000
#
# if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
#     os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
# model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))
#
# for t in range(int(model.Tmax / model.deltat)):
#     model.react()
#     model.time = (t + 1) * model.deltat
#
#     if model.time % deltas == 0:
#         if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
#             os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
#         model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))

# """
# Run - CRT90
#
# """
#
# sdirec = '_temp2'
#
# if os.path.exists(sdirec):
#     shutil.rmtree(sdirec)
# os.mkdir(sdirec)
#
# model = m.Model(**dict(copy.deepcopy(params)))
# model.Tmax = 1000
# model.kPA = 0
#
# if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
#     os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
# model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))
#
# for t in range(int(model.Tmax / model.deltat)):
#     model.react()
#     model.time = (t + 1) * model.deltat
#
#     if model.time % deltas == 0:
#         if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
#             os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
#         model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))
#
# """
# Run - PAR-1 RNAi
#
# """
#
# sdirec = '_temp3'
#
# if os.path.exists(sdirec):
#     shutil.rmtree(sdirec)
# os.mkdir(sdirec)
#
# model = m.Model(**dict(copy.deepcopy(params)))
# model.Tmax = 1000
# model.kAP = 0
#
# if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
#     os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
# model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))
#
# for t in range(int(model.Tmax / model.deltat)):
#     model.react()
#     model.time = (t + 1) * model.deltat
#
#     if model.time % deltas == 0:
#         if not os.path.exists('%s/%s' % (sdirec, '{:05d}'.format(int(model.time)))):
#             os.mkdir('%s/%s' % (sdirec, '{:05d}'.format(int(model.time))))
#         model.save('%s/%s/' % (sdirec, '{:05d}'.format(int(model.time))))

# """
# Slider plot
#
# """
#
# from matplotlib.widgets import Slider
#
# sdirec = '_temp'
#
# direcs = x.direcslist(sdirec)
#
# ax = plt.subplot2grid((1, 2), (0, 0))
# ax2 = plt.subplot2grid((1, 2), (0, 1))
# plt.subplots_adjust(bottom=0.25, wspace=0.5)
# axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
# sframe = Slider(axframe, 'Iteration', 0, len(direcs), valinit=0, valfmt='%d')
#
#
# def update(i):
#     ax.clear()
#     ax2.clear()
#     direc = direcs[int(i)]
#     print(direc)
#
#     # Cortical
#     am = np.loadtxt(direc + '/am.txt')
#     pm = np.loadtxt(direc + '/pm.txt')
#     ax.plot(am, c='k')
#     ax.plot(pm, c='r')
#     ax.plot(params['am_0'], c='k', linestyle='--')
#     ax.plot(params['pm_0'], c='r', linestyle='--')
#     ax.set_ylim(bottom=-0.5, top=5)
#     ax.set_ylabel('Cortical concentration')
#
#     # Cytoplasmic
#     ac = np.loadtxt(direc + '/ac.txt')
#     pc = np.loadtxt(direc + '/pc.txt')
#     ax2.bar(1, ac, color='k', alpha=0.2)
#     ax2.bar(2, pc, color='r', alpha=0.2)
#     ax2.bar(1, params['ac_0'], color='k', alpha=0.2, linewidth=5, linestyle='--')
#     ax2.bar(2, params['pc_0'], color='r', alpha=0.2, linewidth=5, linestyle='--')
#     ax2.set_xticks([])
#     ax2.set_ylabel('Cytoplasmic concentration')
#     ax2.set_ylim(bottom=0, top=0.2)
#
#     sns.despine()
#
#
# sframe.on_changed(update)
# plt.show()

# """
# Animation
#
# """
#
# sdirec = '_temp'
#
# direcs = x.direcslist(sdirec)
# plt.clf()
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
#
# def update(i):
#     ax.clear()
#     direc = direcs[int(i)]
#     print(direc)
#
#     # Cortical
#     am = np.loadtxt(direc + '/am.txt')
#     pm = np.loadtxt(direc + '/pm.txt')
#     time = np.loadtxt(direc + '/time.txt')
#     ax.plot(am, c='r')
#     ax.plot(pm, c='b')
#     ax.plot(params['am_0'], c='r', linestyle='--')
#     ax.plot(params['pm_0'], c='b', linestyle='--')
#     ax.set_ylim(bottom=-0.5, top=5)
#     ax.set_ylabel('Cortical concentration')
#     ax.set_xticks([])
#     ax.set_xlabel('Anterior - Posterior')
#     ax.set_title('Time = ' + str(int(time)) + ' seconds')
#
#     sns.despine()
#
#
# paranim = animation.FuncAnimation(fig, update, frames=np.r_[0, np.arange(len(direcs))])
# writer = animation.writers['ffmpeg']
# writer = writer(fps=20, bitrate=2000)
# paranim.save('animation1.mp4', writer=writer)

"""
Final state

"""

sdirec = '_temp'
direc = x.direcslist(sdirec)[-1]

fig = plt.figure()
ax = fig.add_subplot(111)
am = np.loadtxt(direc + '/am.txt')
pm = np.loadtxt(direc + '/pm.txt')
time = np.loadtxt(direc + '/time.txt')
ax.plot(am, c='r')
ax.plot(pm, c='b')
ax.plot(params['am_0'], c='r', linestyle='--')
ax.plot(params['pm_0'], c='b', linestyle='--')
ax.set_ylim(bottom=-0.5, top=5)
ax.set_ylabel('Cortical concentration')
ax.set_xticks([])
ax.set_xlabel('Anterior - Posterior')
sns.despine()
plt.show()



# """
# ASI over time
#
# """
#
# sdirec = '_temp2'
#
# times = np.zeros([len(x.direcslist(sdirec))])
# asis = np.zeros([len(x.direcslist(sdirec))])
# for i, d in enumerate(x.direcslist(sdirec)):
#     times[i] = np.loadtxt(d + '/time.txt')
#     asis[i] = abs(x.asi(np.loadtxt(d + '/pm.txt')))
# plt.plot(times, asis)
# plt.xlabel('Time')
# plt.ylabel('pPAR ASI')
# sns.despine()
# plt.show()
