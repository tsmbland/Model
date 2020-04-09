import matplotlib
import copy
import multiprocessing
import itertools
import random
import sys
import os
import numpy as np

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')

from Models.SimplePosFeedback import Model
import glob
import matplotlib.pyplot as plt
from matplotlib import animation

def scatter_fig(direc):
    fig, ax = plt.subplots()
    res = np.zeros([len(glob.glob(direc + '/*.txt'))])
    kAP = np.zeros([len(glob.glob(direc + '/*.txt'))])
    kPA = np.zeros([len(glob.glob(direc + '/*.txt'))])
    for i, a in enumerate(glob.glob(direc + '/*.txt')):
        res[i] = np.loadtxt(a)
        name = a.split('/')[-1][:-4]
        kAP[i] = float(name.split('_')[0])
        kPA[i] = float(name.split('_')[1])
    ax.scatter(kAP, kPA, c=res / max(res), s=2, cmap='gray')

    best = np.argmin(res)
    print(res[best])
    print(kAP[best])
    print(kPA[best])


scatter_fig(home_direc + '/Res')
plt.show()

# """
# Animate
#
# """
#
#
# def animate(model, p1, framerate, name, ymax=None):
#     """
#
#     p1: model time between frames
#
#     """
#
#     # Seu up figure
#     plt.clf()
#     fig = plt.figure()
#     fig.set_size_inches(4.5, 3.5)
#     ax = fig.add_subplot(111)
#     frames = range(0, int(model.Tmax / p1))
#
#     def update_anim(i):
#         ax.clear()
#
#         # Plot current timepoint
#         ax.plot(model.am / 2.8, c='r')
#         ax.plot(model.pm / 5, c='b')
#
#         # Configure plot
#         ax.set_title('Time (s): {0:.0f}'.format(model.time))
#         ax.set_ylim(bottom=0, top=ymax)
#         ax.set_xticks([], [])
#         ax.set_xlabel('Anterior                                                  Posterior')
#         ax.set_ylabel('Membrane Concentration (a.u.)')
#
#         # Run model until next save point
#         if i != 0:
#             for t in range(int(p1 / model.deltat)):
#                 model.react()
#                 model.time += model.deltat
#
#     # Animate
#     paranim = animation.FuncAnimation(fig, update_anim, frames=iter(frames), save_count=len(frames))
#     writer = animation.writers['ffmpeg']
#     writer = writer(fps=framerate, bitrate=2000)
#     paranim.save(name, writer=writer, dpi=300)
#
#
# """
# Import data
#
# """
#
# kAP_vals = np.loadtxt(home_direc + '/kAP_vals.txt')
# kPA_vals = np.loadtxt(home_direc + '/kPA_vals.txt')
# res = np.loadtxt(home_direc + '/Res.txt')
#
# a_mem = np.loadtxt(home_direc + '/a_mem.txt')
# a_cyt = 0.17147142857142855
# p_mem = np.loadtxt(home_direc + '/p_mem.txt')
# p_cyt = 0.06860000000000001
#
# """
# Run
#
# """
#
# koffA = 0.0092
# koffP = 0.0073
# psi = 0.10318684114244771
# dosP = 0.294005475663175
# dosA = 1.05143336288 * dosP
#
# base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
#                    koffP=koffP,
#                    kposP=12.7357711156 * koffP, kAP=10 ** -1.203125, kPA=10 ** -1.1484375, ePneg=1, eAneg=1, xsteps=100,
#                    Tmax=1000, deltat=0.01, deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]),
#                    pc_0=dosP)
#
# animate(base_model, p1=3, framerate=24, name='anim.mp4')
