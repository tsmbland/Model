import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import copy
from matplotlib import animation
from Models.WavePinning import Model
import itertools
from scipy import interpolate

"""
Nullcline

"""


def evaluate(func, xrange, yrange, iterations, resolution=100, args=()):
    for iteration in range(iterations):

        if iteration == 0:

            # Set x and y values
            n_sims = resolution
            xvals = np.linspace(xrange[0], xrange[1], n_sims)
            yvals = np.linspace(yrange[0], yrange[1], n_sims)

            # Evaluate
            res = func(np.tile(xvals, (n_sims, 1)).T, np.tile(yvals, (n_sims, 1)), *args)
            res_2d_sign = np.sign(res)

        else:
            # Find boundary regions
            a = np.nonzero(np.nan_to_num(np.diff(res_2d_sign, axis=0)))
            b = np.nonzero(np.nan_to_num(np.diff(res_2d_sign, axis=1)))
            c = np.nonzero(np.nan_to_num(res_2d_sign[:-1, :-1] - res_2d_sign[1:, 1:]))
            xpoints = np.r_[a[0], b[0], c[0]]
            ypoints = np.r_[a[1], b[1], c[1]]

            # Set x and y values
            n_sims = resolution * (n_sims - 1) + 1
            run_bool = np.zeros([n_sims, n_sims])
            for x, y in zip(xpoints, ypoints):
                run_bool[x * resolution:x * resolution + (resolution + 1),
                y * resolution:y * resolution + (resolution + 1)] = 1
            sims_array_ind = np.nonzero(run_bool)
            xvals = xrange[0] + sims_array_ind[0] * (xrange[1] - xrange[0]) / n_sims
            yvals = yrange[0] + sims_array_ind[1] * (yrange[1] - yrange[0]) / n_sims

            # Evaluate
            res = func(xvals, yvals, *args)

            # Organise res
            res_2d = np.nan * np.zeros([n_sims, n_sims])
            for r in range(len(res)):
                res_2d[sims_array_ind[0][r], sims_array_ind[1][r]] = res[r]
            res_2d_sign = np.sign(res_2d)

    a = np.nonzero(np.nan_to_num(np.diff(res_2d_sign, axis=0)))
    b = np.nonzero(np.nan_to_num(np.diff(res_2d_sign, axis=1)))
    c = np.nonzero(np.nan_to_num(res_2d_sign[:-1, :-1] - res_2d_sign[1:, 1:]))
    xpoints = np.r_[a[0], b[0], c[0]]
    ypoints = np.r_[a[1], b[1], c[1]]

    xpoints = xrange[0] + (xpoints / n_sims) * (xrange[1] - xrange[0])
    ypoints = yrange[0] + (ypoints / n_sims) * (yrange[1] - yrange[0])
    return xpoints, ypoints


def func(b, a, k0=0.067, lamda=1, K=1, delta=1):
    return b * (k0 + (lamda * (a ** 2)) / (K ** 2 + a ** 2)) - delta * a


xpoints, ypoints = evaluate(func, (0, 2.5), (0, 2.5), resolution=50, iterations=2, args=())
sign = np.sign(func(xpoints, ypoints + 0.01))
plt.scatter(xpoints, ypoints, s=0.3, c=sign)
plt.xlabel('Cytoplasmic concentration (μm⁻³)')
plt.ylabel('Membrane concentration (μm⁻²)')
sns.despine()
fig = plt.gcf()
fig.set_size_inches(3.5, 3.5)
plt.tight_layout()
plt.show()
# plt.savefig('nullcline.png', dpi=600)

"""
Animation

"""

# base_model = Model(Dm=0.1, Dc=2, k0=0.067, lamda=1, K=1, delta=1, xsteps=100, Tmax=150, deltat=0.001, deltax=0.1,
#                    m_0=np.zeros([100]), c_0=2.3 * np.ones([100]))
#
#
# def animate(model, p1, framerate, name, ymax=None, y2max=None):
#     """
#     p1: model time between frames
#
#     """
#
#     # Set up model
#     model = copy.deepcopy(base_model)
#
#     # Initial equilibration (no antagonism)
#     for t in range(10000):
#         model.update()
#
#     # Polarise
#     model.m *= np.r_[np.ones([95]), 4 * np.ones([5])]
#     for t in range(1000):
#         model.update()
#
#     # Seu up figure
#     plt.clf()
#     fig = plt.figure()
#     fig.set_size_inches(4.5, 3.5)
#     ax2 = fig.add_subplot(111)
#     ax1 = ax2.twinx()
#     frames = range(0, int(model.Tmax / p1))
#
#     def update_anim(i):
#         ax1.clear()
#         ax2.clear()
#
#         # Plot current timepoint
#         ax1.plot(model.c, c=(0.5, 0.5, 0.5))
#         ax2.plot(model.m, c=(0, 0, 0))
#
#         # Configure plot
#         ax1.set_title('Time (s): {0:.0f}'.format(model.time))
#         ax1.set_xticks([], [])
#         ax1.set_xlabel('Anterior                                                  Posterior')
#
#         ax1.set_ylabel('Cytoplasmic concentration (μm⁻³)', c=(0.5, 0.5, 0.5))
#         ax1.set_ylim(bottom=0, top=y2max)
#         # ax1.spines['left'].set_color('c')
#         ax1.tick_params(axis='y', colors=(0.5, 0.5, 0.5))
#
#         ax2.set_ylabel('Membrane concentration (μm⁻²)', c=(0, 0, 0))
#         ax2.set_ylim(bottom=0, top=ymax)
#         # ax2.spines['right'].set_color('m')
#         # ax2.tick_params(axis='y', colors='m')
#         plt.tight_layout()
#
#         # Run model until next save point
#         if i != 0:
#             for t in range(int(p1 / model.deltat)):
#                 model.update()
#                 model.time += model.deltat
#
#     # Animate
#     paranim = animation.FuncAnimation(fig, update_anim, frames=iter(frames), save_count=len(frames))
#     writer = animation.writers['ffmpeg']
#     writer = writer(fps=framerate, bitrate=2000)
#     paranim.save(name, writer=writer, dpi=300)
#
#
# animate(base_model, p1=1, framerate=24, name='anim.mp4', ymax=1.5, y2max=4)
