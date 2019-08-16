import sys
import os
import itertools
from joblib import Parallel, delayed
import multiprocessing

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../../../MembraneQuant')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../../../ImageAnalysis_Scripts')

import numpy as np
import IA2 as ia2
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import InVivo.s1905__par2_ng as iv
from matplotlib.widgets import Slider

#############################################################################################


"""
Functions

"""


def evaluate(func, xrange, yrange, resolution, args):
    res = func(np.tile(np.linspace(xrange[0], xrange[1], resolution), (resolution, 1)).T,
               np.tile(np.linspace(yrange[0], yrange[1], resolution), (resolution, 1)), *args)
    res_sign = np.sign(res)
    a = np.nonzero(np.diff(res_sign, axis=0))
    b = np.nonzero(np.diff(res_sign, axis=1))
    xpoints = xrange[0] + (np.r_[a[0], b[0]] / resolution) * (xrange[1] - xrange[0])
    ypoints = yrange[0] + (np.r_[a[1], b[1]] / resolution) * (yrange[1] - yrange[0])
    return xpoints, ypoints


def mse(xdata, ydata, func, yrange, args):
    resolution = 1000
    errors = np.zeros([len(xdata)])
    for i, x in enumerate(xdata):
        a = np.nonzero(np.diff(np.sign(func(x, np.linspace(yrange[0], yrange[1], resolution), *args))))[0]
        ypoints = yrange[0] + (a / resolution) * (yrange[1] - yrange[0])
        if len(ypoints != 0):
            errors[i] = min(ydata[i] - ypoints) ** 2
        else:
            return np.nan

    return np.mean(errors)


#############################################################################################

"""
Import data

"""

# d = ia2.ddirec + '/Experiments/e1905__par2_ng_rundown_tom9/'
#
# direcs = ia2.direcslist(d + 'lp637_wt', 1) + ia2.direcslist(d + 'nwg201_wt', 1) + ia2.direcslist(d + 'nwg201_rd',
#                                                                                                  1) + ia2.direcslist(
#     d + 'lp637_rd', 1)
# mem = 0.001 * 0.255 * 1.63672480374 * ia2.bounded_mean_2d(
#     np.array(ia2.importall(direcs, 'Quantification_g/mem_spa1000.txt')),
#     [0.9, 0.1])
# cyt = 0.001 * (0.255 ** 2) * ia2.bounded_mean_2d(np.array(ia2.importall(direcs, 'Quantification_g/cyt_spa1000.txt')),
#                                                  [0.9, 0.1])

# plt.scatter(cyt, mem)
# plt.show()

#############################################################################################

"""
Mass action

"""

# def func(x, y, kon, koff, kpos):
#     return kon * x - koff * y + kpos * x * y
#
#
# res = np.zeros([100, 100])
# kons = np.linspace(0, 0.6, 100)
# kposs = np.linspace(0, 0.15, 100)  # x axis
#
# for i, k1 in enumerate(kons):
#     print(i)
#     for j, k2 in enumerate(kposs):
#         res[i, j] = mse(cyt, mem, func, [-5, 20], (k1, 0.0073, k2))
#
# res[np.isnan(res)] = 1000
# a = np.where(res == np.min(res))
# print(np.min(res))
# print(kons[a[0]])
# print(kposs[a[1]])
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.imshow(res, cmap='gist_rainbow', extent=(0, 0.15, 0, 0.6), vmin=0, vmax=2, origin='lower')
# ax.set_aspect(0.25)
# plt.show()
#
# ###############
#
# plt.scatter(cyt, mem)
#
# xpoints, ypoints = evaluate(func, [0, 0.06], [0, 3], 1000, (0.0969697, 0.0073, 0.0969697))
# plt.scatter(xpoints, ypoints, s=0.1, c='k')
#
# xpoints, ypoints = evaluate(func, [0, 0.06], [0, 3], 1000, (0.0969697, 0.0073, 0))
# plt.scatter(xpoints, ypoints, s=0.1, c='g')
#
# plt.show()

#############################################################################################

"""
Hill

"""


def func(x, y, kon, koff, kpos, khill, exp):
    return kon * x - koff * y + kpos * x * ((y ** exp) / ((khill ** exp) + (y ** exp)))


# kon_max = 0.2
# kpos_max = 2
# khill_max = 20
#
# res = np.zeros([20, 20, 20])
# kons = np.linspace(0, kon_max, 20)
# kposs = np.linspace(0, kpos_max, 20)
# khills = np.linspace(0, khill_max, 20)
#
# for i, kon in enumerate(kons):
#     print(i)
#     for j, kpos in enumerate(kposs):
#         for z, khill in enumerate(khills):
#             res[i, j, z] = mse(cyt, mem, func, [-5, 20], (kon, 0.0073, kpos, khill, 1))
#
# res[np.isnan(res)] = 1000
# a = np.where(res == np.min(res))
# kon_opt = kons[a[0]]
# kpos_opt = kposs[a[1]]
# khill_opt = khills[a[2]]
# print(kon_opt)
# print(kpos_opt)
# print(khill_opt)
# print(res[a[0], a[1], a[2]])

##################

# plt.scatter(cyt, mem)
#
# xpoints, ypoints = evaluate(func, [0, 0.06], [0, 3], 1000, (0.05263158, 0.0073, 0.63157895, 3.15789474, 1))
# plt.scatter(xpoints, ypoints, s=0.1, c='k')
#
# xpoints, ypoints = evaluate(func, [0, 0.06], [0, 3], 1000, (0.05263158, 0.0073, 0, 0, 0))
# plt.scatter(xpoints, ypoints, s=0.1, c='g')
#
# plt.show()

##############################################################################################

"""
Import data

"""


def import_mem(direcs):
    return ia2.bounded_mean_2d(np.array(ia2.importall(direcs, 'Quantification_g/mem_spa1000.txt')), [0.9, 0.1])


def import_cyt(direcs):
    return ia2.bounded_mean_2d(np.array(ia2.importall(direcs, 'Quantification_g/cyt_spa1000.txt')), [0.9, 0.1])


d = ia2.ddirec + '/Experiments/e1905__par2_ng_rundown_tom9/'
a_mem = import_mem(ia2.direcslist(d, 2))
a_cyt = import_cyt(ia2.direcslist(d, 2))
a_mem_wt = import_mem(ia2.direcslist(d + 'lp637_wt', 1))
a_cyt_wt = import_cyt(ia2.direcslist(d + 'lp637_wt', 1))
a_mem_uniform = import_mem(ia2.direcslist(d + 'nwg201_wt', 1) + ia2.direcslist(d + 'nwg201_rd', 1))
a_cyt_uniform = import_cyt(ia2.direcslist(d + 'nwg201_wt', 1) + ia2.direcslist(d + 'nwg201_rd', 1))
a_mem_polar = import_mem(ia2.direcslist(d + 'lp637_wt', 1) + ia2.direcslist(d + 'lp637_rd', 1))
a_cyt_polar = import_cyt(ia2.direcslist(d + 'lp637_wt', 1) + ia2.direcslist(d + 'lp637_rd', 1))

d2 = ia2.ddirec + '/Experiments/e1901__par2_cai/'
b_mem = import_mem(ia2.direcslist(d2, 2))
b_cyt = import_cyt(ia2.direcslist(d2, 2))
b_mem_wt = import_mem(ia2.direcslist(d2 + 'kk1273_wt', 1))
b_cyt_wt = import_cyt(ia2.direcslist(d2 + 'kk1273_wt', 1))
b_mem_uniform = import_mem(ia2.direcslist(d2 + 'kk1273_par6', 1) + ia2.direcslist(d2 + 'th415_par6', 1))
b_cyt_uniform = import_cyt(ia2.direcslist(d2 + 'kk1273_par6', 1) + ia2.direcslist(d2 + 'th415_par6', 1))
b_mem_polar = import_mem(ia2.direcslist(d2 + 'kk1273_wt', 1) + ia2.direcslist(d2 + 'th415_wt', 1))
b_cyt_polar = import_cyt(ia2.direcslist(d2 + 'kk1273_wt', 1) + ia2.direcslist(d2 + 'th415_wt', 1))

d3 = ia2.ddirec + '/Experiments/e1905__ring_crispr/nwg214/'
ring_mem = import_mem(ia2.direcslist(d3, 1))
ring_cyt = import_cyt(ia2.direcslist(d3, 1))

b_mem *= (np.mean(a_mem_wt) / np.mean(b_mem_wt))
b_cyt *= (np.mean(a_cyt_wt) / np.mean(b_cyt_wt))
b_mem_uniform *= (np.mean(a_mem_wt) / np.mean(b_mem_wt))
b_cyt_uniform *= (np.mean(a_cyt_wt) / np.mean(b_cyt_wt))
b_mem_polar *= (np.mean(a_mem_wt) / np.mean(b_mem_wt))
b_cyt_polar *= (np.mean(a_cyt_wt) / np.mean(b_cyt_wt))
ring_cyt *= (np.mean(a_cyt_wt) / np.mean(b_cyt_wt))
ring_mem *= (np.mean(a_mem_wt) / np.mean(b_mem_wt))

cyt_polar = np.r_[a_cyt_polar, b_cyt_polar] * 0.001 * (0.255 ** 2)
mem_polar = np.r_[a_mem_polar, b_mem_polar] * 0.001 * 0.255 * 1.63672480374
cyt_uniform = np.r_[a_cyt_uniform, b_cyt_uniform] * 0.001 * (0.255 ** 2)
mem_uniform = np.r_[a_mem_uniform, b_mem_uniform] * 0.001 * 0.255 * 1.63672480374
ring_cyt *= 0.001 * (0.255 ** 2)
ring_mem *= 0.001 * 0.255 * 1.63672480374

plt.scatter(cyt_uniform, mem_uniform)
plt.scatter(cyt_polar, mem_polar)
plt.scatter(ring_cyt, ring_mem)
plt.show()

"""
Saturation model

"""


def func(x, y, kon, exp_on, khill_on, koff, kpos, khill_pos, exp_pos):
    return kon * (x ** exp_on) / ((khill_on ** exp_on) + (x ** exp_on)) - koff * y + kpos * x * (
        (y ** exp_pos) / ((khill_pos ** exp_pos) + (y ** exp_pos)))


# exp_pos = 1
# exp_on = 1
# koff = 0.0073
# kon_max = 5
# khill_on_max = 100
# kpos_max = 0.2
# khill_pos_max = 100
#
# kon_vals = np.linspace(0, kon_max, 10)
# kpos_vals = np.linspace(0, kpos_max, 10)
# khill_on_vals = np.linspace(0, khill_on_max, 10)
# khill_pos_vals = np.linspace(0, khill_pos_max, 10)
#
# res = np.zeros([10, 10, 10, 10])
# for i, kon in enumerate(kon_vals):
#     print(i)
#     for j, kpos in enumerate(kpos_vals):
#         for z, khill_on in enumerate(khill_on_vals):
#             for q, khill_pos in enumerate(khill_pos_vals):
#                 res[i, j, z, q] = mse(cyt_uniform, mem_uniform, func, [-5, 20],
#                                       (kon, 1, khill_on, 0.0073, kpos, khill_pos, 1))
#
# res[np.isnan(res)] = 1000
# a = np.where(res == np.min(res))
# print(np.min(res))
# print(kon_vals[a[0]])
# print(kpos_vals[a[1]])
# print(khill_on_vals[a[2]])
# print(khill_pos_vals[a[3]])

# plt.scatter(cyt_uniform, mem_uniform)
# xpoints, ypoints = evaluate(func, [0, 0.2], [0, 3], 1000,
#                             (1.11111111, 1, 11.11111111, 0.0073, 0.02222222, 1, 100))
# plt.scatter(xpoints, ypoints, s=0.1, c='k')
# plt.show()


# xrange = [0, 0.2]
# yrange = [0, 3.5]
# resolution = 1000
#
# kon = 0.52631579
# khill_on = 5.26315789
# koff = 0.0073
# kpos = 1
# khill_pos = 10
#
# for khill_on in np.linspace(0, 10, 10):
#     xpoints, ypoints = evaluate(func, xrange, yrange, resolution, (kon, 1, khill_on, koff, kpos, khill_pos, 1))
#     plt.scatter(xpoints, ypoints, s=0.1)
# # plt.xlim(xrange[0], xrange[1])
# # plt.ylim(yrange[0], yrange[1])
# plt.scatter(cyt_uniform, mem_uniform)
# plt.show()

"""
Saturation model 2

"""


def func(x, y, kon, koff, kpos, kpos2, e):
    return kon * x - koff * y + kpos * x * (y ** 1) * np.exp(-kpos2 * (y ** e))


koff = 0.0073
kon = 0.0421342134213
kpos_max = 1
kpos2_max = 1
e = 4

kpos_vals = np.linspace(-3, 2, 1000)
kpos2_vals = np.linspace(-5, 2, 1000)
g = list(itertools.product(range(1000), range(1000)))


# res_uniform = np.zeros([1000, 1000])
# res_polar = np.zeros([1000, 1000])


def func_uniform(i):
    kpos = kpos_vals[g[i][0]]
    kpos2 = kpos2_vals[g[i][1]]
    res = mse(cyt_uniform, mem_uniform, func, [0, 20], (kon, koff, np.exp(kpos), np.exp(kpos2), e))
    return g[i][0], g[i][1], res


def func_polar(i):
    kpos = kpos_vals[g[i][0]]
    kpos2 = kpos2_vals[g[i][1]]
    res = mse(cyt_polar, mem_polar, func, [0, 20], (kon, koff, np.exp(kpos), np.exp(kpos2), e))
    return g[i][0], g[i][1], res


def tidy(res):
    res2 = np.zeros([1000, 1000])
    for i in res:
        res2[int(i[0]), int(i[1])] = i[2]
    return res2


# res_uniform = tidy(
#     np.array(
#         Parallel(n_jobs=multiprocessing.cpu_count(), verbose=10)(delayed(func_uniform)(i) for i in range(1000000))))
# np.savetxt('res_uniform_log_e4.txt', res_uniform)
#
# res_polar = tidy(
#     np.array(Parallel(n_jobs=multiprocessing.cpu_count(), verbose=10)(delayed(func_polar)(i) for i in range(1000000))))
# np.savetxt('res_polar_log_e4.txt', res_polar)


#####

filename_uniform = 'res_uniform_log_e3.txt'
filename_polar = 'res_polar_log_e3.txt'
res_uniform = np.loadtxt(filename_uniform)
res_polar = np.loadtxt(filename_polar)

res_uniform[np.isnan(res_uniform)] = 1000
a_uniform = np.where(res_uniform == np.min(res_uniform))
kpos_opt_uniform = kpos_vals[a_uniform[0]]
kpos2_opt_uniform = kpos2_vals[a_uniform[1]]

res_polar[np.isnan(res_polar)] = 1000
a_polar = np.where(res_polar == np.min(res_polar))
kpos_opt_polar = kpos_vals[a_polar[0]]
kpos2_opt_polar = kpos2_vals[a_polar[1]]


####

# def func(name, figname):
#     plt.close()
#     res = np.loadtxt(name)
#     res[np.isnan(res)] = 1000
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     im = ax.imshow(res, cmap='gist_rainbow', vmin=0, vmax=1, origin='lower', extent=(-5, 2, -3, 2))
#     ax.set_ylabel('log10(k_pos)')
#     ax.set_xlabel('log10(k_sat)')
#     fig.colorbar(im)
#     ax.set_aspect(7 / 5)
#
#     a = np.where(res == np.min(res))
#     print(np.min(res))
#     plt.scatter(kpos2_vals[a[1]], kpos_vals[a[0]], c='k')
#
#     plt.savefig(figname)
#     # plt.show()


# # func('res_uniform_log_e1.txt', 'res_uniform_log_e1.png')
# # func('res_uniform_log_e2.txt', 'res_uniform_log_e2.png')
# func('res_uniform_log_e3.txt', 'res_uniform_log_e3.png')
# # func('res_uniform_log_e4.txt', 'res_uniform_log_e4.png')
#
# # func('res_polar_log_e1.txt', 'res_polar_log_e1.png')
# # func('res_polar_log_e2.txt', 'res_polar_log_e2.png')
# func('res_polar_log_e3.txt', 'res_polar_log_e3.png')
# # func('res_polar_log_e4.txt', 'res_polar_log_e4.png')


#####

print(np.exp(kpos_opt_uniform))
print(np.exp(kpos2_opt_uniform))
print(np.exp(kpos_opt_polar))
print(np.exp(kpos2_opt_polar))

plt.scatter(cyt_uniform, mem_uniform, c='none', edgecolors='c', label='Uniform')
xpoints, ypoints = evaluate(func, [0, 0.25], [0, 5], 5000,
                            (0.0421342134213, 0.0073, 0.18, np.exp(kpos2_opt_uniform), 3))
plt.scatter(xpoints, ypoints, s=0.1, c='c')

plt.scatter(cyt_polar, mem_polar, c='none', edgecolors='m', label='Polarised')
xpoints, ypoints = evaluate(func, [0, 0.25], [0, 5], 5000,
                            (0.0421342134213, 0.0073, 0.18, np.exp(kpos2_opt_polar), 3))
plt.scatter(xpoints, ypoints, s=0.1, c='m')

plt.scatter(ring_cyt, ring_mem, c='none', edgecolors='g', label='C56S')
xpoints, ypoints = evaluate(func, [0, 0.25], [0, 5], 5000, (0.0421342134213, 0.0073, 0, 0, 0))
plt.scatter(xpoints, ypoints, s=0.1, c='g')

sns.despine()
plt.xlim(0, 0.065)
plt.ylim(0, 2.9)
plt.xlabel('Cytoplasmic PAR-2')
plt.ylabel('Cortical PAR-2')
# plt.xlim(left=0)
# plt.ylim(bottom=0)
plt.legend(frameon=False)
plt.show()

"""
Linear model: RING mutant
For finding on rate

"""


# def func(x, y, kon, koff):
#     return kon * x - koff * y
#
#
# koff = 0.0073
# kon_max = 0.1
#
# kon_vals = np.linspace(0, kon_max, 10000)
# res = np.zeros([10000])
#
# for i, kon in enumerate(kon_vals):
#     res[i] = mse(ring_cyt, ring_mem, func, [-5, 20], (kon, koff))
#
# print(kon_vals[np.argmin(res)])
# plt.plot(kon_vals, res)
# plt.show()
