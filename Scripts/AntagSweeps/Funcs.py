import matplotlib

matplotlib.use('Agg')

import numpy as np
import os
import copy
import multiprocessing
import itertools
from scipy import interpolate
import matplotlib.pyplot as plt
import glob
import shutil

"""
Parameter sweep

"""


# class ParamSweep:
#     def __init__(self, base_model, p1, p2, p1_range, p2_range, log, cores, resolution, n_iterations, direc):
#         """
#
#         :param base_model:
#         :param p1: name of parameter 1
#         :param p2: name of parameter 2
#         :param p1_range: range for parameter 1 - (lower, upper)
#         :param p2_range: range for parameter 2 - (lower, upper)
#         :param log: if TRUE will apply log10 to parameter values
#         :param cores: number of cores on machine to use in parallel
#         :param resolution: specifies distance of simulations on grid space
#         :param n_iterations: number of iterations
#         :param direc: directory to save results. Must already exist
#         """
#
#         self.base_model = base_model
#         self.p1 = p1
#         self.p2 = p2
#         self.p1_range = p1_range
#         self.p2_range = p2_range
#         self.log = log
#         self.cores = cores
#         self.resolution = resolution
#         self.n_iterations = n_iterations
#         self.direc = direc
#
#         # Results
#         self.res = None
#         self.p1_vals = None
#         self.p2_vals = None
#
#     def run(self, p1_val, p2_val):
#         """
#         Run simulations for given p1 and p2 values, save result
#
#         Final states:
#         -1: unpolarised, P dominant
#         0:  polarised
#         1:  unpolarised, A dominant
#
#
#         """
#
#         name = '%s_%s' % (p1_val, p2_val)
#
#         # Test if simulation has already been performed
#         if not os.path.exists(self.direc + '/' + name + '.txt'):
#             print(name)
#
#             # Set up model
#             model = copy.deepcopy(self.base_model)
#
#             # Remove antagonism
#             kAP = model.kAP
#             kPA = model.kPA
#             model.kAP = 0
#             model.kPA = 0
#
#             # Initial equilibration (no antagonism)
#             for t in range(10000):
#                 model.react()
#
#             # Polarise
#             model.am *= 2 * np.r_[np.ones([model.xsteps // 2]), np.zeros([model.xsteps // 2])]
#             model.pm *= 2 * np.r_[np.zeros([model.xsteps // 2]), np.ones([model.xsteps // 2])]
#
#             # Add antagonism
#             model.kAP = kAP
#             model.kPA = kPA
#
#             # Set parameters
#             if self.log:
#                 setattr(model, self.p1, 10 ** p1_val)
#                 setattr(model, self.p2, 10 ** p2_val)
#             else:
#                 setattr(model, self.p1, p1_val)
#                 setattr(model, self.p2, p2_val)
#
#             # Simulation
#             for t in range(int(model.Tmax / model.deltat)):
#                 model.react()
#                 model.time = (t + 1) * model.deltat
#
#             # Save state
#             if sum(model.am > model.pm) == len(model.am):
#                 np.savetxt(self.direc + '/' + name + '.txt', [1])
#             elif sum(model.am > model.pm) == 0:
#                 np.savetxt(self.direc + '/' + name + '.txt', [-1])
#             else:
#                 np.savetxt(self.direc + '/' + name + '.txt', [0])
#
#     def compile_res(self):
#         """
#         Import results for all simulations that have been run
#         Assume all other values by nearest neighbour interpolation
#
#         """
#
#         self.res = np.zeros([len(self.p1_vals), len(self.p2_vals)])
#
#         # Import data
#         for i, p1_val in enumerate(self.p1_vals):
#             for j, p2_val in enumerate(self.p2_vals):
#                 name = '%s_%s' % (p1_val, p2_val)
#                 if os.path.exists(self.direc + '/' + name + '.txt'):
#                     self.res[i, j] = np.loadtxt(self.direc + '/' + name + '.txt')
#                 else:
#                     self.res[i, j] = np.nan
#
#         # Interpolate nans (nearest neighbour)
#         x = np.arange(0, self.res.shape[1])
#         y = np.arange(0, self.res.shape[0])
#         array = np.ma.masked_invalid(self.res)
#         xx, yy = np.meshgrid(x, y)
#         x1 = xx[~array.mask]
#         y1 = yy[~array.mask]
#         newarr = array[~array.mask]
#         self.res = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method='nearest')
#
#     def sweep(self):
#         """
#         Perform parameter sweep in p1-p2 parameter space
#         Save results and make figures with each iteration
#
#         """
#
#         for iteration in range(self.n_iterations):
#             print(iteration)
#
#             if iteration == 0:
#                 n_sims = self.resolution
#
#                 # Initial broad parameter sweep
#                 self.p1_vals = np.linspace(self.p1_range[0], self.p1_range[1], n_sims)
#                 self.p2_vals = np.linspace(self.p2_range[0], self.p2_range[1], n_sims)
#                 sims_array = np.array(list(itertools.product(self.p1_vals, self.p2_vals)))
#                 pool = multiprocessing.Pool(self.cores)
#                 pool.starmap(self.run, iter(sims_array))
#
#                 # Compile results
#                 self.compile_res()
#
#             else:
#                 n_sims = self.resolution * (n_sims - 1) + 1
#
#                 # Find boundary regions
#                 a = np.nonzero(np.diff(self.res, axis=0))
#                 b = np.nonzero(np.diff(self.res, axis=1))
#                 xpoints = np.r_[a[1], b[1]]
#                 ypoints = np.r_[a[0], b[0]]
#
#                 # Parameters for next iteration
#                 self.p1_vals = np.linspace(self.p1_range[0], self.p1_range[1], n_sims)
#                 self.p2_vals = np.linspace(self.p2_range[0], self.p2_range[1], n_sims)
#                 run_bool = np.zeros([n_sims, n_sims])
#                 for x, y in zip(xpoints, ypoints):
#                     run_bool[y * self.resolution:y * self.resolution + (self.resolution + 1),
#                     x * self.resolution:x * self.resolution + (self.resolution + 1)] = 1
#                 sims_array = np.array(list(itertools.product(self.p1_vals, self.p2_vals)))[run_bool.flatten() == 1, :]
#
#                 # Run next iteration
#                 pool = multiprocessing.Pool(self.cores)
#                 pool.starmap(self.run, iter(sims_array))
#
#                 # Compile results
#                 self.compile_res()
#
#             # Save figures
#             self.im_fig()
#             self.scatter_fig()
#
#         # Save results
#         if os.path.isdir(self.direc + '/Results/'):
#             shutil.rmtree(self.direc + '/Results/')
#         os.mkdir(self.direc + '/Results/')
#         np.savetxt(self.direc + '/Results/Res.txt', self.res)
#         np.savetxt(self.direc + '/Results/%s_vals.txt' % self.p1, self.p1_vals)
#         np.savetxt(self.direc + '/Results/%s_vals.txt' % self.p2, self.p2_vals)
#
#     def im_fig(self):
#         """
#         Phase space plot (shows nearest-neighbour interpolated values for simulations that haven't been run)
#
#         """
#
#         fig, ax = plt.subplots()
#         ax.imshow(self.res.T, cmap='bwr', origin='lower', alpha=0.5,
#                   extent=(self.p1_range[0], self.p1_range[1], self.p2_range[0], self.p2_range[1]))
#         # ax.set_ylabel(r'$\log_{10}(k_{PA})$')
#         # ax.set_xlabel(r'$\log_{10}(k_{AP})$')
#         if self.log:
#             ax.set_xlabel('log(%s)' % self.p1)
#             ax.set_ylabel('log(%s)' % self.p2)
#         else:
#             ax.set_xlabel(self.p1)
#             ax.set_ylabel(self.p2)
#         fig.set_size_inches(4, 4)
#         fig.tight_layout()
#         fig.savefig(self.direc + '/im_fig.png', dpi=300)
#
#     def scatter_fig(self):
#         """
#         Scatter plot of all simulations that have been run
#
#         """
#
#         fig, ax = plt.subplots()
#         res = np.zeros([len(glob.glob(self.direc + '/*.txt'))])
#         p1 = np.zeros([len(glob.glob(self.direc + '/*.txt'))])
#         p2 = np.zeros([len(glob.glob(self.direc + '/*.txt'))])
#         for i, a in enumerate(glob.glob(self.direc + '/*.txt')):
#
#             # Import data
#             res[i] = np.loadtxt(a)
#             name = a.split('/')[-1][:-4]
#             p1[i] = float(name.split('_')[0])
#             p2[i] = float(name.split('_')[1])
#
#             # Plot point
#             if res[i] == -1:
#                 ax.scatter(p1[i], p2[i], c='b', s=2)
#             elif res[i] == 1:
#                 ax.scatter(p1[i], p2[i], c='r', s=2)
#             else:
#                 ax.scatter(p1[i], p2[i], c='0.5', s=2)
#
#         # ax.set_ylabel(r'$\log_{10}(k_{PA})$')
#         # ax.set_xlabel(r'$\log_{10}(k_{AP})$')
#         if self.log:
#             ax.set_xlabel('log(%s)' % self.p1)
#             ax.set_ylabel('log(%s)' % self.p2)
#         else:
#             ax.set_xlabel(self.p1)
#             ax.set_ylabel(self.p2)
#         fig.set_size_inches(4, 4)
#         fig.tight_layout()
#         fig.savefig(self.direc + '/scatter_fig.png', dpi=300)
