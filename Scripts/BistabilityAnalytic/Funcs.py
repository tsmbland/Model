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
import seaborn as sns

"""
Simple PAR model: nullclines

"""


def evaluate(func, xrange, yrange, resolution, args):
    res = func(np.tile(np.linspace(xrange[0], xrange[1], resolution), (resolution, 1)).T,
               np.tile(np.linspace(yrange[0], yrange[1], resolution), (resolution, 1)), *args)
    res_sign = np.sign(res)

    a = np.nonzero((np.diff(res_sign, axis=0, append=[res_sign[-1, :]]) + np.diff(res_sign.T, axis=0,
                                                                                  append=[res_sign[:, -1]]).T).astype(
        bool).astype(int))

    xpoints = xrange[0] + (a[0] / resolution) * (xrange[1] - xrange[0])
    ypoints = yrange[0] + (a[1] / resolution) * (yrange[1] - yrange[0])
    return xpoints, ypoints


def intersections(func1, func2, xrange, yrange, resolution, dist):
    # Solve func 1
    res1 = func1(np.tile(np.linspace(xrange[0], xrange[1], resolution), (resolution, 1)).T,
                 np.tile(np.linspace(yrange[0], yrange[1], resolution), (resolution, 1)))
    res1_sign = np.sign(res1)
    a = np.diff(res1_sign, axis=0, append=[res1_sign[-1, :]]) + np.diff(res1_sign.T, axis=0,
                                                                        append=[res1_sign[:, -1]]).T
    a = a.astype(bool).astype(int)

    # Solve func 2
    res2 = func2(np.tile(np.linspace(xrange[0], xrange[1], resolution), (resolution, 1)).T,
                 np.tile(np.linspace(yrange[0], yrange[1], resolution), (resolution, 1)))
    res2_sign = np.sign(res2)
    b = np.diff(res2_sign, axis=0, append=[res2_sign[-1, :]]) + np.diff(res2_sign.T, axis=0,
                                                                        append=[res2_sign[:, -1]]).T
    b = b.astype(bool).astype(int)

    # Add, threshold
    c = ((a + b) == 2).astype(int)
    # plt.imshow((a + b).T, origin='lower', extent=(xrange[0], xrange[1], yrange[0], yrange[1]))

    # Coordinates
    coors = np.where(c == 1)
    xpoints = xrange[0] + (coors[0] / resolution) * (xrange[1] - xrange[0])
    ypoints = yrange[0] + (coors[1] / resolution) * (yrange[1] - yrange[0])

    # Pairwise distance
    d = ((((np.tile(xpoints, (len(xpoints), 1)) - np.tile(xpoints, (len(xpoints), 1)).T) ** 2) + (
        (np.tile(ypoints, (len(ypoints), 1)) - np.tile(ypoints, (len(ypoints), 1)).T) ** 2)) ** 0.5)
    d = d < dist

    # Link
    points = np.arange(len(xpoints))
    for i in range(len(xpoints)):
        points[d[i, :]] = points[i]

    # Average
    xpoints2 = [np.mean(xpoints[points == e]) for e in np.unique(points)]
    ypoints2 = [np.mean(ypoints[points == e]) for e in np.unique(points)]
    # plt.scatter(xpoints2, ypoints2)

    return xpoints2, ypoints2


class Model:
    def __init__(self, konA, koffA, konP, koffP, kposA, kposP, kAP, kPA, ePneg, eAneg, pA, pP, psi):
        self.konA = konA
        self.koffA = koffA
        self.konP = konP
        self.koffP = koffP
        self.kposA = kposA
        self.kposP = kposP
        self.kAP = kAP
        self.kPA = kPA
        self.ePneg = ePneg
        self.eAneg = eAneg
        self.pA = pA
        self.pP = pP
        self.psi = psi

    def funcA(self, a, p):
        a_cyt = self.pA - self.psi * a
        return (self.konA * a_cyt) + (self.kposA * a_cyt * a) - (self.koffA * a) - (self.kAP * a * (p ** self.ePneg))

    def funcP(self, a, p):
        p_cyt = self.pP - self.psi * p
        return (self.konP * p_cyt) + (self.kposP * p_cyt * p) - (self.koffP * p) - (self.kPA * p * (a ** self.eAneg))


base_model = Model(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=2,
                   psi=0.1, pA=1, pP=1)


# model = copy.deepcopy(base_model)
# model.kAP = 0.1
# model.kPA = 0.1
#
# intersections(model.funcA, model.funcP, (0, 10), (0, 10), 1000, 0.1)
# plt.show()



class AntagSweep:
    def __init__(self, base_model, xrange, yrange, p1_range, p2_range, cores, resolution, n_iterations, direc,
                 parallel=True):
        self.base_model = base_model
        self.p1_range = p1_range
        self.p2_range = p2_range
        self.xrange = xrange
        self.yrange = yrange
        self.cores = cores
        self.resolution = resolution
        self.n_iterations = n_iterations
        self.direc = direc
        self.parallel = parallel

        # Results
        self.res = None
        self.p1_vals = None
        self.p2_vals = None

    def run(self, kAP, kPA):
        name = '%s_%s' % (kAP, kPA)
        print(name)

        # Test if simulation has already been performed
        if not os.path.exists(self.direc + '/' + name + '.txt'):
            # Set up model
            model = copy.deepcopy(self.base_model)

            # Set antagonism rates
            model.kAP = 10 ** kAP
            model.kPA = 10 ** kPA

            # Analysis
            xpoints, ypoints = intersections(model.funcA, model.funcP, self.xrange, self.yrange, 5000, 0.1)

            # Save
            if len(xpoints) > 1:
                np.savetxt(self.direc + '/' + name + '.txt', [1])
            else:
                np.savetxt(self.direc + '/' + name + '.txt', [0])

    def compile_res(self):
        self.res = np.zeros([len(self.p1_vals), len(self.p2_vals)])

        # Import data
        for i, kAP in enumerate(self.p1_vals):
            for j, kPA in enumerate(self.p2_vals):
                name = '%s_%s' % (kAP, kPA)
                if os.path.exists(self.direc + '/' + name + '.txt'):
                    self.res[i, j] = np.loadtxt(self.direc + '/' + name + '.txt')
                else:
                    self.res[i, j] = np.nan

        # Interpolate nans
        x = np.arange(0, self.res.shape[1])
        y = np.arange(0, self.res.shape[0])
        array = np.ma.masked_invalid(self.res)
        xx, yy = np.meshgrid(x, y)
        x1 = xx[~array.mask]
        y1 = yy[~array.mask]
        newarr = array[~array.mask]
        self.res = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method='nearest')

    def sweep(self):
        for iteration in range(self.n_iterations):
            print(iteration)

            if iteration == 0:
                n_sims = self.resolution

                # Initial broad parameter sweep
                self.p1_vals = np.linspace(self.p1_range[0], self.p1_range[1], n_sims)
                self.p2_vals = np.linspace(self.p2_range[0], self.p2_range[1], n_sims)
                sims_array = np.array(list(itertools.product(self.p1_vals, self.p2_vals)))

                if self.parallel:
                    pool = multiprocessing.Pool(self.cores)
                    pool.starmap(self.run, iter(sims_array))

                else:
                    for i in iter(sims_array):
                        self.run(*i)

                # Compile results
                self.compile_res()

            else:
                n_sims = self.resolution * (n_sims - 1) + 1

                # Find boundary regions
                a = np.nonzero(np.diff(self.res, axis=0))
                b = np.nonzero(np.diff(self.res, axis=1))
                xpoints = np.r_[a[1], b[1]]
                ypoints = np.r_[a[0], b[0]]

                # Parameters for next iteration
                self.p1_vals = np.linspace(self.p1_range[0], self.p1_range[1], n_sims)
                self.p2_vals = np.linspace(self.p2_range[0], self.p2_range[1], n_sims)
                run_bool = np.zeros([n_sims, n_sims])
                for x, y in zip(xpoints, ypoints):
                    run_bool[y * self.resolution:y * self.resolution + (self.resolution + 1),
                    x * self.resolution:x * self.resolution + (self.resolution + 1)] = 1
                sims_array = np.array(list(itertools.product(self.p1_vals, self.p2_vals)))[run_bool.flatten() == 1, :]

                # Run next iteration
                if self.parallel:
                    pool = multiprocessing.Pool(self.cores)
                    pool.starmap(self.run, iter(sims_array))
                else:
                    for i in iter(sims_array):
                        self.run(*i)

                # Compile results
                self.compile_res()

                # Save figures
                self.im_fig()
                self.scatter_fig()

        # Save results
        if os.path.isdir(self.direc + '/Results/'):
            shutil.rmtree(self.direc + '/Results/')
        os.mkdir(self.direc + '/Results/')
        np.savetxt(self.direc + '/Results/Res.txt', self.res)
        np.savetxt(self.direc + '/Results/kAP_vals.txt', self.p1_vals)
        np.savetxt(self.direc + '/Results/kPA_vals.txt', self.p2_vals)

    def im_fig(self):
        fig, ax = plt.subplots()
        ax.imshow(self.res.T, cmap='bwr', origin='lower', alpha=0.5,
                  extent=(self.p1_range[0], self.p1_range[1], self.p2_range[0], self.p2_range[1]))
        ax.set_ylabel('log10(kPA)')
        ax.set_xlabel('log10(kAP)')
        fig.set_size_inches(4, 4)
        fig.tight_layout()
        fig.savefig(self.direc + '/im_fig.png', dpi=300)

    def scatter_fig(self):
        fig, ax = plt.subplots()
        res = np.zeros([len(glob.glob(self.direc + '/*.txt'))])
        kAP = np.zeros([len(glob.glob(self.direc + '/*.txt'))])
        kPA = np.zeros([len(glob.glob(self.direc + '/*.txt'))])
        for i, a in enumerate(glob.glob(self.direc + '/*.txt')):
            res[i] = np.loadtxt(a)
            name = a.split('/')[-1][:-4]
            kAP[i] = float(name.split('_')[0])
            kPA[i] = float(name.split('_')[1])

            if res[i] == -1:
                ax.scatter(kPA[i], kAP[i], c='b', s=2)
            elif res[i] == 1:
                ax.scatter(kPA[i], kAP[i], c='r', s=2)
            else:
                ax.scatter(kPA[i], kAP[i], c='0.5', s=2)

        ax.set_ylabel('log10(kPA)')
        ax.set_xlabel('log10(kAP)')
        fig.set_size_inches(4, 4)
        fig.tight_layout()
        fig.savefig(self.direc + '/scatter_fig.png', dpi=300)
