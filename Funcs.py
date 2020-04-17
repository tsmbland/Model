import matplotlib as mpl

mpl.use('Agg')

import numpy as np
import os
import multiprocessing
import matplotlib.pyplot as plt
import itertools

"""
PDE solver

"""


def diffusion(concs, dx):
    d = concs[np.r_[np.array(range(1, len(concs))), len(concs) - 1]] - 2 * concs + concs[
        np.r_[0, np.array(range(len(concs) - 1))]]
    return d / (dx ** 2)


def pdeRK(dxdt, X0, Tmax, deltat, t_eval):
    """
    Adapted from Hubatsch 2019

    :param dxdt:
    :param X0:
    :param Tmax:
    :param deltat:
    :param t_eval:

    To do:
    - way to break if system stabilises before time limit is reached

    """

    time = 0

    # Adaptive step size parameters
    atol = 0.000001
    rtol = 0.000001

    # 5TH ORDER RK COEFFICIENTS for Dormand-Prince
    a21, a31, a32, a41, a42, a43 = 1 / 5, 3 / 40, 9 / 40, 44 / 45, -56 / 15, 32 / 9
    a51, a52, a53, a54 = 19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729
    a61, a62, a63 = 9017 / 3168, -355 / 33, 46732 / 5247
    a64, a65 = 49 / 176, -5103 / 18656
    a71, a72, a73, a74 = 35 / 384, 0, 500 / 1113, 125 / 192
    a75, a76 = -2187 / 6784, 11 / 84

    b1, b2, b3, b4, b5 = 35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784
    b6, b7 = 11 / 84, 0

    bs1, bs2, bs3, bs4 = 5179 / 57600, 0, 7571 / 16695, 393 / 640
    bs5, bs6, bs7 = -92097 / 339200, 187 / 2100, 1 / 40

    # Set up
    X = X0
    testsoln = dxdt(X)
    nvars = len(testsoln)
    spatial_points = len(testsoln[0])
    t_stored = np.zeros([len(t_eval)])
    X_stored = [np.zeros([len(t_eval), spatial_points]) for i in range(nvars)]
    stored_times = 0

    # Run
    updated = None
    terminate = False
    while not terminate:
        if time > Tmax:
            terminate = True

        # Calculate increments for RK45
        if (time == 0) or not updated:
            X1 = dxdt(X)
        else:
            X1 = X7

        X2 = dxdt([X[i] + deltat * (a21 * X1[i]) for i in range(nvars)])
        X3 = dxdt([X[i] + deltat * (a31 * X1[i] + a32 * X2[i]) for i in range(len(X))])
        X4 = dxdt([X[i] + deltat * (a41 * X1[i] + a42 * X2[i] + a43 * X3[i]) for i in range(nvars)])
        X5 = dxdt([X[i] + deltat * (a51 * X1[i] + a52 * X2[i] + a53 * X3[i] + a54 * X4[i]) for i in range(nvars)])
        X6 = dxdt([X[i] + deltat * (a61 * X1[i] + a62 * X2[i] + a63 * X3[i] + a64 * X4[i] + a65 * X5[i]) for i in
                   range(nvars)])
        X7 = dxdt([X[i] + deltat * (a71 * X1[i] + a73 * X3[i] + a74 * X4[i] + a75 * X5[i] + a76 * X6[i]) for i in
                   range(nvars)])

        # Update concentrations using A1-A6 and P1-P6, coefficient for A7 and P7 is 0.
        Xn_new = [X[i] + deltat * (b1 * X1[i] + b3 * X3[i] + b4 * X4[i] + b5 * X5[i] + b6 * X6[i]) for i in
                  range(nvars)]  # b2/7=0

        # Compute difference between fourth and fifth order
        deltaXnerr = [max(abs(
            (b1 - bs1) * X1[i] + (b3 - bs3) * X3[i] + (b4 - bs4) * X4[i] + (b5 - bs5) * X5[i] + (b6 - bs6) * X6[
                i] - bs7 * X7[i])) for i in range(nvars)]  # b7 is zero

        # Get maximum concentrations for An and Pn
        yXn = [np.maximum(max(abs(Xn_new[i])), max(abs(X[i]))) for i in range(nvars)]

        # Get error scale, combining relative and absolute error
        scaleXn = [atol + yXn[i] * rtol for i in range(nvars)]

        # Compute total error as norm of maximum errors for each species scaled by the error scale
        errs = [(deltaXnerr[i] / scaleXn[i]) ** 2 for i in range(nvars)]
        totalerror = np.sqrt(sum(errs) / nvars)

        # Compute new timestep
        dtnew = 0.8 * deltat * abs(1 / totalerror) ** (1 / 5)

        # Upper and lower bound for timestep to avoid changing too fast
        if dtnew > 10 * deltat:
            dtnew = 10 * deltat
        elif dtnew < deltat / 5:
            dtnew = deltat / 5

        # Set timestep for next round
        deltat = dtnew

        # Accept step if error is on the order of error scale or below
        if totalerror < 1:
            time += deltat
            X = Xn_new
            updated = True

            # Save
            if sum(time > t_eval) > stored_times:
                t_stored[stored_times] = time
                for i in range(nvars):
                    X_stored[i][stored_times, :] = X[i]
                stored_times += 1
        else:
            updated = False

    return X_stored, t_stored


"""
Evaluate 

"""


def evaluate(func, xrange, yrange, iterations, resolution=100, args=()):
    """

    :param func: function to solve, takes x and y as arguments
    :param xrange: (lower x, upper x)
    :param yrange: (lower y, upper y)
    :param iterations:
    :param resolution: resolution of grid space for first iteration
    :param args: extra arguments for bistability
    :return:
    """

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
            xvals = xrange[0] + sims_array_ind[0] * (xrange[1] - xrange[0]) / (n_sims - 1)
            yvals = yrange[0] + sims_array_ind[1] * (yrange[1] - yrange[0]) / (n_sims - 1)

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


"""
Parameter sweep

"""


class ParamSweep2D:
    def __init__(self, func, p1_vals, p2_vals, direc, log=False, parallel=True, cores=None):
        """
        Runs func with all combinations of p1_vals and p2_vals
        Func must save the result to be used for later analysis

        :param func:
        :param p1_vals:
        :param p2_vals:
        :param direc:
        :param log:
        :param parallel:
        :param cores:
        """
        self.func = func
        self.p1_vals = p1_vals
        self.p2_vals = p2_vals
        self.cores = cores
        self.direc = direc
        self.log = log
        self.parallel = parallel

    def single(self, p1_val, p2_val):

        # Run function
        if self.log:
            self.func(10 ** p1_val, 10 ** p2_val)
        else:
            self.func(p1_val, p2_val)

    def run(self):

        # Parameter combinations
        sims_array = np.array(list(itertools.product(self.p1_vals, self.p2_vals)))

        # Run
        if self.parallel:
            pool = multiprocessing.Pool(self.cores)
            pool.starmap(self.single, iter(sims_array))
        else:
            for k in iter(sims_array):
                self.single(*k)


class Bifurcation2D:
    def __init__(self, func, p1_range, p2_range, resolution0, resolution_step, n_iterations, direc,
                 parallel=False, colours=None, crange=None, show_boundaries=False,
                 cores=None):
        """

        :param func: function - takes 2 parameters, returns an integer (must not be zero)
        :param p1_range: range for parameter 1 - (lower, upper)
        :param p2_range: range for parameter 2 - (lower, upper)
        :param cores: number of cores on machine to use in parallel
        :param resolution0: n x n points on initial grid
        :param resolution_step: how much resolution increases with each time step
        :param n_iterations: number of iterations
        :param direc: directory to save results. Must already exist


        To do:
        - save to single text file rather than loads of small files
        - less importing
        - use integer arrays? e.g. use 0 instead of np.nan
        - currently breaks if no boundaries are found
        - also save condensed res array

        """

        self.func = func
        self.p1_range = p1_range
        self.p2_range = p2_range
        self.cores = cores
        self.resolution0 = resolution0
        self.resolution_step = resolution_step
        self.n_iterations = n_iterations
        self.direc = direc
        self.parallel = parallel
        self.colours = colours
        self.crange = crange
        self.show_boundaries = show_boundaries

        # Results
        self.res = None
        self.n_sims = None

    def single(self, p1_val, p2_val):
        """
        Evaluate for given p1 and p2 values, save result

        """

        name = '%s_%s' % (p1_val, p2_val)

        # If bistability not already performed for these parameters
        if not os.path.exists(self.direc + '/Res/' + name + '.txt'):
            # Run function
            state = self.func(p1_val, p2_val)

            # Save state
            with open(self.direc + '/Res/' + name + '.txt', 'w') as f:
                f.write('%d' % state)

    def run(self):
        """
        Run parameter sweep

        """

        # Create results folder
        if not os.path.exists(self.direc + '/Res'):
            os.mkdir(self.direc + '/Res')

        for iteration in range(self.n_iterations):
            print(iteration)

            if iteration == 0:
                self.n_sims = self.resolution0
                run_bool = np.ones([self.n_sims, self.n_sims])

            else:
                self.n_sims = self.resolution_step * (self.n_sims - 1) + 1

                # Find boundary regions
                a = np.nonzero(np.nan_to_num(np.diff(self.res, axis=0)))
                b = np.nonzero(np.nan_to_num(np.diff(self.res, axis=1)))
                c = np.nonzero(np.nan_to_num(self.res[:-1, :-1] - self.res[1:, 1:]))
                xpoints = np.r_[a[0], b[0], c[0]]
                ypoints = np.r_[a[1], b[1], c[1]]
                run_bool = np.zeros([self.n_sims, self.n_sims])
                for x, y in zip(xpoints, ypoints):
                    run_bool[x * self.resolution_step:x * self.resolution_step + (self.resolution_step + 1),
                    y * self.resolution_step:y * self.resolution_step + (self.resolution_step + 1)] = 1

            # Parameter combinations
            sims_array_ind = np.nonzero(run_bool)
            p1vals = self.p1_range[0] + sims_array_ind[0] * (self.p1_range[1] - self.p1_range[0]) / (self.n_sims - 1)
            p2vals = self.p2_range[0] + sims_array_ind[1] * (self.p2_range[1] - self.p2_range[0]) / (self.n_sims - 1)

            # Run
            if self.parallel:
                pool = multiprocessing.Pool(self.cores)
                pool.starmap(self.single, iter(np.c_[p1vals, p2vals]))
            else:
                for k in iter(np.c_[p1vals, p2vals]):
                    self.single(*k)

            # Compile results
            self.res = np.nan * np.zeros([self.n_sims, self.n_sims])
            for r in range(len(p1vals)):
                name = '%s_%s' % (p1vals[r], p2vals[r])
                with open(self.direc + '/Res/' + name + '.txt') as g:
                    self.res[sims_array_ind[0][r], sims_array_ind[1][r]] = g.read()

            # Safety check: If any nans are adjacent to more than one type (in all directions), test them
            if iteration != 0:
                j = 1
                while j != 0:

                    # Compare each nan value to all neighbours (this is quite slow, quicker way?)
                    x = np.dstack((self.res[:-2, :-2], self.res[:-2, 1:-1], self.res[:-2, 2:], self.res[1:-1, :-2],
                                   self.res[1:-1, 2:], self.res[2:, :-2], self.res[2:, 1:-1], self.res[2:, 2:]))
                    mx = np.nanmax(x, axis=2)
                    mn = np.nanmin(x, axis=2)
                    run_bool = np.zeros([self.n_sims, self.n_sims])
                    run_bool[1:-1, 1:-1] = (mx == mx) * (mx != mn) * (self.res[1:-1, 1:-1] != self.res[1:-1, 1:-1])
                    sims_array_ind = np.nonzero(run_bool)
                    j = len(sims_array_ind[0])

                    # Parameter combinations (if any)
                    if j != 0:
                        p1vals = self.p1_range[0] + sims_array_ind[0] * (self.p1_range[1] - self.p1_range[0]) / (
                            self.n_sims - 1)
                        p2vals = self.p2_range[0] + sims_array_ind[1] * (self.p2_range[1] - self.p2_range[0]) / (
                            self.n_sims - 1)

                        # Run
                        if self.parallel:
                            pool = multiprocessing.Pool(self.cores)
                            pool.starmap(self.single, iter(np.c_[p1vals, p2vals]))
                        else:
                            for k in iter(np.c_[p1vals, p2vals]):
                                self.single(*k)

                        # Compile results
                        for r in range(len(p1vals)):
                            name = '%s_%s' % (p1vals[r], p2vals[r])
                            with open(self.direc + '/Res/' + name + '.txt') as g:
                                self.res[sims_array_ind[0][r], sims_array_ind[1][r]] = g.read()

        # Interpolate nans - flood fill algorithm
        self.res = np.nan_to_num(self.res).astype(int)
        o = np.argwhere(self.res == 0)
        while len(o) != 0:
            pos = o[0]
            fillval = findzone(self.res, pos[0], pos[1])
            floodfill(self.res, pos[0], pos[1], fillval)
            o = np.argwhere(self.res == 0)

        # Save final results
        np.savetxt(self.direc + '/Results.txt', self.res, fmt='%i')
        np.savetxt(self.direc + '/p1_vals.txt', np.linspace(self.p1_range[0], self.p1_range[1], self.n_sims))
        np.savetxt(self.direc + '/p2_vals.txt', np.linspace(self.p2_range[0], self.p2_range[1], self.n_sims))

        # Figures
        self.im_fig()

    def im_fig(self):
        """
        Phase space plot (shows nearest-neighbour interpolated values for parameter sets that haven't been evaluated)

        """

        # Set up
        plt.close()
        fig, ax = plt.subplots()

        # Colours
        if self.colours is not None:
            cmap = mpl.colors.ListedColormap(list(self.colours.values()))
            norm = mpl.colors.BoundaryNorm(list(self.colours.keys()) + [100], cmap.N)
            ax.imshow(self.res.T, origin='lower', aspect='auto', cmap=cmap, norm=norm,
                      extent=(self.p1_range[0], self.p1_range[1], self.p2_range[0], self.p2_range[1]))
        else:
            ax.imshow(self.res.T, origin='lower', aspect='auto', vmin=self.crange[0], vmax=self.crange[1],
                      extent=(self.p1_range[0], self.p1_range[1], self.p2_range[0], self.p2_range[1]))

        # Boundaries
        if self.show_boundaries:
            a = np.nonzero(np.diff(self.res, axis=0))
            b = np.nonzero(np.diff(self.res, axis=1))
            c = np.nonzero(self.res[:-1, :-1] - self.res[1:, 1:])
            xpoints = np.r_[a[0], b[0], c[0]] + 1
            ypoints = np.r_[a[1], b[1], c[1]] + 1
            xpoints = self.p1_range[0] + xpoints * (self.p1_range[1] - self.p1_range[0]) / self.n_sims
            ypoints = self.p2_range[0] + ypoints * (self.p2_range[1] - self.p2_range[0]) / self.n_sims
            ax.scatter(xpoints, ypoints, s=0.1, c='k')

        # Save figure
        ax.set_xlim(self.p1_range[0], self.p1_range[1])
        ax.set_ylim(self.p2_range[0], self.p2_range[1])
        fig.set_size_inches(4, 4)
        fig.tight_layout()
        fig.savefig(self.direc + '/im_fig.png', dpi=300)

        # Close
        plt.close()


def floodfill(array, x, y, newval):
    """

    :param array:
    :param x:
    :param y:
    :param newval:
    :return:
    """
    oldval = array[x, y]
    if oldval == newval:
        return
    array[x, y] = newval
    Q = [(x, y)]
    while len(Q) != 0:
        x, y = Q.pop(0)

        if x > 0:
            if array[x - 1, y] == oldval:
                array[x - 1, y] = newval
                Q.append((x - 1, y))
        if x < len(array[:, 0]) - 1:
            if array[x + 1, y] == oldval:
                array[x + 1, y] = newval
                Q.append((x + 1, y))
        if y > 0:
            if array[x, y - 1] == oldval:
                array[x, y - 1] = newval
                Q.append((x, y - 1))
        if y < len(array[0, :]) - 1:
            if array[x, y + 1] == oldval:
                array[x, y + 1] = newval
                Q.append((x, y + 1))


def findzone(array, x, y):
    """

    :param array:
    :param x:
    :param y:
    :return:
    """

    tested = np.zeros(array.shape)
    tested[x, y] = 1
    Q = [(x, y)]
    val = 0
    while val == 0:
        x, y = Q.pop(0)

        vals = []
        if x > 0:
            if tested[x - 1, y] == 0:
                tested[x - 1, y] = 1
                vals.append(array[x - 1, y])
                Q.append((x - 1, y))
        if x < len(array[:, 0]) - 1:
            if tested[x + 1, y] == 0:
                tested[x + 1, y] = 1
                vals.append(array[x + 1, y])
                Q.append((x + 1, y))
        if y > 0:
            if tested[x, y - 1] == 0:
                tested[x, y - 1] = 1
                vals.append(array[x, y - 1])
                Q.append((x, y - 1))
        if y < len(array[0, :]) - 1:
            if tested[x, y + 1] == 0:
                tested[x, y + 1] = 1
                vals.append(array[x, y + 1])
                Q.append((x, y + 1))

        if len(vals) != 0:
            val = max(vals)
    return val


"""
Misc

"""


# def direcslist(dest, levels=0):
#     """
#
#     Gives a list of directories in a given directory (full path)
#     Excludes directories that contain !
#
#     :param dest:
#     :return:
#     """
#     lis = glob.glob('%s/*/' % dest)
#     for level in range(levels):
#         newlis = []
#         for e in lis:
#             newlis.extend(glob.glob('%s/*/' % e))
#         lis = newlis
#     lis = [x[:-1] for x in lis if '!' not in x]
#     return lis
#
#
# def animate(self, direc):
#     """
#     direc: directory to results
#
#     """
#
#     times = np.loadtxt(direc + '/times.txt')
#     A = np.loadtxt(direc + '/A.txt')
#     P = np.loadtxt(direc + '/P.txt')
#
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     plt.subplots_adjust(bottom=0.25, wspace=0.5)
#     axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
#     sframe = Slider(axframe, 'Time (s)', 0, times[-1], valinit=0, valfmt='%d')
#
#     def update(i):
#         ax.clear()
#         tpoint = np.argmin(abs(times - int(i)))
#         a = A[tpoint, :]
#         p = P[tpoint, :]
#         ax.plot(np.linspace(0, self.L, self.xsteps), a, c='tab:red')
#         ax.plot(np.linspace(0, self.L, self.xsteps), p, c='tab:blue')
#         ax.set_ylim(bottom=0)
#         ax.set_ylabel('Cortical concentration (a.u.)')
#         ax.set_xlabel('Position (μm)')
#
#     sframe.on_changed(update)
#     plt.show()
