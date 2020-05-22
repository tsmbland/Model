import numpy as np
import os
import multiprocessing
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

"""
Evaluate functions

"""


def evaluate(func, xrange, yrange, resolution0=1000, iterations=1, resolution_step=2, args=()):
    """
    Solve a function in the form 0 = f(x, y)

    :param func: function to solve (0 = f(x,y,*args))
    :param xrange: (lower x, upper x)
    :param yrange: (lower y, upper y)
    :param resolution0: resolution of grid space for first iteration
    :param iterations: if > 1 will iteratively refine solution by increasing resolution by factor of resolution_step
    :param resolution_step:
    :param args: additional arguments for func
    :return: two 1D arrays of points corresponding to x and y coordinates. NB these are UNORDERED, plot using
        plt.scatter(xpoints,ypoints)


    """

    for iteration in range(iterations):

        if iteration == 0:

            # Set x and y values
            n_sims = resolution0
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
            n_sims = resolution_step * (n_sims - 1) + 1

            run_bool = np.zeros([n_sims, n_sims])
            for x, y in zip(xpoints, ypoints):
                run_bool[x * resolution_step:x * resolution_step + (resolution_step + 1),
                y * resolution_step:y * resolution_step + (resolution_step + 1)] = 1
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
PDE solver

"""


def diffusion(concs, dx=1):
    """
    Simulate single diffusion time step, with reflective boundary conditions

    :param concs: 1D array of concentrations across space
    :param dx: optional, spatial distance between points
    :return:
    """

    d = concs[np.r_[np.array(range(1, len(concs))), len(concs) - 1]] - 2 * concs + concs[
        np.r_[0, np.array(range(len(concs) - 1))]]
    return d / (dx ** 2)


def pdeRK(dxdt, X0, Tmax, deltat, t_eval, killfunc=None, stabilitycheck=False):
    """

    Function for solving system of PDEs using adaptive Runge-Kutta method
    Adapted from Hubatsch 2019

    :param dxdt: function, takes list of 1D arrays corresponding to concentrations over space, returns list of gradients
    :param X0: starting conditions, list of 1D arrays
    :param Tmax: maximum time
    :param deltat: initial timestep
    :param t_eval: array of timepoints to save
    :param killfunc: optional func, takes same input as dxdt, integration will terminate when func returns True
    :param stabilitycheck: if True, will terminate integration when system changes by less that 1% per minute

    Returns
    X: final state
    time: final time
    X_stored: saved states according to t_eval
    t_stored: times corresponding to X_stored. Will npt be exactly the same as t_eval but should be close

    """

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
    X_stored = [np.zeros([len(t_eval), spatial_points]) for _ in range(nvars)]
    n_stored_times = 0

    # Run
    time = 0
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
        # sometimes see "RuntimeWarning: divide by zero encountered in double_scalars". Need to look into
        dtnew = 0.8 * deltat * abs(1 / totalerror) ** (1 / 5)

        # Upper and lower bound for timestep to avoid changing too fast
        if dtnew > 10 * deltat:
            dtnew = 10 * deltat
        elif dtnew < deltat / 5:
            dtnew = deltat / 5

        # Compute max percentage change
        change = max(max(abs(X[i] - Xn_new[i]) / Xn_new[i]) * (60 / dtnew) for i in range(nvars))

        # Set timestep for next round
        deltat = dtnew

        # Accept step if error is on the order of error scale or below
        if totalerror < 1:
            time += deltat
            X = Xn_new
            updated = True

            # Store
            if sum(time > t_eval) > n_stored_times:
                t_stored[n_stored_times] = time
                for i in range(nvars):
                    X_stored[i][n_stored_times, :] = X[i]
                n_stored_times += 1

            # Kill function
            if killfunc is not None:
                brk = killfunc(X)
                if brk:
                    break

            # Check stability:
            if stabilitycheck:
                if change < 0.001:
                    break

        else:
            updated = False

    return X, time, X_stored, t_stored


"""
Parameter space

"""


class ParamSpace1D:
    def __init__(self, func, p_vals, direc, parallel=False, cores=None):
        """
        Runs func with all combinations of p1_vals and p2_vals, saves results to csv file

        :param func: function, takes two parameter values, returns a float
        :param p_vals: array of parameter values
        :param direc: directory to save results
        :param parallel: if True, will run in parallel using number of cores specified
        :param cores: number of cores to use, if parallel=True

        To do:
        - check that this works
        - ability to stop and resume
        - make more like 2D class, i.e. ability to explore boundaries

        """
        self.func = func
        self.p_vals = p_vals
        self.cores = cores
        self.direc = direc
        self.parallel = parallel

    def single_eval(self, p_val):

        # Run function
        res = self.func(p_val)

        # Save results
        with open(self.direc + '/Res.csv', 'a') as f:
            f.write("{:.12f}".format(p_val) + ',' + str(res) + '\n')

    def run(self):

        # Make results directory, if it doesn't already exist
        if not os.path.isdir(self.direc):
            os.mkdir(self.direc)

        # Run
        if self.parallel:
            pool = multiprocessing.Pool(self.cores)
            pool.map(self.single_eval, self.p_vals)
        else:
            for p in iter(self.p_vals):
                self.single_eval(p)


class ParamSpace2D:
    def __init__(self, func, p1_range, p2_range, resolution0, direc, resolution_step=2, n_iterations=1,
                 explore_boundaries=True, parallel=False, cores=None, crange=None, cmap=None, save_fig=True):
        """

        Class to perform parameter space sweeps

        Run by calling run function
        - performs analysis, saves results with each iteration to .csv files, and saves final figure
        - progress is saved, if interrupted can be resumed without loss by calling the run function again

        Computation parameters
        :param func: function - takes 2 parameters, returns an integer (must not be zero)
        :param p1_range: range for parameter 1 (lower, upper)
        :param p2_range: range for parameter 2 (lower, upper)
        :param resolution0: n x n points on initial grid
        :param resolution_step: how much resolution increases with each iteration
        :param n_iterations: number of iterations
        :param explore_boundaries: if True, will focus parameter search to regions where func result varies
        :param parallel: if True, will run in parallel using number of cores specified
        :param cores: number of cores on machine to use in parallel

        # Saving parameters
        :param direc: directory to save results. Directory must already exist

        # Figure parameters
        :param save_fig: if True, will save figure
        :param crange: (lowest value, highest value)
        :param cmap: if None, pyplot will use 'viridis'

        To do:
        - adjust so it can be run from multiple machines simultaneously
        - test

        """

        # Computation
        self.func = func
        self.p1_range = p1_range
        self.p2_range = p2_range
        self.resolution0 = resolution0
        self.resolution_step = resolution_step
        self.n_iterations = n_iterations
        self.explore_boundaries = explore_boundaries
        self.parallel = parallel

        if cores is None:
            self.cores = multiprocessing.cpu_count()
        else:
            self.cores = cores

        # Saving
        self.direc = direc
        self.save_fig = save_fig

        # Figure
        self.crange = crange
        self.cmap = cmap

        # Results
        self.iteration = None
        self.res = None
        self.n_sims = None

    def single_eval(self, p1val_p2val):
        """
        Single function call for given p1 and p2 values, save result

        """

        # Run function
        state = self.func(*[float(i) for i in p1val_p2val.split(',')])

        # Save state
        with open(self.direc + '/' + str(self.iteration) + '.csv', 'a') as f:
            f.write(p1val_p2val + ',' + str(state) + '\n')

    def batch_eval(self, pcombs):
        """
        Evaluate parameter sets in bulk
        pcombs is list of strings 'p1val,p2val'

        """
        if self.parallel:
            pool = multiprocessing.Pool(self.cores)
            pool.map(self.single_eval, pcombs)
        else:
            for k in iter(pcombs):
                self.single_eval(k)

    def import_res(self):
        """
        Import all results from current iteration, load into self.res

        """

        with open(self.direc + '/' + str(self.iteration) + '.csv') as g:
            for line in g:
                p1, p2, val = line[:-1].split(',')
                xind = ((float(p1) - self.p1_range[0]) * (self.n_sims - 1)) / (self.p1_range[1] - self.p1_range[0])
                yind = ((float(p2) - self.p2_range[0]) * (self.n_sims - 1)) / (self.p2_range[1] - self.p2_range[0])
                if '.' in val:
                    self.res[round(xind), round(yind)] = float(val)
                else:
                    self.res[round(xind), round(yind)] = int(val)

    def run(self):
        """
        Run algorithm, save figure

        """

        # Make results directory, if it doesn't already exist
        if not os.path.isdir(self.direc):
            os.mkdir(self.direc)

        for iteration in range(self.n_iterations):
            print(iteration)
            self.iteration = iteration

            # First iteration, initial grid
            if self.iteration == 0:
                self.n_sims = self.resolution0
                run_bool = np.ones([self.n_sims, self.n_sims])

            # Subsequent iteration, explore boundary regions
            else:
                self.n_sims = self.resolution_step * (self.n_sims - 1) + 1

                if self.explore_boundaries:
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

                else:
                    run_bool = np.ones([self.n_sims, self.n_sims])

            # Parameter combinations
            sims_array_ind = np.nonzero(run_bool)
            p1vals = self.p1_range[0] + sims_array_ind[0] * (self.p1_range[1] - self.p1_range[0]) / (self.n_sims - 1)
            p2vals = self.p2_range[0] + sims_array_ind[1] * (self.p2_range[1] - self.p2_range[0]) / (self.n_sims - 1)
            pcombs = ["{:.12f}".format(p1vals[i]) + ',' + "{:.12f}".format(p2vals[i]) for i in range(len(p1vals))]

            # Remove parameters already tested (if algorithm run before)
            if os.path.isfile(self.direc + '/' + str(self.iteration) + '.csv'):
                with open(self.direc + '/' + str(self.iteration) + '.csv') as f:
                    for line in f:
                        p = line.split(',')[0] + ',' + line.split(',')[1]
                        if p in pcombs:
                            pcombs.remove(p)

            # Carry over combinations from previous iteration
            if self.iteration != 0:
                with open(self.direc + '/' + str(self.iteration - 1) + '.csv') as f:
                    with open(self.direc + '/' + str(self.iteration) + '.csv', 'a') as g:
                        for line in f:
                            p = line.split(',')[0] + ',' + line.split(',')[1]
                            if p in pcombs:
                                pcombs.remove(p)
                                g.write(line)

            # Run
            self.batch_eval(pcombs)

            # Import results
            self.res = np.nan * np.zeros([self.n_sims, self.n_sims])
            self.import_res()

            # Safety check: test any untested points directly adjacent to boundaries
            # (often required if resolution0 or resolution_step are too small)
            if self.explore_boundaries and self.iteration != 0:
                j = 1
                while j != 0:

                    # Compare each nan value to all neighbours
                    # (this is quite slow, quicker way?)
                    # throws warning: "All-NaN axis encountered" - this is fine
                    res_padded = np.nan * np.zeros([self.n_sims + 2, self.n_sims + 2])
                    res_padded[1:-1, 1:-1] = self.res
                    x = np.dstack(
                        (res_padded[:-2, :-2], res_padded[:-2, 1:-1], res_padded[:-2, 2:], res_padded[1:-1, :-2],
                         res_padded[1:-1, 2:], res_padded[2:, :-2], res_padded[2:, 1:-1], res_padded[2:, 2:]))

                    mx = np.nanmax(x, axis=2)
                    mn = np.nanmin(x, axis=2)
                    run_bool = (mx == mx) * (mx != mn) * (self.res != self.res)
                    sims_array_ind = np.nonzero(run_bool)
                    j = len(sims_array_ind[0])

                    # Parameter combinations (if any)
                    if j != 0:
                        p1vals = self.p1_range[0] + sims_array_ind[0] * (self.p1_range[1] - self.p1_range[0]) / (
                            self.n_sims - 1)
                        p2vals = self.p2_range[0] + sims_array_ind[1] * (self.p2_range[1] - self.p2_range[0]) / (
                            self.n_sims - 1)
                        pcombs = ["{:.12f}".format(p1vals[i]) + ',' + "{:.12f}".format(p2vals[i]) for i in
                                  range(len(p1vals))]

                        # Run
                        self.batch_eval(pcombs)

                        # Compile results
                        self.import_res()

        # Interpolate missing values
        if self.explore_boundaries:

            # If no boundary (i.e. uniform parameter space), reload output from first iteration
            if np.sum(~np.isnan(self.res)) == 0:
                with open(self.direc + '/0.csv') as f:
                    self.res[:, :] = int(f.readline().split(',')[2])

            # Interpolate nans by flood fill algorithm
            else:
                o = np.argwhere(self.res != self.res)
                while len(o) != 0:
                    pos = o[0]
                    fillval = findzone(self.res, pos[0], pos[1])
                    floodfill(self.res, pos[0], pos[1], fillval)
                    o = np.argwhere(self.res != self.res)

        # Save figure
        if self.save_fig:
            fig, ax = self.fig()
            fig.set_size_inches(4, 4)
            fig.tight_layout()
            plt.savefig(self.direc + '/fig.png', dpi=300)
            plt.close()

        # Specify int/float format
        f = (self.res % 1) == 0

        # Save row by row
        with open(self.direc + '/Res.txt', 'w') as fh:
            for i, row in enumerate(self.res):
                formats = f[i]
                line = ' '.join(
                    "{:.0f}".format(value) if formats[j] else "{:.12f}".format(value) for j, value in enumerate(row))
                fh.write(line + '\n')

    def fig(self):
        """
        Parameter space plot

        """

        # Set up
        plt.close()
        fig, ax = plt.subplots()

        # Colours
        ax.imshow(self.res.T, origin='lower', aspect='auto', vmin=self.crange[0], vmax=self.crange[1],
                  cmap=self.cmap, extent=(self.p1_range[0], self.p1_range[1], self.p2_range[0], self.p2_range[1]))

        # Figure adjustments
        ax.set_xlim(self.p1_range[0], self.p1_range[1])
        ax.set_ylim(self.p2_range[0], self.p2_range[1])
        return fig, ax


def floodfill(array, x, y, newval):
    """
    Queue based flood fill algorithm
    Edits array in place

    :param array:
    :param x: x position
    :param y: y position
    :param newval: new value
    :return:
    """

    array[x, y] = newval
    Q = [(x, y)]
    while len(Q) != 0:
        x, y = Q.pop(0)

        if x > 0:
            if np.isnan(array[x - 1, y]):
                array[x - 1, y] = newval
                Q.append((x - 1, y))
        if x < len(array[:, 0]) - 1:
            if np.isnan(array[x + 1, y]):
                array[x + 1, y] = newval
                Q.append((x + 1, y))
        if y > 0:
            if np.isnan(array[x, y - 1]):
                array[x, y - 1] = newval
                Q.append((x, y - 1))
        if y < len(array[0, :]) - 1:
            if np.isnan(array[x, y + 1]):
                array[x, y + 1] = newval
                Q.append((x, y + 1))


def findzone(array, x, y):
    """
    For 2D phase space class.
    Finds zone corresponding to a given location
    Based on flood fill algorithm

    :param array:
    :param x: x position
    :param y: y position
    :return: zone value
    """

    tested = np.zeros(array.shape)
    tested[x, y] = 1
    Q = [(x, y)]
    val = np.nan
    while val != val:
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


def animatePAR(direc):
    """
    Display an animation for PAR model

    direc: directory to results

    """

    times = np.loadtxt(direc + '/times.txt')
    A = np.loadtxt(direc + '/A.txt')
    P = np.loadtxt(direc + '/P.txt')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.25, wspace=0.5)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Time (s)', 0, times[-1], valinit=0, valfmt='%d')

    def update(i):
        ax.clear()
        tpoint = np.argmin(abs(times - int(i)))
        a = A[tpoint, :]
        p = P[tpoint, :]
        ax.plot(a, c='tab:red')
        ax.plot(p, c='tab:blue')
        ax.set_ylim(bottom=0)
        ax.set_ylabel('Cortical concentration (a.u.)')
        ax.set_xlabel('Position (μm)')

    sframe.on_changed(update)
    plt.show()
