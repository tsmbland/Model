# 180723


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import animation
import seaborn as sns
from joblib import Parallel, delayed
import os
import pickle
import itertools
import random
import multiprocessing
import sys

sns.set()
sns.set_style("ticks")


############################### MODELS #############################


class Model0:
    """Originial model with mass action antagonism and no positive feedback"""

    class Params:
        def __init__(self, Da, Dp, konA, koffA, konP, koffP, kAP, kPA, ePneg, eAneg, pA, pP, pgen):
            # Diffusion
            self.Da = Da  # um2 s-1
            self.Dp = Dp  # um2 s-1

            # Membrane exchange
            self.konA = konA  # um s-1
            self.koffA = koffA  # s-1
            self.konP = konP  # um s-1
            self.koffP = koffP  # s-1

            # Antagonism
            self.kAP = kAP  # um2 s-1
            self.kPA = kPA  # um4 s-1
            self.ePneg = ePneg
            self.eAneg = eAneg

            # Pools
            self.pA = pA  # um-3
            self.pP = pP  # um-3

            # Misc
            self.L = pgen.L  # um
            self.xsteps = pgen.xsteps
            self.psi = pgen.psi  # um-1
            self.Tmax = pgen.Tmax  # s
            self.deltat = pgen.deltat  # s

            # Equilibration
            self.Aeqmin = pgen.Aeqmin
            self.Aeqmax = pgen.Aeqmax
            self.Peqmin = pgen.Peqmin
            self.Peqmax = pgen.Peqmax

    def __init__(self, p):
        self.aco = np.zeros([p.xsteps])
        self.pco = np.zeros([p.xsteps])

    def update_aco(self, p):
        diff = diffusion(self.aco, p.Da, p)
        off = (p.koffA * self.aco)
        on = (p.konA * (p.pA - p.psi * np.mean(self.aco)))
        ant = p.kAP * (self.pco ** p.ePneg) * self.aco
        self.aco += ((diff + on - off - ant) * p.deltat)

    def update_pco(self, p):
        diff = diffusion(self.pco, p.Dp, p)
        off = (p.koffP * self.pco)
        on = (p.konP * (p.pP - p.psi * np.mean(self.pco)))
        ant = p.kPA * (self.aco ** p.eAneg) * self.pco
        self.pco += ((diff + on - off - ant) * p.deltat)


class Model1:
    """Model with Hill functions for antagonism"""

    class Params:
        def __init__(self, Da, Dp, konA, koffA, konP, koffP, kAP, kPA, kAneg, kPneg, ePneg, eAneg, pA, pP, pgen):
            # Diffusion
            self.Da = Da  # um2 s-1
            self.Dp = Dp  # um2 s-1

            # Membrane exchange
            self.konA = konA  # um s-1
            self.koffA = koffA  # s-1
            self.konP = konP  # um s-1
            self.koffP = koffP  # s-1

            # Antagonism
            self.kAP = kAP  # s-1
            self.kPA = kPA  # s-1
            self.kAneg = kAneg  # um-2
            self.kPneg = kPneg  # um-2
            self.ePneg = ePneg  #
            self.eAneg = eAneg  #

            # Pools
            self.pA = pA  # um-3
            self.pP = pP  # um-3

            # Misc
            self.L = pgen.L  # um
            self.xsteps = pgen.xsteps
            self.psi = pgen.psi  # um-1
            self.Tmax = pgen.Tmax  # s
            self.deltat = pgen.deltat  # s

            # Equilibration
            self.Aeqmin = pgen.Aeqmin
            self.Aeqmax = pgen.Aeqmax
            self.Peqmin = pgen.Peqmin
            self.Peqmax = pgen.Peqmax

    def __init__(self, p):
        self.aco = np.zeros([p.xsteps])
        self.pco = np.zeros([p.xsteps])

    def update_aco(self, p):
        diff = diffusion(self.aco, p.Da, p)
        off = (p.koffA * self.aco)
        on = (p.konA * (p.pA - p.psi * np.mean(self.aco)))
        ant = p.kAP * self.aco * (self.pco ** p.ePneg) / (
            (p.kPneg ** p.ePneg) + (self.pco ** p.ePneg))
        self.aco += ((diff + on - off - ant) * p.deltat)

    def update_pco(self, p):
        diff = diffusion(self.pco, p.Dp, p)
        off = (p.koffP * self.pco)
        on = (p.konP * (p.pP - p.psi * np.mean(self.pco)))
        ant = p.kPA * self.pco * (self.aco ** p.eAneg) / (
            (p.kAneg ** p.eAneg) + (self.aco ** p.eAneg))
        self.pco += (diff + on - off - ant) * p.deltat


class Model2:
    """Model with Hill functions for antagonism and positive feedback"""

    class Params:
        def __init__(self, Da, Dp, konA1, koffA, konP1, koffP, kAP, kPA, kAneg, kPneg, ePneg, eAneg, konA2, konP2,
                     kApos, kPpos, eApos, ePpos, pA, pP, pgen):
            # Diffusion
            self.Da = Da  # um2 s-1
            self.Dp = Dp  # um2 s-1

            # Membrane exchange
            self.konA1 = konA1  # um s-1
            self.koffA = koffA  # s-1
            self.konP1 = konP1  # um s-1
            self.koffP = koffP  # s-1

            # Antagonism
            self.kAP = kAP  # s-1
            self.kPA = kPA  # s-1
            self.kAneg = kAneg  # um-2
            self.kPneg = kPneg  # um-2
            self.ePneg = ePneg  #
            self.eAneg = eAneg  #

            # Positive feedback
            self.konA2 = konA2  # um2 s-1
            self.konP2 = konP2  # um2 s-1
            self.kApos = kApos  # um-2
            self.kPpos = kPpos  # um-2
            self.eApos = eApos
            self.ePpos = ePpos

            # Pools
            self.pA = pA  # um-3
            self.pP = pP  # um-3

            # Misc
            self.L = pgen.L  # um
            self.xsteps = pgen.xsteps
            self.psi = pgen.psi  # um-1
            self.Tmax = pgen.Tmax  # s
            self.deltat = pgen.deltat  # s

            # Equilibration
            self.Aeqmin = pgen.Aeqmin
            self.Aeqmax = pgen.Aeqmax
            self.Peqmin = pgen.Peqmin
            self.Peqmax = pgen.Peqmax

    def __init__(self, p):
        self.aco = np.zeros([p.xsteps])
        self.pco = np.zeros([p.xsteps])

    def update_aco(self, p):
        diff = diffusion(self.aco, p.Da, p)
        off = (p.koffA * self.aco)
        on = (p.konA1 * (p.pA - p.psi * np.mean(self.aco))) + (
            p.konA2 * (p.pA - p.psi * np.mean(self.aco)) * (self.aco ** p.eApos) / (
                (p.kApos ** p.eApos) + (self.aco ** p.eApos)))
        ant = p.kAP * self.aco * (self.pco ** p.ePneg) / (
            (p.kPneg ** p.ePneg) + (self.pco ** p.ePneg))
        self.aco += ((diff + on - off - ant) * p.deltat)

    def update_pco(self, p):
        diff = diffusion(self.pco, p.Dp, p)
        off = (p.koffP * self.pco)
        on = (p.konP1 * (p.pP - p.psi * np.mean(self.pco))) + (
            p.konP2 * (p.pP - p.psi * np.mean(self.pco)) * (self.pco ** p.ePpos) / (
                (p.kPpos ** p.ePpos) + (self.pco ** p.ePpos)))
        ant = p.kPA * self.pco * (self.aco ** p.eAneg) / (
            (p.kAneg ** p.eAneg) + (self.aco ** p.eAneg))
        self.pco += (diff + on - off - ant) * p.deltat


class Model3:
    """Model with antagonism on positive feedback"""

    class Params:
        def __init__(self, Da, Dp, konA1, koffA, konP1, koffP, kAneg, kPneg, ePneg, eAneg, konA2, konP2, kApos, kPpos,
                     eApos, ePpos, pA, pP, pgen):
            # Diffusion
            self.Da = Da  # um2 s-1
            self.Dp = Dp  # um2 s-1

            # Membrane exchange
            self.konA1 = konA1  # um s-1
            self.koffA = koffA  # s-1
            self.konP1 = konP1  # um s-1
            self.koffP = koffP  # s-1

            # Positive feedback
            self.konA2 = konA2  # um2 s-1
            self.konP2 = konP2  # um2 s-1
            self.kApos = kApos  # um2
            self.kPpos = kPpos  # um2
            self.eApos = eApos
            self.ePpos = ePpos

            # Antagonism
            self.kPneg = kPneg  # um2
            self.kAneg = kAneg  # um2
            self.ePneg = ePneg
            self.eAneg = eAneg

            # Pools
            self.pA = pA  # um-3
            self.pP = pP  # um-3

            # Misc
            self.L = pgen.L  # um
            self.xsteps = pgen.xsteps
            self.psi = pgen.psi  # um-1
            self.Tmax = pgen.Tmax  # s
            self.deltat = pgen.deltat  # s

            # Equilibration
            self.Aeqmin = pgen.Aeqmin
            self.Aeqmax = pgen.Aeqmax
            self.Peqmin = pgen.Peqmin
            self.Peqmax = pgen.Peqmax

    def __init__(self, p):
        self.aco = np.zeros([p.xsteps])
        self.pco = np.zeros([p.xsteps])

    def update_aco(self, p):
        diff = diffusion(self.aco, p.Da, p)
        off = (p.koffA * self.aco)
        int_on = (p.konA1 * (p.pA - p.psi * np.mean(self.aco)))
        pf_on = (p.konA2 * (p.pA - p.psi * np.mean(self.aco)) * (self.aco ** p.eApos) / (
            (p.kApos ** p.eApos) + (self.aco ** p.eApos)))
        ant = (p.kPneg ** p.ePneg) / ((p.kPneg ** p.ePneg) + (self.pco ** p.ePneg))
        self.aco += ((diff + int_on - off + pf_on * ant) * p.deltat)

    def update_pco(self, p):
        diff = diffusion(self.pco, p.Dp, p)
        off = (p.koffP * self.pco)
        int_on = (p.konP1 * (p.pP - p.psi * np.mean(self.pco)))
        pf_on = (p.konP2 * (p.pP - p.psi * np.mean(self.pco)) * (self.pco ** p.ePpos) / (
            (p.kPpos ** p.ePpos) + (self.pco ** p.ePpos)))
        ant = (p.kAneg ** p.eAneg) / ((p.kAneg ** p.eAneg) + (self.aco ** p.eAneg))
        self.pco += ((diff + int_on - off + pf_on * ant) * p.deltat)


class Model4:
    """Model with pPAR state switch"""

    class Params:
        def __init__(self, Da, Dp, konA, koffA, konP, koffP, kAP, kPA, kPneg1, kPneg2, kAneg1, kAneg2, ePneg, eAneg,
                     kApos, kPpos, eApos, ePpos, pA, pP, pgen):
            # Diffusion
            self.Da = Da  # um2 s-1
            self.Dp = Dp  # um2 s-1

            # Membrane exchange
            self.konA = konA  # um s-1
            self.koffA = koffA  # s-1
            self.konP = konP  # um s-1
            self.koffP = koffP  # s-1

            # Antagonism
            self.kAP = kAP  # s-1
            self.kPA = kPA  # s-1
            self.kPneg1 = kPneg1  # um2
            self.kPneg2 = kPneg2  # um2
            self.kAneg1 = kAneg1  # um2
            self.kAneg2 = kAneg2  # um2
            self.ePneg = ePneg
            self.eAneg = eAneg

            # Positive feedback
            self.kApos = kApos  # um2
            self.kPpos = kPpos  # um2
            self.eApos = eApos
            self.ePpos = ePpos

            # Pools
            self.pA = pA  # um-3
            self.pP = pP  # um-3

            # Misc
            self.L = pgen.L  # um
            self.xsteps = pgen.xsteps
            self.psi = pgen.psi  # um-1
            self.Tmax = pgen.Tmax  # s
            self.deltat = pgen.deltat  # s

            # Equilibration
            self.Aeqmin = pgen.Aeqmin
            self.Aeqmax = pgen.Aeqmax
            self.Peqmin = pgen.Peqmin
            self.Peqmax = pgen.Peqmax

    def __init__(self, p):
        self.aco = np.zeros([p.xsteps])
        self.pco = np.zeros([p.xsteps])

    def update_aco(self, p):
        diff = diffusion(self.aco, p.Da, p)
        off = (p.koffA * self.aco)
        on = (p.konA * (p.pA - p.psi * np.mean(self.aco)))
        k = p.kPneg1 + p.kPneg2 * ((self.aco ** p.eApos) / ((p.kApos ** p.eApos) + (self.aco ** p.eApos)))
        ant = p.kAP * self.aco * ((self.pco ** p.ePneg) / ((k ** p.ePneg) + (self.pco ** p.ePneg)))
        self.aco += ((diff + on - off - ant) * p.deltat)

    def update_pco(self, p):
        diff = diffusion(self.pco, p.Dp, p)
        off = (p.koffP * self.pco)
        on = (p.konP * (p.pP - p.psi * np.mean(self.pco)))
        k = p.kAneg1 + p.kAneg2 * (self.pco ** p.ePpos) / ((p.kPpos ** p.ePpos) + (self.pco ** p.ePpos))
        ant = p.kPA * self.pco * (self.aco ** p.eAneg) / ((k ** p.eAneg) + (self.aco ** p.eAneg))
        self.pco += ((diff + on - off - ant) * p.deltat)


class MiscParams:
    """
    General parameters shared between all models

    """

    def __init__(self, L, xsteps, psi, Tmax, deltat, Aeqmin, Aeqmax, Peqmin, Peqmax):
        # Misc
        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s

        # Equilibration
        self.Aeqmin = Aeqmin
        self.Aeqmax = Aeqmax
        self.Peqmin = Peqmin
        self.Peqmax = Peqmax


def diffusion(concs, coeff, p):
    if hasattr(concs, '__len__'):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (p.L / p.xsteps)
    else:
        diff = 0
    return diff


def flow(concs, flows, p):
    fl = (concs[np.append(np.array(range(1, p.xsteps)), [p.xsteps - 2])] * flows[
        np.append(np.array(range(1, p.xsteps)), [p.xsteps - 2])] - concs[
              np.append([1], np.array(range(p.xsteps - 1)))] * flows[
              np.append([1], np.array(range(p.xsteps - 1)))]) / (2 * (p.L / p.xsteps))
    return fl


#########################  ALGORITHM FUNCTIONS ##########################


class Res:
    # Structure containing simulation results and the parameters used

    def __init__(self, p):
        self.p = p
        self.aco = np.zeros([int(p.Tmax / p.deltat) + 1, p.xsteps])
        self.pco = np.zeros([int(p.Tmax / p.deltat) + 1, p.xsteps])


def run_model(m, p):
    x = m(p)
    res = Res(p)

    # Equilibrate
    x.aco = equilibrate_aco(m, p)
    x.pco = equilibrate_pco(m, p)
    res.aco[0, :] = x.aco
    res.pco[0, :] = x.pco

    # Run model
    for t in range(int(p.Tmax / p.deltat)):
        x.update_aco(p)
        x.update_pco(p)
        res.aco[t + 1, :] = x.aco
        res.pco[t + 1, :] = x.pco

    return res


def equilibrate_aco(m, p):
    x = m(p)
    for t in range(5000):
        x.update_aco(p)
        x.aco[:int(p.xsteps * p.Aeqmin)] = 0
        x.aco[int(p.xsteps * p.Aeqmax):] = 0
    return x.aco


def equilibrate_pco(m, p):
    x = m(p)
    for t in range(5000):
        x.update_pco(p)
        x.pco[:int(p.xsteps * p.Peqmin)] = 0
        x.pco[int(p.xsteps * p.Peqmax):] = 0
    return x.pco


def nullcline(m, p):
    """
    need to paralelise

    :param m:
    :param p:
    :return:
    """

    p.Da = 0
    p.Dp = 0
    x = m(p)

    xlim = 2
    ylim = 2
    res1 = 100
    res2 = 20
    tmax = 1000

    x2 = np.linspace(0, xlim, res1)
    y1 = np.linspace(0, ylim, res1)
    x1 = np.zeros([res1, res2])
    y2 = np.zeros([res1, res2])

    # P nullcline
    x.aco = x2
    for a in range(res2):
        x.pco = np.linspace(0, ylim, res2)[a] * np.ones([res1])
        for t in range(tmax):
            x.update_pco(p)
        y2[:, a] = x.pco

    # A nullcline
    x.pco = y1
    for a in range(res2):
        x.aco = np.linspace(0, xlim, res2)[a] * np.ones([res1])
        for t in range(tmax):
            x.update_aco(p)
        x1[:, a] = x.aco

    return x1, y1, x2, y2


# def solve_equation(func):
#     """
#     Function to solve equations in to form f(x,y) = 0
#
#     :param func:
#     :return:
#     """
#
#     xlim = 2
#     ylim = 2
#     res1 = 10000
#     res2 = 10
#     tmax = 10000
#
#     x2 = np.linspace(0, xlim, res1)
#     y1 = np.linspace(0, ylim, res1)
#     x1 = np.zeros([res1, res2])
#     y2 = np.zeros([res1, res2])
#
#     # 1
#     x = x2
#     for a in range(res2):
#         y = np.linspace(0, ylim, res2)[a] * np.ones([res1])
#         for t in range(tmax):
#             y += func(x, y)
#         y2[:, a] = y
#
#     # 2
#     y = y1
#     for a in range(res2):
#         x = np.linspace(0, xlim, res2)[a] * np.ones([res1])
#         for t in range(tmax):
#             x += func(x, y)
#         x1[:, a] = x
#
#     return x1, y1, x2, y2


def trajectory(m, p, a0, p0):
    # Init
    x = m(p)
    aco = np.zeros([int(p.Tmax / p.deltat) + 1])
    pco = np.zeros([int(p.Tmax / p.deltat) + 1])
    x.aco = a0
    x.pco = p0

    # Run
    for t in range(int(p.Tmax / p.deltat)):
        x.update_aco(p)
        x.update_pco(p)
        aco[t + 1] = x.aco
        pco[t + 1] = x.pco

    return aco, pco


def savedata(res, jobid, subjobid, simid, compression):
    """

    :param res: results from simulation
    :param jobid:
    :param subjobid:
    :param compression:
    :return:
    """

    # Data compression
    if compression == 0:
        pass
    elif compression == 1:
        res.aco = np.asarray([res.aco[-1, :], ])
        res.pco = np.asarray([res.pco[-1, :], ])

    # Create job directory
    if not os.path.isdir('../PolarityData/%s' % ('{0:04}'.format(jobid))):
        os.makedirs('../PolarityData/%s' % ('{0:04}'.format(jobid)))

    # Create subjob directory
    if not os.path.isdir('../PolarityData/%s/%s' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid))):
        os.makedirs('../PolarityData/%s/%s' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid)))

    # Save data
    file = open(
        '../PolarityData/%s/%s/%s.pkl' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid), '{0:04}'.format(simid)),
        'wb')
    pickle.dump(res, file)


def params_array_rand(params, ranges, seeds, nsims):
    """
    Creates a random set of parameter value arrays based on ranges

    :param params: free parameters
    :param ranges: constraints for free parameters
    :param seeds: array of seeds same length as params
    :param nsims: number of parameter value arrays to generate
    :return: valsarray
    """

    # Create empty array for parameter values
    valsarray = np.zeros([nsims, len(params)])

    # Set seeds
    if seeds is None:
        seeds = np.array(range(len(params)))
    seeds = seeds * random.random()

    # Fill array
    for param in range(len(params)):
        for simid in range(nsims):
            random.seed(simid * seeds[param])
            valsarray[simid, param] = random.uniform(ranges[param][0], ranges[param][1])

    return valsarray


def edit_paramset(p, params, valsarray):
    """
    Creates new parameter set based on base parameter set and values to change

    :param p: base parameter set
    :param params: free parameters
    :param valsarray: array of values for free parameters
    :return: updated parameter set, p
    """

    if len(params) > 1:
        for param in range(len(params)):
            setattr(p, params[param], valsarray[param])
    else:
        setattr(p, params[0], valsarray)
    return p


def func_parsim(m, p, params, vals, jobid, subjobid, simid, nsims, compression):
    """

    :param m: model
    :param p: parameter set
    :param params: params to vary
    :param vals: set of values for parameters
    :param jobid: jobid
    :param subjobid: subjobid
    :param simid: simid
    :param nsims: number of simulations
    :param compression:
    :return:
    """

    if subjobid < nsims:

        print('1')

        # Set parameters
        p = edit_paramset(p, params, vals)

        print('2')

        # Run model
        res = run_model(m, p)

        print('3')

        # Save data
        savedata(res, jobid, subjobid, simid, compression)

        print('4')


def func_genalg(m, p, params, vals, jobid, subjobid, simid, pop, compression):
    """
    Performs simulation(s) with given parameters, calculates score based on performance

    :param m: model
    :param p: base parameter set
    :param params: free parameters
    :param vals: array of values for free parameters
    :param jobid: job id
    :param subjobid: subjob id
    :param pop: genetic algorithm population size
    :param compression: 1 = save last time point only, 0 = save all time points
    :return: score representing success of simulation(s)
    """

    if subjobid < pop:
        # Set parameters
        p = edit_paramset(p, params, vals)

        # Run model
        res = run_model(m, p)

        # Assess performance
        [data, labels] = stats(res)
        score = abs(data.asi[-1])

        # Save data
        res.score = score
        savedata(res, jobid, subjobid, simid, compression=compression)

    return score


def func_genalg_180222(m, p, params, vals, jobid, subjobid, simid, pop, compression):
    """
    Performs simulation(s) with given parameters, calculates score based on performance

    :param m: model
    :param p: base parameter set
    :param params: free parameters
    :param vals: array of values for free parameters
    :param jobid: job id
    :param subjobid: subjob id
    :param pop: genetic algorithm population size
    :param compression: results compression
    :return: score representing success of simulation(s)
    """

    if subjobid < pop:
        # Set parameters
        p = edit_paramset(p, params, vals)

        # Run model
        res = run_model(m, p)

        # Assess performance
        base = loaddata(9999, 0, 0)
        mse_a = np.mean(((res.aco - base.aco) ** 2))
        mse_p = np.mean(((res.pco - base.pco) ** 2))
        score = np.mean([mse_a, mse_p])

        # Save data
        res.score = score
        savedata(res, jobid, subjobid, simid, compression=compression)

    return score


def genalg_newgen(valsarray, scores, ranges, seeds):
    """
    Creates next generation of parameter values for genetic algorithm
    Should change to have mutation, crossover and survival rates

    :param valsarray: array of parameter values from previous generaton
    :param scores: performance scores for parameter sets
    :param pop: population size
    :param seeds: for randomisation
    :return: valsarray
    """

    # Set seeds
    rand = random.random()
    if seeds is None:
        seeds = np.array(range(1, len(valsarray[0, :]) + 1)) * rand
    else:
        seeds = seeds * rand

    # Ranking
    print(np.mean(scores))
    valsarraytop = valsarray[np.array(scores).argsort()[:int(len(valsarray[:, 0]) * 0.2)], :]

    # Create empty array for new parameter values
    valsarray[:, :] = 0

    # Generate next generation of parameter values
    for param in range(len(valsarray[0, :])):

        # Mixing
        np.random.seed(int(seeds[param] * 1000))
        valsarray[:, param] = np.random.choice(valsarraytop[:, param], len(valsarray[:, 0]))

        # Mutation
        for simid in range(len(valsarray[:, 0])):
            random.seed(simid * seeds[param])
            if random.uniform(0, 1) < 0.05:
                valsarray[simid, param] = random.uniform(ranges[param][0], ranges[param][1])

    return valsarray


############################ ALGORITHMS ##########################


def alg_singlesim(m, p, jobid=0, subjobid=0, simid=0, compression=1):
    """

    :param m: model
    :param p: parameter set
    :param jobid: job id
    :param subjobid: subjob id
    :param simid: simulation id
    :param compression: specify compression
    :return: saves results (simid.res) in directory jobid/subjobid
    """

    # Run model
    res = run_model(m, p)

    # Save data
    savedata(res, jobid, subjobid, simid, compression)


def alg_parsim(m, p, params, vals, jobid=0, subjobid=0, cores=multiprocessing.cpu_count(), compression=1):
    """

    :param m: model
    :param p: parameter set
    :param params: parameters to vary
    :param vals: set of values for parameters
    :param jobid: jobid
    :param subjobid: subjobid
    :param cores: cores for parallel programming
    :param compression: compression of saved data
    """

    # Create array of parameter values
    if len(params) > 1:
        valsarray = list(itertools.product(*vals))
        nsims = len(valsarray)
    else:
        valsarray = vals[0]
        nsims = len(vals[0])

    # Run simulations
    Parallel(n_jobs=cores, verbose=51)(
        delayed(func_parsim)(m, p, params, valsarray[simid], jobid, subjobid, simid, nsims, compression) for simid
        in range(nsims))


def alg_parsim_clust(m, p, params, vals, jobid, subjobid, cores, node, compression=1):
    """

    :param m: model
    :param p: parameter set
    :param params: parameters to vary
    :param vals: set of values for parameters
    :param jobid: jobid
    :param subjobid: subjobid
    :param cores: cores for parallel programming
    :param node: node id (e.g. int(sys.argv[1]))
    :param compression: compression of saved data
    """


    # Allocate subjobs to node
    subjobmin = node * cores
    subjobmax = (node + 1) * cores

    # Create array of parameter values
    if len(params) > 1:
        valsarray = list(itertools.product(*vals))
        nsubjobs = len(valsarray)
    else:
        valsarray = vals[0]
        nsubjobs = len(vals[0])

    # Run simulations
    Parallel(n_jobs=cores)(
        delayed(func_parsim)(m, p, params, valsarray[simid], jobid, subjobid, simid, nsubjobs, compression) for simid
        in range(subjobmin, subjobmax))


def alg_parsim_rand(m, p, params, ranges, nsims, jobid, subjobid, cores=multiprocessing.cpu_count(), seeds=None,
                    compression=1):
    """

    :param m: model
    :param p: parameter set
    :param params: parameters to vary
    :param ranges: constraints for free parameters
    :param nsims: number of simulations to perform
    :param jobid: jobid
    :param subjobid: subjobid
    :param cores: cores for parallel programming
    :param seeds: array of seeds same length as params
    :param compression: compression of saved data
    """

    # Set seeds
    if seeds is None:
        seeds = range(len(params))

    # Create array of parameter values
    valsarray = params_array_rand(params, ranges, seeds, nsims)

    # Run simulations
    Parallel(n_jobs=cores, verbose=51)(
        delayed(func_parsim)(m, p, params, valsarray[simid], jobid, subjobid, simid, nsims, compression) for simid in
        range(nsims))


def alg_parsim_rand_clust(m, p, params, ranges, nsims, jobid, subjobid, cores, node, seeds=None, compression=1):
    """

    :param m: model
    :param p: parameter set
    :param params: parameters to vary
    :param ranges: constraints for free parameters
    :param nsims: number of simulations to perform
    :param jobid: jobid
    :param subjobid: subjobid
    :param cores: cores for parallel programming
    :param node: node id (e.g. int(sys.argv[1]))
    :param seeds: array of seeds same length as params
    :param compression: compression of saved data
    """

    # Set seeds
    if seeds is None:
        seeds = range(len(params))

    # Create array of parameter values
    valsarray = params_array_rand(params, ranges, seeds, nsims)

    # Allocate subjobs to node
    subjobmin = node * cores
    subjobmax = (node + 1) * cores

    # Run simulations
    Parallel(n_jobs=cores)(
        delayed(func_parsim)(m, p, params, valsarray[simid], jobid, subjobid, simid, nsims, compression) for simid in
        range(subjobmin, subjobmax))


def gen_alg(m, p, func, params, ranges, pop, gens, jobid, cores=multiprocessing.cpu_count(), seeds=None, compression=1):
    """
    Genetic algorithm

    :param m: model
    :param p: base parameter set
    :param func: function to generate success score
    :param params: free parameters
    :param ranges: constraints for free parameters
    :param pop: genetic algorithm population size
    :param gens: genetic algorithm generations
    :param jobid: job id
    :param cores: cores to use
    :param seeds: optional, to fix certain parameters together if needed
    :param compression: results compression
    return: none
    """

    # Create array of parameter values
    valsarray = params_array_rand(params, ranges, seeds, pop)

    # Optimisation
    for gen in range(gens):
        # Simulations
        scores = Parallel(n_jobs=cores, verbose=51)(
            delayed(func)(m=m, p=p, params=params, vals=valsarray[simid], jobid=jobid, subjobid=gen, simid=simid,
                          pop=pop,
                          compression=compression) for simid in range(pop))

        # Next generation
        valsarray = genalg_newgen(valsarray, scores, ranges, seeds)

    # Print results
    print(np.array(range(pop))[[np.array(scores).argsort()]], np.array(scores)[[np.array(scores).argsort()]])


def gen_alg_clust(m, p, func, params, ranges, jobid, cores, nodes, node, innergens=1, seeds=None, compression=1):
    """
    Genetic algorithm to run on cluster

    :param m: model
    :param p: base parameter set
    :param func: function to generate success score
    :param params: free parameters
    :param ranges: constraints for free parameters
    :param innergens: generations of algorithm between pooling events
    :param jobid: job id
    :param cores: cores per node (=population for that node). NB pop = cores * nodes
    :param nodes: number of parallel nodes (e.g. int(sys.argv[2])). Set by array in run.sh
    :param node: node id (e.g. int(sys.argv[1]))
    :param seeds:
    :param compression: results compression
    :return: none
    """

    # Import/create parameter set
    gen = countsubjobs(jobid)
    if gen != 0:
        valsarray = np.zeros([cores * nodes, len(params)])
        scores = np.zeros([cores * nodes])
        for simid in range(cores * nodes):
            res = loaddata(jobid, gen - 1, simid)
            scores[simid] = res.score
            for param in range(len(params)):
                valsarray[simid, param] = getattr(res.p, params[param])
        valsarray = genalg_newgen(valsarray, scores, ranges, seeds)

    else:
        valsarray = params_array_rand(params, ranges, seeds, cores)

    # Optimisation
    for innergen in range(innergens):
        # Simulations
        scores = Parallel(n_jobs=cores, verbose=0)(
            delayed(func)(m=m, p=p, params=params, vals=valsarray[sim], jobid=jobid, subjobid=gen,
                          simid=(node * cores) + sim, pop=cores, compression=compression) for sim in range(cores))

        # Next generation
        valsarray = genalg_newgen(valsarray, scores, ranges, seeds)


############################## ANALYSIS ##########################


def loaddata(jobid, subjobid, simid):
    data = open(
        '../PolarityData/%s/%s/%s.pkl' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid), '{0:04}'.format(simid)),
        'rb')
    res = pickle.load(data)
    return res


def countsims(jobid, subjobid):
    count = 0
    for root, dirs, files in os.walk('../PolarityData/%s/%s' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid))):
        for file in files:
            if file.endswith('.pkl'):
                count += 1
    return count


def countsubjobs(jobid):
    count = 0
    for root, dirs, files in os.walk('../PolarityData/%s' % ('{0:04}'.format(jobid))):
        for dirs in dirs:
            count += 1
    return count


def stats_batch(job, subjobs):
    nsubjobs = len(subjobs)

    res = loaddata(jobid=job, subjobid=0)
    tsteps = int(res.Tmax / res.deltat) + 1

    class data:
        asi = np.zeros([nsubjobs, tsteps])
        a_cyt = np.zeros([nsubjobs, tsteps])
        a_mem = np.zeros([nsubjobs, tsteps])
        a_mem_cyt = np.zeros([nsubjobs, tsteps])
        a_size = np.zeros([nsubjobs, tsteps])
        p_cyt = np.zeros([nsubjobs, tsteps])
        p_mem = np.zeros([nsubjobs, tsteps])
        p_mem_cyt = np.zeros([nsubjobs, tsteps])
        p_size = np.zeros([nsubjobs, tsteps])
        subjob = np.zeros([nsubjobs, tsteps])


    class labels:
        asi = 'Asymmetry index'
        a_cyt = 'A cytoplasmic concentration [μm⁻³]'
        a_mem = 'A domain concentration [μm⁻²]'
        a_mem_cyt = 'A membrane:cytoplasmic ratio'
        a_size = 'A domain size [μm]'
        p_cyt = 'P cytoplasmic concentration [μm⁻³]'
        p_mem = 'P domain concentration [μm⁻²]'
        p_mem_cyt = 'P membrane;cytoplasmic ratio'
        p_size = 'P domain size [μm]'

    for subjob in range(nsubjobs):
        res = loaddata(jobid=job, subjobid=subjobs[subjob])
        data.asi[subjob, :] = (2 * np.sum((np.sign(res.aco - res.pco) + 1) / 2, axis=1) - res.p.xsteps) / res.p.xsteps
        data.a_mem[subjob, :] = np.amax(res.aco, axis=1)
        data.a_cyt[subjob, :] = (res.p.pA - res.p.psi * np.mean(res.aco, axis=1))
        data.a_mem_cyt[subjob, :] = data.a_mem[subjob] / data.a_cyt[subjob]
        data.a_size[subjob, :] = np.sum(res.aco.transpose() > (0.5 * np.tile(data.a_mem[subjob], [res.p.xsteps, 1])),
                                        axis=0) * res.p.L / res.p.xsteps
        data.p_mem[subjob, :] = np.amax(res.pco, axis=1)
        data.p_cyt[subjob, :] = (res.p.pA - res.p.psi * np.mean(res.aco, axis=1))
        data.p_mem_cyt[subjob, :] = data.p_mem[subjob] / data.p_cyt[subjob]
        data.p_size[subjob, :] = np.sum(res.pco.transpose() > (0.5 * np.tile(data.p_mem[subjob], [res.p.xsteps, 1])),
                                        axis=0) * res.p.L / res.p.xsteps
        data.subjob[subjob, :] = subjob

    return data, labels


def stats(res):
    tsteps = int(res.p.Tmax / res.p.deltat) + 1

    class data:
        asi = np.zeros([tsteps])
        a_cyt = np.zeros([tsteps])
        a_mem = np.zeros([tsteps])
        a_mem_cyt = np.zeros([tsteps])
        a_size = np.zeros([tsteps])
        p_cyt = np.zeros([tsteps])
        p_mem = np.zeros([tsteps])
        p_mem_cyt = np.zeros([tsteps])
        p_size = np.zeros([tsteps])
        subjob = np.zeros([tsteps])

    class labels:
        asi = 'Asymmetry index'
        a_cyt = 'A cytoplasmic concentration [μm⁻³]'
        a_mem = 'A domain concentration [μm⁻²]'
        a_mem_cyt = 'A membrane:cytoplasmic ratio'
        a_size = 'A domain size [μm]'
        p_cyt = 'P cytoplasmic concentration [μm⁻³]'
        p_mem = 'P domain concentration [μm⁻²]'
        p_mem_cyt = 'P membrane;cytoplasmic ratio'
        p_size = 'P domain size [μm]'

    data.asi = (2 * np.sum((np.sign(res.aco - res.pco) + 1) / 2, axis=1) - res.p.xsteps) / res.p.xsteps
    data.a_mem = np.amax(res.aco, axis=1)
    data.a_cyt = (res.p.pA - res.p.psi * np.mean(res.aco, axis=1))
    data.a_mem_cyt = data.a_mem / data.a_cyt
    data.a_size = np.sum(res.aco.transpose() > (0.5 * np.tile(data.a_mem, [res.p.xsteps, 1])),
                         axis=0) * res.p.L / res.p.xsteps
    data.p_mem = np.amax(res.pco, axis=1)
    data.p_cyt = (res.p.pA - res.p.psi * np.mean(res.aco, axis=1))
    data.p_mem_cyt = data.p_mem / data.p_cyt
    data.p_size = np.sum(res.pco.transpose() > (0.5 * np.tile(data.p_mem, [res.p.xsteps, 1])),
                         axis=0) * res.p.L / res.p.xsteps

    return data, labels


def asirank(job, subjobs):
    # Print ranked list of best subjobs along with their asi
    [data, labels] = stats_batch(job, subjobs)
    print(data.subjob[abs(data.asi[:, -1]).argsort()], data.asi[abs(data.asi[:, -1]).argsort()])


############################## PLOTS ##############################


def parplot(i, ax, aco, pco, p):
    ax.plot(np.array(range(p.xsteps)) / (p.xsteps - 1) * p.L, aco[i, :], label='A', c='r')
    ax.plot(np.array(range(p.xsteps)) / (p.xsteps - 1) * p.L, pco[i, :], label='P', c='dodgerblue')
    ax.set_xlabel('x [μm]')
    ax.set_ylabel('Concentration [a.u.]')
    concmax = max(aco.max(), pco.max())
    ax.set_ylim(0, 1.1 * concmax)
    # ax.legend()


def plot_singlesim(jobid=0, subjobid=0, simid=0):
    res = loaddata(jobid, subjobid, simid)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    parplot(-1, ax, res.aco, res.pco, res.p)
    plt.show()


def sliderplot(jobid=0, subjobid=0, simid=0):
    res = loaddata(jobid, subjobid, simid)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Time (s)', 0, res.p.Tmax, valinit=0, valfmt='%d')

    def update_slider(val):
        ax.clear()
        i = int(sframe.val / res.p.deltat)
        parplot(i, ax, res.aco, res.pco, res.p)
        ax.set_title('Time (s): {0:.0f}'.format(i * res.p.deltat))

    parplot(0, ax, res.aco, res.pco, res.p)
    ax.set_title('Time (s): {0:.0f}'.format(0 * res.p.deltat))
    sframe.on_changed(update_slider)
    plt.show()


def anim(jobid=0, subjobid=0, simid=0, animrate=100, framerate=24):
    """
    :param animrate: seconds of model time per second of animation
    :param framerate: frames per second

    """

    res = loaddata(jobid, subjobid, simid)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    frames = range(0, int(res.p.Tmax / res.p.deltat), int(animrate / (framerate * res.p.deltat)))

    def update_anim(i):
        ax.clear()
        parplot(i, ax, res.aco, res.pco, res.p)
        ax.set_title('Time (s): {0:.0f}'.format(i * res.p.deltat))

    paranim = animation.FuncAnimation(fig, update_anim, frames=iter(frames), save_count=len(frames))
    writer = animation.writers['ffmpeg']
    writer = writer(fps=framerate, bitrate=2000)
    paranim.save('animation.mp4', writer=writer)


def paramsliderplot(jobid, subjobid, param):
    sims = countsims(jobid, subjobid)
    xsteps = loaddata(jobid=jobid, subjobid=subjobid, simid=0).p.xsteps
    params = np.zeros(sims)
    aco = np.zeros([sims, xsteps])
    pco = np.zeros([sims, xsteps])

    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        params[simid] = getattr(res.p, param)
        aco[simid, :] = res.aco[-1, :]
        pco[simid, :] = res.pco[-1, :]

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, param, 0, max(params), valinit=0)

    def update_slider(val):
        ax.clear()
        val = sframe.val
        subjob = (np.abs(params - val)).argmin()
        parplot(subjob, ax, aco, pco, res.p)

    parplot(0, ax, aco, pco, res.p)
    sframe.on_changed(update_slider)
    plt.show()


def paramsliderplot2(jobid, subjobid, param1, param2):
    sims = countsims(jobid, subjobid)
    xsteps = loaddata(jobid, subjobid, 0).p.xsteps
    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)
    aco = np.zeros([sims, xsteps])
    pco = np.zeros([sims, xsteps])

    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)
        aco[simid, :] = res.aco[-1, :]
        pco[simid, :] = res.pco[-1, :]

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    param1_slider_ax = plt.axes([0.25, 0.1, 0.65, 0.03])
    param1_slider = Slider(param1_slider_ax, param1, 0, max(param1vals), valinit=0)
    param2_slider_ax = plt.axes([0.25, 0.05, 0.65, 0.03])
    param2_slider = Slider(param2_slider_ax, param2, 0, max(param2vals), valinit=0)

    def update_slider(val):
        ax.clear()
        param1_val = param1_slider.val
        param2_val = param2_slider.val
        subjob = (np.abs(param1vals - param1_val) * np.abs(param2vals - param2_val)).argmin()
        parplot(subjob, ax, aco, pco, res.p)

    parplot(0, ax, aco, pco, res.p)
    param1_slider.on_changed(update_slider)
    param2_slider.on_changed(update_slider)
    plt.show()


# def paramsliderplot_multi(job, subjobs, params):
#     # Should count the number of subjobs
#
#     xsteps = loaddata(job=job, subjob=0).xsteps
#
#     paramsvals = np.zeros([subjobs, len(params)])
#     aco = np.zeros([subjobs, xsteps])
#     pco = np.zeros([subjobs, xsteps])
#
#     for subjob in range(subjobs):
#         p = loaddata(job=job, subjob=subjob)
#         for param in range(len(params)):
#             paramsvals[subjob, param] = getattr(p, params[param])
#         aco[subjob, :] = p.aco[-1, :]
#         pco[subjob, :] = p.pco[-1, :]
#
#     plt.clf()
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     plt.subplots_adjust(left=0.25, bottom=0.25)
#
#     for param in range(len(params)):
#         param_slider_ax = plt.axes([0.25, 0.1, 0.65, 0.03])
#         param_slider = Slider(param_slider_ax, params[param], 0, max(paramsvals[:, param]), valinit=0)
#
#
#     def update_slider(val):
#         ax.clear()
#         param1_val = param1_slider.val
#         param2_val = param2_slider.val
#         subjob = (np.abs(param1vals - param1_val) * np.abs(param2vals - param2_val)).argmin()
#         parplot(subjob, ax, aco, pco, p)
#
#     parplot(0, ax, aco, pco, p)
#     param1_slider.on_changed(update_slider)
#     param2_slider.on_changed(update_slider)
#     plt.show()


def dataplot(jobid, subjobid, param, var):
    sims = countsims(jobid, subjobid)
    (d, labels) = stats_batch(jobid, range(sims))
    paramvals = np.zeros(sims)

    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        paramvals[simid] = getattr(res.p, param)

    plt.plot(paramvals, getattr(d, var)[:, -1])
    plt.ylabel(getattr(labels, var))
    plt.xlabel(param)
    plt.show()


def dataplot2(jobid, subjobid, param1, param2, var):
    sims = countsims(jobid, subjobid)
    (d, labels) = stats_batch(jobid, range(sims))
    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)

    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    vals = np.unique(param2vals)
    for param2val in range(len(vals)):
        print(param2vals[param2val])
        ax.plot(param1vals[param2vals == vals[param2val]],
                getattr(d, var)[:, -1][param2vals == vals[param2val]])
    # cbar = plt.colorbar(ax)
    # cbar.set_label(param2)
    plt.ylabel(getattr(labels, var))
    plt.xlabel(param1)
    plt.show()


def dataplot2_heatmap(jobid, subjobid, param1, param2, var):
    sims = countsims(jobid, subjobid)
    (d, labels) = stats_batch(jobid, range(sims))
    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)
    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)

    sc = plt.scatter(param1vals, param2vals, c=getattr(d, var)[:, -1])
    cbar = plt.colorbar(sc)
    cbar.set_label(getattr(labels, var))
    plt.xlabel(param1)
    plt.ylabel(param2)

    plt.show()


def dataplot2_heatmap_180220(jobid, subjobid, param1, param2, var):
    sims = countsims(jobid, subjobid)
    (d, labels) = stats_batch(jobid, range(sims))

    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)
    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)

    asi = getattr(d, var)[:, -1]

    sc = plt.scatter(param1vals[abs(asi) < 0.25], param2vals[abs(asi) < 0.25])
    # cbar = plt.colorbar(sc)
    # cbar.set_label(getattr(labels, var))
    plt.xlabel(param1)
    plt.ylabel(param2)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.show()


def dataplot2_heatmap_180529(jobid, subjobid, param1, param2):
    sims = countsims(jobid, subjobid)

    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)
    pols = np.zeros(sims)
    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)
        pols[simid] = test_180529(res)

    plt.scatter(param1vals[abs(pols) == 0], param2vals[abs(pols) == 0], c='silver', marker='s')
    plt.scatter(param1vals[abs(pols) == 1], param2vals[abs(pols) == 1], c='dodgerblue', marker='s')
    plt.scatter(param1vals[abs(pols) == 2], param2vals[abs(pols) == 2], c='r', marker='s')

    plt.xlabel(param1)
    plt.ylabel(param2)
    # plt.xlim([0, 1])
    # plt.ylim([0, 1])
    plt.show()


def genalg_plot1(jobid, base):
    nsubjobs = countsubjobs(jobid)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Generation', 0, nsubjobs, valinit=0, valfmt='%d')

    def update_slider(val):
        ax.clear()
        subjobid = int(sframe.val)
        nsims = countsims(jobid, subjobid)

        # Plot base data
        res = loaddata(base[0], base[1], base[2])
        ax.plot(np.array(range(res.p.xsteps)) / (res.p.xsteps - 1) * res.p.L, res.aco[-1, :], c='g')
        ax.plot(np.array(range(res.p.xsteps)) / (res.p.xsteps - 1) * res.p.L, res.pco[-1, :], c='r')
        ax.set_xlabel('Anterior - Posterior [μm]')
        ax.set_ylabel('Concentration [μm⁻²]')
        concmax = max(res.aco.max(), res.pco.max())
        ax.set_ylim(0, 1.5 * concmax)

        # Plot genentic algorithm simulations
        for simid in range(nsims):
            res = loaddata(jobid, subjobid, simid)
            ax.plot(np.array(range(res.p.xsteps)) / (res.p.xsteps - 1) * res.p.L, res.aco[-1, :], c='g', alpha=0.1)
            ax.plot(np.array(range(res.p.xsteps)) / (res.p.xsteps - 1) * res.p.L, res.pco[-1, :], c='r', alpha=0.1)

    sframe.on_changed(update_slider)
    plt.show()


def genalg_spider(jobid, subjobid, params):
    # Set up plot
    angles = [n / float(len(params)) * 2 * np.pi for n in range(len(params))]
    angles += angles[:1]
    ax = plt.subplot(111, polar=True)
    plt.xticks(angles[:-1], params)
    ax.set_rlabel_position(0)

    # Plot data
    nsims = countsims(jobid, subjobid)
    for simid in range(nsims):
        values = np.zeros(len(params) + 1)
        res = loaddata(jobid, subjobid, simid)
        for param in range(len(params)):
            values[param] = getattr(res.p, params[param])
        values[-1] = values[0]
        ax.plot(angles, values)
    plt.show()


def genalg_spider_slider(jobid, params):
    nsubjobs = countsubjobs(jobid)

    plt.clf()
    fig = plt.figure()
    ax = plt.subplot(111, polar=True)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Generation', 0, nsubjobs, valinit=0, valfmt='%d')

    angles = [n / float(len(params)) * 2 * np.pi for n in range(len(params))]
    angles += angles[:1]

    def update_slider(val):
        ax.clear()
        subjobid = int(sframe.val)
        ax.set_xticks(angles[:-1], params)
        ax.set_rlabel_position(0)

        # Plot data
        nsims = countsims(jobid, subjobid)
        for simid in range(nsims):
            values = np.zeros(len(params) + 1)
            res = loaddata(jobid, subjobid, simid)
            for param in range(len(params)):
                values[param] = getattr(res.p, params[param])
            values[-1] = values[0]
            ax.plot(angles, values)

    sframe.on_changed(update_slider)
    plt.show()


########################## STABILITY ANALYSIS #########################


def init_phaseplot():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('A [μm⁻²]')
    ax.set_ylabel('P [μm⁻²]')
    return fig, ax


def trajplot(i, ax, aco, pco, p):
    ax.scatter(aco[i, :], pco[i, :], s=1)
    ax.set_xlabel('aPAR [μm⁻²]')
    ax.set_ylabel('pPAR [μm⁻²]')
    ax.set_title('Time (s): {0:.0f}'.format(i * p.deltat))
    xmax = aco.max()
    ymax = pco.max()
    ax.set_xlim(0, 1.1 * xmax)
    ax.set_ylim(0, 1.1 * ymax)


def trajectories_slider(job, subjob, ax, sframe):
    [aco, pco, p] = loaddata(job, subjob)

    def update_slider(val):
        ax.clear()
        i = int(sframe.val / p.deltat)
        trajplot(i, ax, aco, pco, p)

    trajplot(0, ax, aco, pco, p)
    sframe.on_changed(update_slider)


def alg_stab1(m, p):
    [fig, ax] = init_phaseplot()
    [x1, y1, x2, y2] = nullcline(m, p)
    for a in range(len(x1[0, :])):
        ax.scatter(x1[:, a], y1, c='g', s=5)
        ax.scatter(x2, y2[:, a], c='r', s=5)
    plt.show()


def alg_stab2(m, p, a0=0, p0=0):
    [fig, ax] = init_phaseplot()
    [x1, y1, x2, y2] = nullcline(m, p)
    for a in range(len(x1[0, :])):
        ax.scatter(x1[:, a], y1, c='g', s=5)
        ax.scatter(x2, y2[:, a], c='r', s=5)
    [aco, pco] = trajectory(m, p, a0, p0)
    ax.plot(aco, pco)
    plt.show()


def alg_stab3(m, p):
    [fig, ax] = init_phaseplot()
    [x1, y1, x2, y2] = nullcline(m, p)
    for a in range(len(x1[0, :])):
        ax.scatter(x1[:, a], y1, c='g', s=5)
        ax.scatter(x2, y2[:, a], c='r', s=5)

    # a0eq = np.amax(m.equilibrate_aco(p))
    # [aco1, pco1] = trajectory(m, p, a0eq, 0)
    # ax.plot(aco1, pco1)
    # p0eq = np.amax(m.equilibrate_pco(p))
    # [aco2, pco2] = trajectory(m, p, 0, p0eq)
    # ax.plot(aco2, pco2)
    plt.show()


# def solve_equation_plot(func):
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     [x1, y1, x2, y2] = solve_equation(func)
#     for a in range(len(x1[0, :])):
#         ax.scatter(x1[:, a], y1, c='g', s=5)
#         # ax.scatter(x2, y2[:, a], c='r', s=5)
#     plt.show()
#
#
# def func(x, y):
#     change = x + x * (y ** 2) - 2 * y
#     return change





######################## PARAMETER SETS #######################


gparams0 = MiscParams(L=50, xsteps=500, psi=0.3, Tmax=1000, deltat=0.01, Aeqmin=0, Aeqmax=0.5, Peqmin=0.5, Peqmax=1)

gparams1 = MiscParams(L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.1, Aeqmin=0, Aeqmax=0.5, Peqmin=0.5, Peqmax=1)

paramset0_0 = Model0.Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2,
                            ePneg=1, eAneg=2, pA=1.56, pP=1, pgen=gparams1)

paramset0_1 = Model0.Params(Da=1, Dp=1, konA=1, koffA=0.1, konP=1, koffP=0.1, kAP=1, kPA=1,
                            ePneg=2, eAneg=2, pA=1, pP=1, pgen=gparams0)

paramset0_2 = Model0.Params(Da=0.1, Dp=0.1, konA=0.006, koffA=0.005, konP=0.006, koffP=0.005, kAP=1, kPA=1,
                            ePneg=2, eAneg=2, pA=0, pP=1, pgen=gparams0)

paramset4_1 = Model4.Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0, kPA=0.8,
                            kPneg1=1.5, kPneg2=0, kAneg1=0.5, kAneg2=3.5, ePneg=10, eAneg=10, kApos=1, kPpos=0.5,
                            eApos=1, ePpos=10, pA=1.56, pP=1, pgen=gparams1)

paramset4_0 = Model4.Params(Da=1, Dp=1, konA=1, koffA=0.3, konP=1, koffP=0.3, kAP=3, kPA=3,
                            kPneg1=1, kPneg2=1, kAneg1=1, kAneg2=1, ePneg=20, eAneg=20, kApos=1, kPpos=1, eApos=20,
                            ePpos=20, pA=1, pP=1, pgen=gparams0)

paramset2_1 = Model2.Params(Da=0.28, Dp=0.15, konA1=0.0085, koffA=0.0054, konP1=0.0474, koffP=0.0073, kAP=0, kPA=0.15,
                            kAneg=0.3, kPneg=0.3, ePneg=1, eAneg=10, konA2=0, konP2=1, kApos=0, kPpos=2, eApos=0,
                            ePpos=10, pA=1.56, pP=1, pgen=gparams1)

paramset2_2 = Model2.Params(Da=0.28, Dp=0.15, konA1=0.0085, koffA=0.0054, konP1=0.02539162, koffP=0.0073, kAP=0,
                            kPA=0.15, kAneg=0.3, kPneg=0.3, ePneg=1, eAneg=10, konA2=0, konP2=0.00257667, kApos=0,
                            kPpos=1, eApos=0, ePpos=10, pA=1.56, pP=1, pgen=gparams1)


#######################################################################

# To do:

# Mutation and crossover rates for genalg
# Genetic algorithm optimisation
# A way to plot the results from genetic algorithm e.g. spider plot or PC plot
# Fix empty/double plot issue
# Function that finds position of stable fixed points
# Test models 2 and 3 with 'realistic' parameters. Use genetic algorithm
# Stochastic model
# Paramsliderplot with arbitrary number of sliders
# Change parplot function
# Paramsliderplots don't like the ends of the sliders
# Better way to deal with axis labels?
# Better definition of asymmetry index
# Maybe autoname jobs (with more characters) to avoid overwriting data. Can be manually overwritten
# Animations not starting from zero
# Make functions compatible with 4 component model
# More comments to clarify things
# Should autosave a file in job directory indicating the type of simulation and the arguments
# More shared code between plot functions
# Change jobid, subjodid, simid arguments to single directory argument?
# Function to make animation from any slider plot???


#########################################################################
