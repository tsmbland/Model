import numpy as np
from joblib import Parallel, delayed
import os
import pickle
import itertools
import random
import multiprocessing
import sys
import copy
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import animation
import seaborn as sns
import glob
import pandas as pd

sns.set()
sns.set_style("ticks")

#########################  ALGORITHM FUNCTIONS ##########################


datadirec = None

"""
From local to server: '../../../../../../../Volumes/lab-goehringn/working/Tom/ModelData'

From server to server: '../working/Tom/ModelData'

From local to local: '../../ModelData'
"""


class Compression2:
    def __init__(self, res):
        self.scores = res.scores
        self.params = res.params


def subjobslist(jobid):
    dlist = glob.glob('%s/%s/*/' % (datadirec, '{0:04}'.format(jobid)))
    dlist = [os.path.basename(os.path.normpath(x))[:-4] for x in dlist if '!' not in x]
    return dlist


def simidlist(jobid, subjobid):
    slist = glob.glob('%s/%s/%s/*.pkl' % (datadirec, '{0:04}'.format(jobid), '{0:04}'.format(subjobid)))
    slist = [int(os.path.basename(os.path.normpath(x))[:-4]) for x in slist if '!' not in x]
    return slist


def savedata(res, jobid, subjobid, simid, compression):
    """

    :param res: results from simulation
    :param jobid:
    :param subjobid:
    :param simid
    :return:
    """

    # Create job directory
    if not os.path.isdir('%s/%s' % (datadirec, '{0:04}'.format(jobid))):
        os.makedirs('%s/%s' % (datadirec, '{0:04}'.format(jobid)))

    # Create subjob directory
    if not os.path.isdir('%s/%s/%s' % (datadirec, '{0:04}'.format(jobid), '{0:04}'.format(subjobid))):
        os.makedirs('%s/%s/%s' % (datadirec, '{0:04}'.format(jobid), '{0:04}'.format(subjobid)))

    # Create pickle file
    file = open(
        '%s/%s/%s/%s.pkl' % (datadirec, '{0:04}'.format(jobid), '{0:04}'.format(subjobid), '{0:04}'.format(simid)),
        'wb')

    # Compress
    if compression == 0:
        pickle.dump(res, file)
    if compression == 1:
        res.compress()
        pickle.dump(res, file)
    elif compression == 2:
        pickle.dump(Compression2(res), file)


def loaddata(jobid, subjobid, simid):
    data = open(
        '%s/%s/%s/%s.pkl' % (datadirec, '{0:04}'.format(jobid), '{0:04}'.format(subjobid), '{0:04}'.format(simid)),
        'rb')
    res = pickle.load(data)
    return res


def params_array_rand(params, ranges, seeds, nsims):
    """
    Creates a random set of parameter value arrays based on ranges

    :param params: free parameters (list of names)
    :param ranges: constraints for free parameters
    :param seeds: array of seeds same length as params
    :param nsims: number of parameter value arrays to generate
    :return: valsarray
    """

    # Create empty array for parameter values
    valsarray = np.zeros([nsims, len(params)])

    # Set seeds
    if seeds is None:
        seeds = np.random.rand(len(params))
    else:
        seeds *= random.random()

    # Fill array
    for param in range(len(params)):
        for simid in range(nsims):
            random.seed((1 + simid) * seeds[param])
            valsarray[simid, param] = random.uniform(ranges[param][0], ranges[param][1])

    return valsarray


def edit_paramset(p, params, valsarray):
    """
    Creates new parameter set based on base parameter set and values to change

    :param p: base parameter set (class)
    :param params: free parameters (list of names)
    :param valsarray: array of values for free parameters
    :return: updated parameter set, p
    """

    if len(params) > 1:
        for param in range(len(params)):
            setattr(p, params[param], valsarray[param])
    else:
        setattr(p, params[0], valsarray)
    return p


def func_parsim(m, params, vals, jobid, subjobid, simid, nsims, compression, funcs=[]):
    """

    :param m: model
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
        # Initiate new model
        m = copy.deepcopy(m)

        # Set parameters
        p = edit_paramset(m.params, params, vals)
        m.params = p

        # Run model
        alg_singlesim(m, jobid, subjobid, simid, compression, funcs)


def func_genalg_0000(m, params, vals, jobid, subjobid, simid, pop):
    """
    Performs simulation(s) with given parameters, calculates score based on performance

    :param m: model
    :param params: free parameters
    :param vals: array of values for free parameters
    :param jobid: job id
    :param subjobid: subjob id
    :param pop: genetic algorithm population size
    :return: score representing success of simulation(s)
    """

    if subjobid < pop:
        # Initiate new model
        m = copy.deepcopy(m)

        # Set parameters
        p = edit_paramset(m.params, params, vals)
        m.params = p

        # Run model
        r = m.run()

        # Calculate score
        asi1 = (np.mean(r.aco[0, 0:len(r.aco[0]) // 2]) - np.mean(r.aco[0, len(r.aco[0]) // 2:])) / (
            2 * (np.mean(r.aco[0, 0:len(r.aco[0]) // 2]) + np.mean(r.aco[0, len(r.aco[0]) // 2:])))
        asi2 = (np.mean(r.pco[0, 0:len(r.pco[0]) // 2]) - np.mean(r.pco[0, len(r.pco[0]) // 2:])) / (
            2 * (np.mean(r.pco[0, 0:len(r.pco[0]) // 2]) + np.mean(r.pco[0, len(r.pco[0]) // 2:])))
        score = 0.5 - (abs(asi1) + abs(asi2)) / 2  # low = polarised

        r.score = score
        savedata(r, jobid, subjobid, simid, 1)

        return score


def func_genalg_0001(m, params, vals, jobid, subjobid, simid, pop):
    """
    Performs simulation(s) with given parameters, calculates score based on performance

    :param m: model
    :param p: base parameter set
    :param params: free parameters
    :param vals: array of values for free parameters
    :param jobid: job id
    :param subjobid: subjob id
    :param pop: genetic algorithm population size
    :return: score representing success of simulation(s)
    """

    if subjobid < pop:
        # Initiate new model
        m = copy.deepcopy(m)

        # Set parameters
        p = edit_paramset(m.params, params, vals)
        m.params = p

        # Run model
        r = m.run()

        # Assess performance
        base = loaddata(9999, 0, 0)
        mse_a = np.mean(((r.aco - base.aco) ** 2))
        mse_p = np.mean(((r.pco - base.pco) ** 2))
        score = np.mean([mse_a, mse_p])

        r.score = score
        savedata(r, jobid, subjobid, simid, 1)

        return score


def genalg_newgen(valsarray, scores, ranges, seeds):
    """
    Creates next generation of parameter values for genetic algorithm
    Should change to have mutation, crossover and survival rates

    :param valsarray: array of parameter values from previous generaton
    :param scores: performance scores for parameter sets (lower is better?)
    :param ranges:
    :param seeds: for randomisation
    :return: valsarray
    """

    # Set seeds
    if seeds is None:
        seeds = np.array(range(1, len(valsarray[0, :]) + 1)) * random.random()
    else:
        seeds = seeds * random.random()

    # Ranking
    print(np.mean(scores))
    valsarraytop = valsarray[np.array(scores).argsort()[:int(len(valsarray[:, 0]) * 0.2)], :]

    # Create empty array for new parameter values
    valsarray[:, :] = 0

    # Generate next generation of parameter values
    for param in range(len(valsarray[0, :])):

        # Mixing
        np.random.seed(int(seeds[param]))
        valsarray[:, param] = np.random.choice(valsarraytop[:, param], len(valsarray[:, 0]))

        # Mutation
        for simid in range(len(valsarray[:, 0])):
            random.seed((1 + simid) * seeds[param])
            if random.uniform(0, 1) < 0.05:
                valsarray[simid, param] = random.uniform(ranges[param][0], ranges[param][1])

    return valsarray


############################ ALGORITHMS ##########################


def alg_singlesim(m, jobid=0, subjobid=0, simid=0, compression=1, funcs=[]):
    """
    Runs a single simulation from a model, and saves that data

    :param m: model
    :param jobid: job id
    :param subjobid: subjob id
    :param simid: simulation id
    :param compression: specify compression
    :param args: functions to run on r. func(r)
    :return: saves results (simid.res) in directory jobid/subjobid
    """

    # Run model
    r = m.run()

    # Run additional functions
    for func in funcs:
        func(r)

    # Compress and save
    savedata(r, jobid, subjobid, simid, compression)


def alg_parsim(m, params, vals, jobid=0, subjobid=0, cores=multiprocessing.cpu_count(), compression=1, funcs=[]):
    """
    Runs multiple simulations in parallel, saves the data
    Specify free parameters and list of values for each, will simulate all combinations

    :param m: model
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
        delayed(func_parsim)(m, params, valsarray[simid], jobid, subjobid, simid, nsims, compression, funcs) for simid
        in range(nsims))


def alg_parsim_clust(m, params, vals, jobid, subjobid, cores, node, compression=2, funcs=[]):
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
        delayed(func_parsim)(m, params, valsarray[simid], jobid, subjobid, simid, nsubjobs, compression, funcs) for
        simid
        in range(subjobmin, subjobmax))


def alg_parsim_rand(m, params, ranges, nsims, jobid=0, subjobid=0, cores=multiprocessing.cpu_count(), seeds=None,
                    compression=1, funcs=[]):
    """


    :param m: model
    :param params: parameters to vary
    :param ranges: constraints for free parameters
    :param nsims: number of simulations to perform
    :param jobid: jobid
    :param subjobid: subjobid
    :param cores: cores for parallel programming
    :param seeds: array of seeds same length as params
    :param compression: compression of saved data
    """

    # Create array of parameter values
    valsarray = params_array_rand(params, ranges, seeds, nsims)

    # Run simulations
    Parallel(n_jobs=cores, verbose=51)(
        delayed(func_parsim)(m, params, valsarray[simid], jobid, subjobid, simid, nsims, compression, funcs) for simid
        in
        range(nsims))


def alg_parsim_rand_clust(m, params, ranges, nsims, jobid, subjobid, cores, node, seeds=None, compression=2, funcs=[]):
    """

    :param m: model
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

    # Create array of parameter values
    valsarray = params_array_rand(params, ranges, seeds, nsims)  # <- wasteful, creating full array on each core

    # Allocate subjobs to node
    subjobmin = node * cores
    subjobmax = (node + 1) * cores

    # Run simulations
    Parallel(n_jobs=cores)(
        delayed(func_parsim)(m, params, valsarray[simid], jobid, subjobid, simid, nsims, compression, funcs) for simid
        in
        range(subjobmin, subjobmax))


def gen_alg(m, func, params, ranges, pop, gens, jobid, cores=multiprocessing.cpu_count(), seeds=None):
    """
    Genetic algorithm

    :param m: model
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

        print(valsarray)
        scores = Parallel(n_jobs=cores, verbose=51)(
            delayed(func)(m, params, valsarray[simid], jobid, gen, simid, pop) for simid in
            range(pop))

        # Next generation
        valsarray = genalg_newgen(valsarray, scores, ranges, seeds)

    # Print results
    print(np.array(range(pop))[[np.array(scores).argsort()]], np.array(scores)[[np.array(scores).argsort()]])


def gen_alg_clust(m, func, params, ranges, jobid, cores, nodes, node, innergens=1, seeds=None):
    """
    Genetic algorithm to run on cluster

    :param m: model
    :param func: function to generate success score
    :param params: free parameters
    :param ranges: constraints for free parameters
    :param innergens: generations of algorithm between pooling events
    :param jobid: job id
    :param cores: cores per node (=population for that node). NB pop = cores * nodes
    :param nodes: number of parallel nodes (e.g. int(sys.argv[2])). Set by array in run.sh
    :param node: node id (e.g. int(sys.argv[1]))
    :param seeds:
    :return: none
    """

    # Import/create parameter set
    gen = len(subjobslist(jobid))
    if gen != 0:
        valsarray = np.zeros([cores * nodes, len(params)])
        scores = np.zeros([cores * nodes])
        for simid in range(cores * nodes):
            res = loaddata(jobid, gen - 1, simid)
            scores[simid] = res.score
            for param in range(len(params)):
                valsarray[simid, param] = getattr(res.params, params[param])
        valsarray = genalg_newgen(valsarray, scores, ranges, seeds)

    else:
        valsarray = params_array_rand(params, ranges, seeds, cores)

    # Optimisation
    for innergen in range(innergens):
        # Simulations
        scores = Parallel(n_jobs=cores, verbose=0)(
            delayed(func)(m, params, valsarray[sim], jobid, gen, (node * cores) + sim, cores) for sim
            in range(cores))

        # Next generation
        valsarray = genalg_newgen(valsarray, scores, ranges, seeds)


########################## ANALYSIS FUNCTIONS ##########################

def mse(res):
    base = loaddata(9999, 0, 0)
    mse_a = np.mean(((res.aco[-1, :] - base.aco[-1, :]) ** 2))
    mse_p = np.mean(((res.pco[-1, :] - base.pco[-1, :]) ** 2))
    score = np.mean([mse_a, mse_p])
    res.scores['mse'] = score


def asi_a(res):
    ant = np.mean(res.aco[-1, 0:len(res.aco[-1, :]) // 2])
    post = np.mean(res.aco[-1, len(res.aco[-1, :]) // 2:])
    asi = (ant - post) / (2 * (ant + post))
    res.scores['asi_a'] = asi


def asi_p(res):
    ant = np.mean(res.pco[-1, 0:len(res.pco[-1, :]) // 2])
    post = np.mean(res.pco[-1, len(res.pco[-1, :]) // 2:])
    asi = (ant - post) / (2 * (ant + post))
    res.scores['asi_p'] = asi


def print_scores_batch(jobid, subjobid):
    for simid in simidlist(jobid, subjobid):
        r = loaddata(jobid, subjobid, simid)
        print([simid, r.scores])


def save_scores_batch(jobid, subjobid):
    df = pd.DataFrame()
    for simid in simidlist(jobid, subjobid):
        try:
            print(simid)
            r = loaddata(jobid, subjobid, simid)
            columns = ['simid']
            vals = [simid]
            for key, value in r.scores.items():
                columns.append(key)
                vals.append(value)
            row = pd.DataFrame([vals], columns=columns)
            df = df.append(row)
        except:
            pass

    df.to_csv('%s/res.csv' % direc_to(jobid, subjobid))


def direc_to(jobid, subjobid):
    return '%s/%s/%s' % (datadirec, '{0:04}'.format(jobid), '{0:04}'.format(subjobid))


############################## ANALYSIS ##########################


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
    sframe = Slider(axframe, 'Time (s)', 0, res.params.Tmax, valinit=0, valfmt='%d')

    def update_slider(val):
        ax.clear()
        i = int(sframe.val / res.params.deltat)
        parplot(i, ax, res.aco, res.pco, res.params)
        ax.set_title('Time (s): {0:.0f}'.format(i * res.params.deltat))

    parplot(0, ax, res.aco, res.pco, res.params)
    ax.set_title('Time (s): {0:.0f}'.format(0 * res.params.deltat))
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
    sims = len(simidlist(jobid, subjobid))
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
    sims = len(simidlist(jobid, subjobid))
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


def dataplot(jobid, subjobid, param, var):
    sims = len(simidlist(jobid, subjobid))
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
    sims = len(simidlist(jobid, subjobid))
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
    sims = len(simidlist(jobid, subjobid))
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
    sims = len(simidlist(jobid, subjobid))
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
    sims = len(simidlist(jobid, subjobid))

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


def genalg_plot1(jobid, base=None):
    nsubjobs = len(subjobslist(jobid))

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Generation', 0, nsubjobs, valinit=0, valfmt='%d')

    def update_slider(val):
        ax.clear()
        subjobid = int(sframe.val)
        nsims = len(simidlist(jobid, subjobid))

        # Plot base data
        if base is not None:
            res = loaddata(base[0], base[1], base[2])
            ax.plot(np.array(range(res.params.xsteps)) / (res.params.xsteps - 1) * res.params.L, res.aco[-1, :], c='g')
            ax.plot(np.array(range(res.params.xsteps)) / (res.params.xsteps - 1) * res.params.L, res.pco[-1, :], c='r')
            concmax = max(res.aco.max(), res.pco.max())
            ax.set_ylim(0, 1.5 * concmax)

        # Plot genentic algorithm simulations
        for simid in range(nsims):
            res = loaddata(jobid, subjobid, simid)
            ax.plot(np.array(range(res.params.xsteps)) / (res.params.xsteps - 1) * res.params.L, res.aco[-1, :], c='g',
                    alpha=0.1)
            ax.plot(np.array(range(res.params.xsteps)) / (res.params.xsteps - 1) * res.params.L, res.pco[-1, :], c='r',
                    alpha=0.1)

        ax.set_xlabel('Anterior - Posterior [μm]')
        ax.set_ylabel('Concentration [μm⁻²]')

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
    nsims = len(simidlist(jobid, subjobid))
    for simid in range(nsims):
        values = np.zeros(len(params) + 1)
        res = loaddata(jobid, subjobid, simid)
        for param in range(len(params)):
            values[param] = getattr(res.p, params[param])
        values[-1] = values[0]
        ax.plot(angles, values)
    plt.show()


def genalg_spider_slider(jobid, params):
    nsubjobs = len(subjobslist(jobid))

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
        nsims = len(simidlist(jobid, subjobid))
        for simid in range(nsims):
            values = np.zeros(len(params) + 1)
            res = loaddata(jobid, subjobid, simid)
            for param in range(len(params)):
                values[param] = getattr(res.params, params[param])
            values[-1] = values[0]
            ax.plot(angles, values)

    sframe.on_changed(update_slider)
    plt.show()
