import numpy as np
from joblib import Parallel, delayed
import os
import pickle
import itertools
import random
import multiprocessing
import sys
import copy

os.chdir(os.path.expanduser('../../ModelData'))


#########################  ALGORITHM FUNCTIONS ##########################


class Compression2:
    def __init__(self, res):
        self.scores = res.scores
        self.params = res.params


def savedata(res, jobid, subjobid, simid, compression):
    """

    :param res: results from simulation
    :param jobid:
    :param subjobid:
    :param simid
    :return:
    """

    # Create job directory
    if not os.path.isdir('%s' % ('{0:04}'.format(jobid))):
        os.makedirs('%s' % ('{0:04}'.format(jobid)))

    # Create subjob directory
    if not os.path.isdir('%s/%s' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid))):
        os.makedirs('%s/%s' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid)))

    # Create pickle file
    file = open(
        '%s/%s/%s.pkl' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid), '{0:04}'.format(simid)), 'wb')

    # Compress
    if compression == 0:
        pickle.dump(res, file)
    if compression == 1:
        res.compress()
        pickle.dump(res, file)
    elif compression == 2:
        pickle.dump(Compression2(res), file)

    # Save data
    file = open(
        '%s/%s/%s.pkl' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid), '{0:04}'.format(simid)), 'wb')
    pickle.dump(res, file)


def loaddata(jobid, subjobid, simid):
    data = open(
        '%s/%s/%s.pkl' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid), '{0:04}'.format(simid)),
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
    gen = countsubjobs(jobid)
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


#################################################

def countsubjobs(jobid):
    count = 0
    for root, dirs, files in os.walk('%s' % ('{0:04}'.format(jobid))):
        for d in dirs:
            count += 1
    return count


########################## ANALYSIS FUNCTIONS ##########################

def mse(res):
    base = loaddata(9999, 0, 0)
    mse_a = np.mean(((res.aco - base.aco) ** 2))
    mse_p = np.mean(((res.pco - base.pco) ** 2))
    score = np.mean([mse_a, mse_p])
    res.scores['mse'] = score
