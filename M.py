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
import csv
from pyDOE import lhs

"""


"""

# ddirec = '/Users/blandt/Desktop/ModelData/'


ddirec = os.path.dirname(os.path.realpath(__file__)) + '/../working/Tom/ModelData/'


# ddirec = os.path.dirname(
#     os.path.realpath(__file__)) + '/../../../../../../../Volumes/lab-goehringn/working/Tom/ModelData/'


####################### FUNCTIONS ####################


def init_single(p, direc, name=0):
    """
    Takes parameters dict, saves file in new directory specified by direc and name

    :param p: parameters class
    :param direc: parent directory
    :param name: simulation folder
    :return:
    """

    # Make folder
    if not os.path.exists('%s/%s' % (direc, name)):
        os.makedirs('%s/%s' % (direc, name))

    pickle.dump(p, open('%s/%s/_params.pkl' % (direc, name), 'wb'))


def init_multi(p, direc, params, vals):
    # Find all unique parameter combinations
    valsarray = list(itertools.product(*vals))

    # For each simulation (one per parameter combination)
    for j, s in enumerate(valsarray):

        # Create copy of original parameter dict
        pcopy = copy.deepcopy(p)

        # Make changes
        for param in range(len(params)):
            pcopy[params[param]] = s[param]

        # Save parameters file
        init_single(pcopy, direc, '{:02d}'.format(j))


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
            valsarray[simid, param] = 10 ** random.uniform(np.log10(ranges[param][0]), np.log10(ranges[param][1]))

    return valsarray


def params_array_rand_lhs(params, ranges, nsims):
    """
    Like above but uses latin hypercube sampling

    :param params:
    :param ranges:
    :param nsims:
    :return:
    """

    rands = lhs(len(params), nsims)

    # Create empty array for parameter values
    valsarray = np.zeros([nsims, len(params)])

    # Fill array
    for param in range(len(params)):
        lb = np.log10(ranges[param][0])
        ub = np.log10(ranges[param][1])

        for simid in range(nsims):
            valsarray[simid, param] = 10 ** (lb + (ub - lb) * rands[simid, param])

    return valsarray


def init_rand(p, direc, params, ranges, seeds, nsims):
    """
    :param p: dictionary
    :param direc: parent directory
    :param params: free parameters
    :param ranges:
    :param seeds:
    :param nsims:
    :return:
    """

    # Create array of unique parameter combinations
    valsarray = params_array_rand_lhs(params, ranges, nsims)

    # For each simulation (one per parameter combination)
    for j, s in enumerate(valsarray):

        # Create copy of original parameter dict
        pcopy = copy.deepcopy(p)

        # Make changes
        for param in range(len(params)):
            pcopy[params[param]] = s[param]

        # Save parameters file
        init_single(pcopy, direc, j)


def direcslist(dest, levels=0):
    """

    Gives a list of directories in a given directory (full path)
    Excludes directories that contain !

    :param dest:
    :return:
    """
    lis = glob.glob('%s/*/' % dest)
    for level in range(levels):
        newlis = []
        for e in lis:
            newlis.extend(glob.glob('%s/*/' % e))
        lis = newlis
    lis = [x[:-1] for x in lis if '!' not in x]
    return lis


def bounded_mean_1d(array, bounds, weights=None):
    """
    Averages 1D array over region specified by bounds

    Should add interpolation step first

    Array and weights should be same length

    :param array:
    :param bounds:
    :return:
    """

    if weights is None:
        weights = np.ones([len(array)])
    if bounds[0] < bounds[1]:
        mean = np.average(array[int(len(array) * bounds[0]): int(len(array) * bounds[1] + 1)],
                          weights=weights[int(len(array) * bounds[0]): int(len(array) * bounds[1] + 1)])
    else:
        mean = np.average(np.hstack((array[:int(len(array) * bounds[1] + 1)], array[int(len(array) * bounds[0]):])),
                          weights=np.hstack(
                              (weights[:int(len(array) * bounds[1] + 1)], weights[int(len(array) * bounds[0]):])))
    return mean


def asi(signals):
    ant = bounded_mean_1d(signals, (0, 0.5))
    post = bounded_mean_1d(signals, (0.5, 1))
    return (ant - post) / (2 * (ant + post))
