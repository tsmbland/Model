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

sns.set()
sns.set_style("ticks")

datadirec = None

"""
From local to server: '../../../../../../../Volumes/lab-goehringn/working/Tom/ModelData'

From server to server: '../working/Tom/ModelData'

From local to local: '../../ModelData'
"""


########################## MODEL INITIATION ######################


def init_single(p, direc, name=0):
    """
    Takes parameters class, saves file in new directory specified by direc and name

    :param p: parameters class
    :param direc: parent directory
    :param name: simulation folder
    :return:
    """

    # Make folder
    if not os.path.exists('%s/%s' % (direc, name)):
        os.makedirs('%s/%s' % (direc, name))

    # Create csv file
    with open('%s/%s/_params.csv' % (direc, name), 'w') as f:
        w = csv.DictWriter(f, vars(p).keys())
        w.writeheader()
        w.writerow(vars(p))


def init_multi(p, direc, params, vals):

    # Find all unique parameter combinations
    valsarray = list(itertools.product(*vals))

    # For each simulation (one per parameter combination)
    for j, s in enumerate(valsarray):

        # Create copy of original parameter class
        pcopy = copy.deepcopy(p)

        # Make changes
        for param in range(len(params)):
            setattr(pcopy, params[param], s[param])

        # Save parameters file
        init_single(pcopy, direc, j)


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


def init_rand(p, direc, params, ranges, seeds, nsims):

    # Create array of unique parameter combinations
    valsarray = params_array_rand(params, ranges, seeds, nsims)

    # For each simulation (one per parameter combination)
    for j, s in enumerate(valsarray):

        # Create copy of original parameter class
        pcopy = copy.deepcopy(p)

        # Make changes
        for param in range(len(params)):
            setattr(pcopy, params[param], s[param])

        # Save parameters file
        init_single(pcopy, direc, j)


#########################  ALGORITHM FUNCTIONS ####################


def import_params(direc, model):
    a = pd.read_csv('%s/_params.csv' % direc)
    a = a.to_dict('records')
    return model.Params(**a[0])


def clean_directory(direc):
    flist = glob.glob('%s/*' % direc)
    p = glob.glob('%s/_params.csv' % direc)[0]
    flist.remove(p)
    for f in flist:
        os.remove(f)


def loaddata(direc, m, t):
    # Load parameters
    p = import_params(direc, m)

    # Import time csv, find row
    x = pd.read_csv('%s/time.csv' % direc, header=None)
    row = x[x.values == t].index

    # Create data class
    s = m.Species(p)

    # Loop through species, fill in class
    for key in vars(s):
        x = pd.read_csv('%s/%s.csv' % (direc, key), header=None)
        setattr(s, key, np.array(x.iloc[row])[0])

    return s


############################ ALGORITHMS ##########################


def alg_singlesim(m, direc):
    clean_directory(direc)
    p = import_params(direc, m)
    m = m.Model(p)
    m.run(direc)


def alg_multisim(m, direc):
    for d in glob.glob('%s/*' % direc):
        alg_singlesim(m, d)


def alg_parsim(m, direc, cores):
    Parallel(n_jobs=cores, verbose=51)(
        delayed(alg_singlesim)(m, d) for d in glob.glob('%s/*' % direc))


def alg_parsim_clust(m, direc, cores, node):
    subjobmin = node * cores
    subjobmax = (node + 1) * cores
    Parallel(n_jobs=cores)(
        delayed(alg_singlesim)(m, d) for d in glob.glob('%s/*' % direc)[subjobmin: subjobmax])


########################## ANALYSIS FUNCTIONS ##########################


def mse_a_0(s):
    base = loaddata(...)
    return {'mse_a_0': np.mean(((s.aco - base.aco) ** 2))}


def mse_p_0(s):
    base = loaddata(...)
    return {'mse_p_p': np.mean(((s.pco - base.pco) ** 2))}


def mse_a_1(s):
    base = loaddata(...)
    return {'mse_a_1': np.mean(((s.aco - base.aco) ** 2))}


def mse_p_1(s):
    base = loaddata(...)
    return {'mse_p_1': np.mean(((s.pco - base.pco) ** 2))}


def asi_a(s):
    ant = np.mean(s.aco[len(s.aco) // 2])
    post = np.mean(s.aco[len(s.aco) // 2:])
    return {'asi_a': (ant - post) / (2 * (ant + post))}


def asi_p(s):
    ant = np.mean(s.pco[len(s.pco) // 2])
    post = np.mean(s.pco[len(s.pco) // 2:])
    return {'asi_p': (ant - post) / (2 * (ant + post))}


def domainsize_a(s):
    return {'domainsize_a': np.sum(s.aco > (0.2 * max(s.aco))) / len(s.aco)}


def domainsize_p(s):
    return {'domainsize_p': np.sum(s.pco > (0.2 * max(s.pco))) / len(s.pco)}


all_analysis = [mse_a_0, mse_p_0, mse_a_1, mse_p_1, asi_a, asi_p, domainsize_a, domainsize_p]
a2 = [asi_a, asi_p, domainsize_a, domainsize_p]


def batch_analysis(direc, m, funcs):
    # Remove existing file
    if os.path.exists('%s/_analysis.csv' % direc):
        os.remove('%s/_analysis.csv' % direc)

    # Loop through time points
    times = pd.read_csv('%s/time.csv' % direc, header=None).values
    for time in times:

        # Load data
        s = loaddata(direc, m, time)

        # Create row
        a = {'time': time}
        for f in funcs:
            a.update(f(s))
        row = pd.DataFrame(a)

        # Add to file
        try:
            b = pd.read_csv('%s/_analysis.csv' % direc)
            b = b.append(row)
            b.to_csv('%s/_analysis.csv' % direc, index=False)
        except:
            row.to_csv('%s/_analysis.csv' % direc, index=False)


############################## PLOTS ##############################


def plot_singlesim(direc, m, t):
    d = loaddata(direc, m, t)
    p = import_params(direc, m)
    plt.plot(np.array(range(p.xsteps)) / (p.xsteps - 1) * p.L, d.aco, label='A', c='r')
    plt.plot(np.array(range(p.xsteps)) / (p.xsteps - 1) * p.L, d.pco, label='P', c='dodgerblue')
    plt.show()


def sliderplot(direc, m):
    p = import_params(direc, m)
    times = pd.read_csv('%s/time.csv' % direc, header=None)[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Time (s)', 0, max(times), valinit=0, valfmt='%d')

    def update_slider(val):
        ax.clear()
        t = min(times, key=lambda x: abs(x - int(sframe.val)))
        d = loaddata(direc, m, t)
        ax.plot(np.array(range(p.xsteps)) / (p.xsteps - 1) * p.L, d.aco, label='A', c='r')
        ax.plot(np.array(range(p.xsteps)) / (p.xsteps - 1) * p.L, d.pco, label='P', c='dodgerblue')

    update_slider(0)
    sframe.on_changed(update_slider)
    plt.show()
