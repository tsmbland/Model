import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
import Models.FlowsReview as m
import copy
from joblib import Parallel, delayed
import numpy as np

"""
A: no exchange, vary D

"""

base_model = m.Model(Dm=0, Dc=0, Vm=0.05, Vc=0, kon=0, koff=0, xsteps=100, Tmax=1000, deltat=0.001, deltax=0.1, c_0=0,
                     m_0=1, flowtype=3)


def func(mod, d):
    name = 'Res3/%s_%s' % ('A', d)
    print(name)

    # Set up model
    model = copy.deepcopy(mod)
    model.Dm = d

    # Simulate
    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat

    # Save
    if not os.path.exists(name):
        os.mkdir(name)
    model.save(name + '/')
    print(name + ' Done!')


D_m = [0.01, 0.05, 1]
Parallel(n_jobs=3)(delayed(func)(base_model, D_m[i]) for i in range(len(D_m)))

"""
B: exchanging, vary kon/koff

"""

base_model = m.Model(Dm=0, Dc=1, Vm=0.05, Vc=0, kon=0, koff=0, xsteps=100, Tmax=1000, deltat=0.001, deltax=0.1,
                     c_0=1, m_0=1, flowtype=3)


def func(mod, k):
    name = 'Res3/%s_%s' % ('B', k)
    print(name)

    # Set up model
    model = copy.deepcopy(mod)
    model.koff = k
    model.kon = k

    # Simulate
    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat

    # Save
    if not os.path.exists(name):
        os.mkdir(name)
    model.save(name + '/')
    print(name + ' Done!')


k = [0.007, 0.01, 0.1]
Parallel(n_jobs=3)(delayed(func)(base_model, k[i]) for i in range(len(k)))

"""
C: retrograde flow

"""

base_model = m.Model(Dm=0.01, Dc=1, Vm=0.05, Vc=0, kon=np.r_[0 * np.ones([80]), 5 * np.ones([20])], koff=0,
                     xsteps=100, Tmax=1000, deltat=0.001, deltax=0.1, c_0=1, m_0=1, flowtype=3)


def func(mod, koff):
    name = 'Res3/%s_%s' % ('C', koff)
    print(name)

    # Set up model
    model = copy.deepcopy(mod)
    model.koff = koff
    model.kon *= koff

    # Simulate
    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat

    # Save
    if not os.path.exists(name):
        os.mkdir(name)
    model.save(name + '/')
    print(name + ' Done!')


koff = [0.001, 0.01, 0.1]
Parallel(n_jobs=7)(delayed(func)(base_model, koff[i]) for i in range(len(koff)))
