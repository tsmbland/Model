import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
import Models.FlowsReview as m
import copy
from joblib import Parallel, delayed

base_model = m.Model(Dm=0, Dc=10, Vm=0, Vc=0, kon=0, koff=0, xsteps=100, Tmax=1000, deltat=0.001, deltax=0.1, c_0=1)


def func(mod, d, koff):
    name = 'Res/%s_%s' % (d, koff)
    print(name)

    # Set up model
    model = copy.deepcopy(mod)
    model.Dm = d
    model.koff = koff
    model.kon = 2 * koff

    # Simulate without flow
    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat

    # Simulate with flow
    model.time = 0
    # model.Vm = 0.003  # realistic flow
    model.Vm = 0.01  # linear flow

    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat

    # Save
    if not os.path.exists(name):
        os.mkdir(name)
    model.save(name + '/')
    print(name + ' Done!')


"""
Diffusion

"""

d_params = [0.1, 1, 10]
Parallel(n_jobs=8)(delayed(func)(base_model, d_params[i], 0.001) for i in range(len(d_params)))

"""
Off rates

"""

koff_params = [0.0001, 0.001, 0.01, 0.1]
Parallel(n_jobs=8)(delayed(func)(base_model, 1, koff_params[i]) for i in range(len(koff_params)))
