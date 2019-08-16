import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
import M as x
import Models.m0008 as m
import numpy as np
import pickle
import copy
import pandas as pd

# pdirec = x.ddirec + 's1904__m0008_optimise/'
#
# print(sys.argv)
# group = int(sys.argv[1])
# dlist = x.direcslist(pdirec)[group * 150: (group + 1) * 150]
#
#
# def func(direc):
#     # Import parameters
#     p = pickle.load(open('%s/_params.pkl' % direc, 'rb'))
#
#     # Run model
#     model = m.Model(**dict(copy.deepcopy(p)))
#     model.run()
#
#     # Calculate MSE
#     pm = model.pm1 + model.pm2s + model.pm2d
#     mse = np.mean([np.mean((p['am_0'] - model.am) ** 2), np.mean((p['pm1_0'] - pm) ** 2)])
#     np.savetxt(direc + '/mse.txt', [mse])
#
#
# for d in dlist:
#     func(d)

"""


"""

pdirec = x.ddirec + 's1904__m0008_optimise/'
params = ['kon_a', 'ra', 'kon_p', 'kon_p_2', 'kd_f', 'kd_b', 'rp']

group = int(sys.argv[1])
dlist = range(group * 100, (group + 1) * 100)
db = pd.read_csv('db.csv')


def func(i):
    # Import base parameter set
    p = pickle.load(open('%s/_base.pkl' % pdirec, 'rb'))
    model = m.Model(**dict(copy.deepcopy(p)))

    # Set parameters
    d = dict(db.loc[i])
    for p in params:
        setattr(model, p, d[p])

    """
    Simulation A: Full model
    
    """

    # Run model
    model.run()

    # Calculate MSE
    pm = model.pm1 + model.pm2s + model.pm2d
    mse = np.mean((p['pm1_0'] - pm) ** 2)
    np.savetxt(pdirec + '/mse_%s_A.txt' % '{:05d}'.format(i), [mse])

    """
    Simulation B: PAR-1 RNAi
    
    """

    # Set parameters

    # Run model
    model.run()

    # Calculate MSE
    pm = model.pm1 + model.pm2s + model.pm2d
    mse = np.mean((p['pm1_0'] - pm) ** 2)
    np.savetxt(pdirec + '/mse_%s_B.txt' % '{:05d}'.format(i), [mse])

    """
    Simulation C: PAR-3 -/-
    
    """
    # Set parameters

    # Run model
    model.run()

    # Calculate MSE
    pm = model.pm1 + model.pm2s + model.pm2d
    mse = np.mean((p['pm1_0'] - pm) ** 2)
    np.savetxt(pdirec + '/mse_%s_C.txt' % '{:05d}'.format(i), [mse])

    """
    Simulation D: CRT90
    
    """
    # Set parameters

    # Run model
    model.run()

    # Calculate MSE
    pm = model.pm1 + model.pm2s + model.pm2d
    mse = np.mean((p['pm1_0'] - pm) ** 2)
    np.savetxt(pdirec + '/mse_%s_D.txt' % '{:05d}'.format(i), [mse])


for d in dlist:
    func(d)
