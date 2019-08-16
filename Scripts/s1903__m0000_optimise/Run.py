import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
import M as x
import Models.m0000 as m
import numpy as np
import pickle
import copy


def func(direc):
    # Import parameters
    p = pickle.load(open('%s/_params.pkl' % direc, 'rb'))

    # Run model
    model = m.Model(**dict(copy.deepcopy(p)))
    model.run()

    # Calculate MSE
    mse = np.mean([np.mean((p['am_0'] - model.am) ** 2), np.mean((p['pm_0'] - model.pm) ** 2)])
    np.savetxt(direc + '/mse.txt', [mse])


dlist = x.direcslist(x.ddirec + 's1903__m0000_optimise')
func(dlist[int(sys.argv[1])])
