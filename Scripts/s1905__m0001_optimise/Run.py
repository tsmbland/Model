import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
import M as x
import Models.m0001 as m
import numpy as np
import pickle
import copy
import pandas as pd

"""


"""

pdirec = x.ddirec + 's1905__m0001_optimise/'
params = ['konA', 'kAP', 'kPA']

group = int(sys.argv[1])
dlist = range(group * 10, (group + 1) * 10)
db = pd.read_csv(pdirec + 'db.csv')


def func(i):
    # Import base parameter set
    p = pickle.load(open('%s/_base.pkl' % pdirec, 'rb'))
    model = m.Model(**dict(copy.deepcopy(p)))

    # Set parameters
    d = dict(db.loc[i])
    for p in params:
        setattr(model, p, d[p])

    # Run model
    model.run()

    # Calculate MSE
    mse = np.mean([np.mean((p['am_0'] - model.am) ** 2), np.mean((p['pm_0'] - model.pm) ** 2)])
    np.savetxt(pdirec + '/mse_%s.txt' % '{:05d}'.format(i), [mse])


for d in dlist:
    print(d)
    func(d)
