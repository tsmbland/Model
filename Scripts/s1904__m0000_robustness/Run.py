import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
import M as x
import Models.m0000 as m
import numpy as np
import pickle
import copy

pdirec = x.ddirec + 's1904__m0000_robustness/'

print(sys.argv)
group = int(sys.argv[1])
dlist = x.direcslist(pdirec)[group * 150: (group + 1) * 150]


def func(direc):
    # Import parameters
    p = pickle.load(open('%s/_params.pkl' % direc, 'rb'))

    # Run model
    model = m.Model(**dict(copy.deepcopy(p)))
    model.run()

    # Calculate score
    score = sum(model.am > model.pm)
    if score == 500:  # A dominant
        score = 1
    elif score == 0:  # P dominant
        score = - 1
    else:  # Polarised
        score = 0
    np.savetxt(direc + '/score.txt', [score])


for d in dlist:
    func(d)
