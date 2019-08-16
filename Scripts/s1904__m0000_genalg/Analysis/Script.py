import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import M as x
import numpy as np

direc = x.ddirec + 's1904__m0000_optimise/g001'

popsize = len(x.direcslist(direc))

scores = np.zeros([popsize])
for p in range(popsize):
    scores[p] = np.loadtxt(direc + '/%s/mse.txt' % int(p))
print(np.argsort(scores))
print(np.sort(scores))

