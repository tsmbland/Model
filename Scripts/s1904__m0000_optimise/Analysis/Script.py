import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import M as x
import numpy as np
import pandas as pd

direc = x.ddirec + 's1904__m0000_optimise_B'
direcs = x.direcslist(direc, 1)
scores = np.zeros([len(direcs)])

for i, d in enumerate(direcs):
    print(d)
    try:
        scores[i] = np.loadtxt(d + '/mse.txt')
    except:
        pass
scores[scores == 0] = np.nan

print(np.array(direcs)[scores.argsort()])
print(np.sort(scores))






