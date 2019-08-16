import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')

import M as x
import numpy as np
import pandas as pd
import pickle

direc = x.ddirec + 's1904__m0000_optimise_B'
direcs = x.direcslist(direc, 1)
scores = np.zeros([len(direcs)])
params = ['konA', 'konP', 'kAP', 'kPA']
db = []
for d in direcs:
    score = np.loadtxt(d + '/mse.txt')
    p = pickle.load(open('%s/_params.pkl' % d, 'rb'))

    dic = {'Directory': d, 'MSE': score}
    for param in params:
        dic[param] = p[param]
    db.append(dic)

pd.DataFrame(db).to_csv(direc + '/db.csv')
